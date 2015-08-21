/*=========================================================================

  Module                  : vtkParticlePartitionFilter.cxx

  Copyright (C) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing
  1) This copyright notice appears on all copies of source code
  2) An acknowledgment appears with any substantial usage of the code
  3) If this code is contributed to any other open source project, it
  must not be reformatted such that the indentation, bracketing or
  overall style is modified significantly.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
//
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetAttributes.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"
#include "vtkIdTypeArray.h"
#include "vtkBoundingBox.h"
#include "vtkMath.h"
#include "vtkPointLocator.h"
#include "vtkPKdTree.h"
//
// For PARAVIEW_USE_MPI
#include "vtkPVConfig.h"
#ifdef PARAVIEW_USE_MPI
  #include "vtkMPI.h"
  #include "vtkMPIController.h"
  #include "vtkMPICommunicator.h"
#endif
#include "vtkDummyController.h"
//
#include "vtkBoundsExtentTranslator.h"
#include "vtkParticlePartitionFilter.h"
//
#include <sstream>
//
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <numeric>
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkParticlePartitionFilter);
//----------------------------------------------------------------------------
// vtkParticlePartitionFilter :: implementation 
//----------------------------------------------------------------------------
vtkParticlePartitionFilter::vtkParticlePartitionFilter()
{
  this->GhostCellOverlap = 0.0;
  this->NumberOfGhostLevels= 0;
  this->GridSpacing      = 0.0;
  this->GridOrigin[0]    = this->GridOrigin[1] = this->GridOrigin[2] = 0.0;
}
//----------------------------------------------------------------------------
vtkParticlePartitionFilter::~vtkParticlePartitionFilter()
{
}
//----------------------------------------------------------------------------
int vtkParticlePartitionFilter::RequestData(vtkInformation* info,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  //
  // Calculate even distribution of points across processes
  // This step only performs the load balance analysis, 
  // no actual sending of data takes place yet.
  //
    this->PartitionPoints(info, inputVector, outputVector);

  if (this->UpdateNumPieces==1) {
    // input has been copied to output during PartitionPoints
    return 1;
  }

  //
  // Initialize the halo regions around the computed bounding boxes of the distribution
  //
  this->AddHaloToBoundingBoxes();

  //
  // based on the point partition, decide which particles are ghost particles
  // and need to be sent/kept in addition to the default load balance
  //
  // Concept: One algorithm would be to 
  // a) compute load balance (done in PartitionPoints above)
  // b) ascertain ghost particles based on the expected load balance
  // c) perform a particle exchange of all a+b) together in one go
  // however when attempting this, it was found during step b) some particles are being sent as
  // non-ghost particles to one process and ghost particles others - this means that the c) exchange must
  // send different ghost flags to different processes -  we can't just send/receive a normal list
  // but instead have to maintain a map of process/ghost/ID flags. This makes the send/receive
  // more expensive as the maps must be checked for every particle sent (to each process). 
  //
  // A second algorithm is to
  // a) compute load balance (done in PartitionPoints above)
  // b) ascertain ghost particles based on the expected load balance
  // c) compute the inverse lists for ghost exchange
  // d) perform a particle exchange for just the main load balance step from a), but allocating
  //    all space determined by a+c) so that only one final list is allocated.
  //    Note that some points marked for transfer will stay on this process as ghosts
  //    which does require extra bookkeeping.
  // e) perform a particle exchange using the lists from c), during the exchange we know that
  //    all send/received are ghosts so ghost flags can be set unilaterally
  //
  MigrationLists ghost_info;
  this->FindPointsInHaloRegions(this->ZoltanCallbackData.Input->GetPoints(), this->MigrateLists.known, this->LoadBalanceData, ghost_info.known);

  //
  // create the inverse map of who sends/receives from who : Ghost particles
  //
  this->ComputeInvertLists(ghost_info);
  std::cout << "Rank " << this->UpdatePiece << " Imports " << ghost_info.num_found << " ghost particles " << std::endl;

  //
  // create the inverse map of who sends/receives from who : Core particles
  //
  this->ComputeInvertLists(this->MigrateLists);
  std::cout << "Rank " << this->UpdatePiece << " Imports " << this->MigrateLists.num_found << " core particles " << std::endl;

  // we want to reserve some extra space for the ghost particles when they are sent in
  this->MigrateLists.num_reserved = ghost_info.num_found;
  
  //
  // Based on the original load balance step perform the point exchange for core particles
  // pass in ghost info so that space can be allocated for the final 
  //
  // NB : not yet supporting retaining of migration lists for later point data exchange
  this->ManualPointMigrate(this->MigrateLists, false); 

  // we have now allocated the output and filled the point data for non ghost Ids
  vtkIdType N = this->ZoltanCallbackData.Output->GetNumberOfPoints();

  // create a ghost array
  vtkSmartPointer<vtkUnsignedCharArray> GhostArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  GhostArray->SetName("vtkGhostLevels");
  GhostArray->SetNumberOfComponents(1);
  GhostArray->SetNumberOfTuples(N);
  unsigned char *ghost = GhostArray->GetPointer(0);

  vtkIdType i=0;
  for (; i<(N-ghost_info.num_found); i++) {
    ghost[i] = 0;
  }
  
  for (; i<N; i++) {
    ghost[i] = 1;
  }
  
//  for (i=0; i<ghost_info.known.GlobalIds.size(); i++) {
//    ghost[ghost_info.known.GlobalIds[i] - this->ZoltanCallbackData.ProcessOffsetsPointId[this->ZoltanCallbackData.ProcessRank]] = 1;
//  }
  
  // some local points were kept as they were inside the local ghost region, we need to mark them
  for (std::vector<vtkIdType>::iterator it = this->MigrateLists.known.LocalIdsToKeep.begin(); it!=this->MigrateLists.known.LocalIdsToKeep.end(); ++it) {
    vtkIdType Id = this->ZoltanCallbackData.LocalToLocalIdMap[*it];
    ghost[Id] = 2;
  }
  
  // now exchange ghost cells too
  this->ZoltanPointMigrate(ghost_info, false);


  //
  // clean up arrays that zoltan passed back to us
  //
  if (this->LoadBalanceData.importGlobalGids) {
//    Zoltan_LB_Free_Part(&this->LoadBalanceData.importGlobalGids, &this->LoadBalanceData.importLocalGids, &this->LoadBalanceData.importProcs, &this->LoadBalanceData.importToPart);
//    Zoltan_LB_Free_Part(&this->LoadBalanceData.exportGlobalGids, &this->LoadBalanceData.exportLocalGids, &this->LoadBalanceData.exportProcs, &this->LoadBalanceData.exportToPart);
  }
  
  //
  // If polydata create Vertices for each final point
  //
  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  if (vtkPolyData::SafeDownCast(this->ZoltanCallbackData.Output)) {
    vtkIdType N = this->ZoltanCallbackData.Output->GetNumberOfPoints();
    vtkIdType *arraydata = cells->WritePointer(N, 2*N);
    for (int i=0; i<this->ZoltanCallbackData.Output->GetNumberOfPoints(); i++) {
      arraydata[i*2]   = 1;
      arraydata[i*2+1] = i;
    }
    vtkPolyData::SafeDownCast(this->ZoltanCallbackData.Output)->SetVerts(cells);
  }

  //
  // build a tree of bounding boxes to use for rendering info/hints or other spatial tests
  //
//  this->CreatePkdTree();
//  this->ExtentTranslator->SetKdTree(this->GetKdtree());

  //
  //*****************************************************************
  // Free the arrays allocated by Zoltan_LB_Partition, and free
  // the storage allocated for the Zoltan structure.
  //*****************************************************************
  //
  Zoltan_Destroy(&this->ZoltanData);

  this->ZoltanCallbackData.Output->GetPointData()->AddArray(GhostArray);

  this->Controller->Barrier();
  this->Timer->StopTimer();
  vtkDebugMacro(<<"Particle partitioning : " << this->Timer->GetElapsedTime() << " seconds");
  return 1;
}
//----------------------------------------------------------------------------
vtkBoundingBox *vtkParticlePartitionFilter::GetPartitionBoundingBoxWithHalo(int partition)
{
  if (partition<this->BoxListWithHalo.size()) {
    return &this->BoxListWithHalo[partition];
  }
  vtkErrorMacro(<<"Partition not found in Bounding Box list");
  return NULL;
}
//----------------------------------------------------------------------------
void vtkParticlePartitionFilter::AddHaloToBoundingBoxes()
{
  //
  // Set the halo/ghost regions we need around each process bounding box
  //  
  this->BoxListWithHalo.clear();
  if (this->InputExtentTranslator && this->InputExtentTranslator->GetBoundsHalosEnabled()) {
    this->ExtentTranslator->SetBoundsHalosEnabled(1);
    for (int p=0; p<this->UpdateNumPieces; p++) {
      vtkBoundingBox box;  
      box.SetBounds(this->InputExtentTranslator->GetBoundsHaloForPiece(p));
      this->BoxListWithHalo.push_back(box);
      this->ExtentTranslator->SetBoundsHaloForPiece(p,this->InputExtentTranslator->GetBoundsHaloForPiece(p));
    }
  }
  else {
    this->ExtentTranslator->SetBoundsHalosEnabled(1);
    // @todo : extend this to handle AMR ghost regions etc.
    std::vector<double> ghostOverlaps(this->UpdateNumPieces,this->GhostCellOverlap);
    for (int p=0; p<this->UpdateNumPieces; p++) {
      vtkBoundingBox box = this->BoxList[p];      box.Inflate(ghostOverlaps[p]*(this->NumberOfGhostLevels+1));
      this->BoxListWithHalo.push_back(box);
      this->ExtentTranslator->SetBoundsHaloForPiece(p, box);
    }
  }
}
//-------------------------------------------------------------------------
void vtkParticlePartitionFilter::FindPointsInHaloRegions(
  vtkPoints *pts, PartitionInfo &point_partitioninfo, ZoltanLoadBalanceData &loadBalanceData, PartitionInfo &ghost_info)
{
  // create N maps to store ghost assignments
  std::vector< std::map<int, vtkIdType> > GhostProcessMap(this->UpdateNumPieces);
  //
  // What is the bounding box of points we have been given to start with
  //
  vtkIdType numPts = pts->GetNumberOfPoints();
  double bounds[6];
  pts->GetBounds(bounds);
  vtkBoundingBox localPointsbox(bounds);

  // we know that some points on this process will be exported to remote processes
  // so build a point to process map to quickly lookup the process Id from the point Id
  // 1) initially set all points as local to this process
  std::vector<int> localId_to_process_map(numPts, this->UpdatePiece); 
  std::vector<int> ghost_flag(numPts, 0); 
  // 2) loop over all to be exported and note the destination
  int offset = this->ZoltanCallbackData.ProcessOffsetsPointId[this->ZoltanCallbackData.ProcessRank];
  for (vtkIdType i=0; i<loadBalanceData.numExport; i++) {
    vtkIdType id               = loadBalanceData.exportGlobalGids[i] - offset;
    // cout<<"i: "<<i<<" \tlocal_id:"<<id<<"\tproc:"<<loadBalanceData.exportProcs[i]<<">>"<<this->UpdatePiece<<endl;
    localId_to_process_map[id] = loadBalanceData.exportProcs[i];
    // points marked for export by load balance are not ghost points
    ghost_flag[i] = 0;
  }
  // Since we already have a list of points to export, we don't want to
  // duplicate them, so traverse the points list once per process
  // skipping those already flagged for export
  vtkIdType N = pts->GetNumberOfPoints(), pE=0;
  for (int proc=0; proc<this->UpdateNumPieces; proc++) {
    vtkBoundingBox &b = this->BoxListWithHalo[proc];
    int pc = 0;
    //
    // any remote process box (+halo) which does not overlap our local points does not need to be tested
    //
    if (localPointsbox.Intersects(b)) {      
      for (vtkIdType i=0; i<N; i++) {
        vtkIdType gID = i + this->ZoltanCallbackData.ProcessOffsetsPointId[this->ZoltanCallbackData.ProcessRank];
        // if this ID is already marked as exported to the process then we don't need to send it again
        // But, if it's marked for export and we need a local copy, we must add it to our keep list
        if (/*pE<loadBalanceData.numExport && */loadBalanceData.exportGlobalGids[pE]==gID && loadBalanceData.exportProcs[pE]==proc) {
          pE++;
          continue;
        }
        double *pt = pts->GetPoint(i);
        if (b.ContainsPoint(pt)) {
          // if the bounding box is actually our local box
          if (proc==this->UpdatePiece) {
            if (localId_to_process_map[i]==this->UpdatePiece) {
              // this point is due to stay on this process anyway
            }
            else {
              // this point has already been flagged for export but we need it as a ghost locally
              point_partitioninfo.LocalIdsToKeep.push_back(i);
//              ghost_flag[i] = 1;
            }
          }
          // the bounding box is a remote one
          else {
            if (localId_to_process_map[i]==proc) {
              // this point is already marked for export to the process so it is not a ghost cell
//              ghost_flag[i] = 0;
//              point_partitioninfo.LocalIdsToSend.push_back(i);
            }
            else {
              // this point is due to be exported to one process as a non ghost 
              // but another copy must be sent to a different process as a ghost
              ghost_info.GlobalIds.push_back(gID);
              ghost_info.Procs.push_back(proc);
//              point_partitioninfo.LocalIdsToKeep.push_back(i);
//              point_partitioninfo.LocalIdsToSend.push_back(i);
//              ghost_flag[i] = 1;
//              GhostProcessMap[proc][i] = 1;

            }
            pc++;
          }
        }
      }
      std::cout << "Rank " << this->UpdatePiece << " exporting " << pc << " ghost particles to rank " << proc << std::endl;
    }
  }

  for (int i=0; i<this->UpdateNumPieces; i++) {
    std::cout << "Rank " << this->UpdatePiece << " exporting " << GhostProcessMap[i].size() << " ghost particles to rank " << i << std::endl;
  }
  std::cout << "Rank " << this->UpdatePiece << " LocalIdsToKeep " << point_partitioninfo.LocalIdsToKeep.size() << std::endl;
  
  //
  // Set the send list to the points from original zoltan load balance 
  //
  point_partitioninfo.nIDs         = loadBalanceData.numExport;
  point_partitioninfo.GlobalIdsPtr = loadBalanceData.exportGlobalGids;
  point_partitioninfo.ProcsPtr     = loadBalanceData.exportProcs;
  
  vtkDebugMacro(<<"FindPointsInHaloRegions "  << 
    " numImport : " << this->LoadBalanceData.numImport <<
    " numExport : " << point_partitioninfo.GlobalIds .size()
  );

}
//----------------------------------------------------------------------------
