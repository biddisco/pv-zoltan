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
}
//----------------------------------------------------------------------------
vtkParticlePartitionFilter::~vtkParticlePartitionFilter()
{
}
//----------------------------------------------------------------------------
void vtkParticlePartitionFilter::InitializeGhostFlags(vtkPointSet *input)
{
  vtkSmartPointer<vtkUnsignedCharArray> GhostArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  vtkIdType N = input->GetNumberOfPoints();
  GhostArray->SetName("vtkGhostLevels");
  GhostArray->SetNumberOfComponents(1);
  GhostArray->SetNumberOfTuples(N);
  unsigned char *ghost = GhostArray->GetPointer(0);
  for (vtkIdType i=0; i<N; i++) {
    ghost[i]  = 0;
  }
  input->GetPointData()->AddArray(GhostArray);
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
  this->FindPointsInHaloRegions(this->ZoltanCallbackData.Input->GetPoints(), this->MigrateLists.known, this->LoadBalanceData);

  //
  // Based on the original partition and our extra ghost allocations
  // perform the main point exchange between all processes
  //
  this->ManualPointMigrate(this->MigrateLists, false, this->KeepInversePointLists==1);
  

  // add the ghost points to the original points
/*
  vtkDebugMacro(<<"Entering BuildCellToProcessList");
  if (this->ZoltanCallbackData.PointType==VTK_FLOAT) {
    this->BuildCellToProcessList<float>(this->ZoltanCallbackData.Input, 
      cell_partitioninfo,       // lists of which cells to send to which process
      this->MigrateLists.known, // list of which points to send to which process
      this->LoadBalanceData     // the partition information generated during PartitionPoints
    );
  }
  else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
    this->BuildCellToProcessList<double>(this->ZoltanCallbackData.Input, 
      cell_partitioninfo,       // lists of which cells to send to which process
      this->MigrateLists.known, // list of which points to send to which process
      this->LoadBalanceData     // the partition information generated during PartitionPoints
    );
  }

*/




  //
  // clean up arrays that zoltan passed back to us
  //
  if (this->LoadBalanceData.importGlobalGids) {
    Zoltan_LB_Free_Part(&this->LoadBalanceData.importGlobalGids, &this->LoadBalanceData.importLocalGids, &this->LoadBalanceData.importProcs, &this->LoadBalanceData.importToPart);
    Zoltan_LB_Free_Part(&this->LoadBalanceData.exportGlobalGids, &this->LoadBalanceData.exportLocalGids, &this->LoadBalanceData.exportProcs, &this->LoadBalanceData.exportToPart);
  }
  
  //
  // Find points in halo regions and sent them to remote processes
  //
//  if (this->UpdateNumPieces>1) {
//    this->ExchangeHaloPoints(info, inputVector, outputVector);
//  }

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
  //*****************************************************************
  // Free the arrays allocated by Zoltan_LB_Partition, and free
  // the storage allocated for the Zoltan structure.
  //*****************************************************************
  //
  Zoltan_Destroy(&this->ZoltanData);

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
  if (this->InputExtentTranslator && this->InputExtentTranslator->GetBoundsHalosPresent()) {
    this->ExtentTranslator->SetBoundsHalosPresent(1);
    for (int p=0; p<this->UpdateNumPieces; p++) {
      vtkBoundingBox box;  
      box.SetBounds(this->InputExtentTranslator->GetBoundsHaloForPiece(p));
      this->BoxListWithHalo.push_back(box);
      this->ExtentTranslator->SetBoundsHaloForPiece(p,this->InputExtentTranslator->GetBoundsHaloForPiece(p));
    }
  }
  else {
    this->ExtentTranslator->SetBoundsHalosPresent(1);
    // @todo : extend this to handle AMR ghost regions etc.
    std::vector<double> ghostOverlaps(this->UpdateNumPieces,this->GhostCellOverlap);
    for (int p=0; p<this->UpdateNumPieces; p++) {
      vtkBoundingBox box = this->BoxList[p];  
      box.Inflate(ghostOverlaps[p]);
      this->BoxListWithHalo.push_back(box);
      this->ExtentTranslator->SetBoundsHaloForPiece(p, box);
    }
  }
}
//-------------------------------------------------------------------------
void vtkParticlePartitionFilter::FindPointsInHaloRegions(
  vtkPoints *pts, PartitionInfo &point_partitioninfo, ZoltanLoadBalanceData &loadBalanceData)
{
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
  // 2) loop over all to be exported and note the destination
  for (vtkIdType i=0; i<loadBalanceData.numExport; i++) {
    vtkIdType id               = loadBalanceData.exportGlobalGids[i] - this->ZoltanCallbackData.ProcessOffsetsPointId[this->ZoltanCallbackData.ProcessRank];
    localId_to_process_map[id] = loadBalanceData.exportProcs[i];
  }

  //
  // This vector will store {Id,process} pairs wilst the list is being scanned
  typedef std::pair<vtkIdType, int> process_tuple;
  std::vector<process_tuple> process_vector;

  // Since we already have a list of points to export, we don't want to
  // duplicate them, so traverse the points list once per process
  // skipping those already flagged for export
  vtkIdType N = pts->GetNumberOfPoints(), pE=0;
  for (int proc=0; proc<this->UpdateNumPieces; proc++) {
    vtkBoundingBox &b = this->BoxListWithHalo[proc];
    int pc = 0;
    //
    // any remote process box (+halo) which does not overlap our local box can be ignored
    //
    if (localPointsbox.Intersects(b)) {      
      for (vtkIdType i=0; i<N; i++) {
        bool marked = false;
        vtkIdType gID = i + this->ZoltanCallbackData.ProcessOffsetsPointId[this->ZoltanCallbackData.ProcessRank];
        // if this ID is already marked as exported to the process then we don't need to send it again
        // But, if it's marked for export and we need a local copy, we must add it to our keep list
        if (/*pE<loadBalanceData.numExport && */loadBalanceData.exportGlobalGids[pE]==gID && loadBalanceData.exportProcs[pE]==proc) {
          pE++;
          continue;
        }
        double *pt = pts->GetPoint(i);
        if (b.ContainsPoint(pt)) {
          if (proc==this->UpdatePiece) {
            if (localId_to_process_map[i]==this->UpdatePiece) {
            }
            else {
              point_partitioninfo.LocalIdsToKeep.push_back(i);
            }
          }
          else {
            process_vector.push_back( process_tuple(gID, proc) );
            pc++;
          }
        }
      }
      std::cout << "Rank " << this->UpdatePiece << " exporting " << pc << " ghost particles to rank " << proc << std::endl;
      std::cout << "Rank " << this->UpdatePiece << " LocalIdsToKeep " << point_partitioninfo.LocalIdsToKeep.size() << std::endl;
    }
  }

  //
  // After examining ghosts, we found some points we need to send away,
  // we must add these to the points the original load balance already declared in export lists 
  //
  point_partitioninfo.GlobalIds.reserve(loadBalanceData.numExport + process_vector.size());
  point_partitioninfo.Procs.reserve(loadBalanceData.numExport + process_vector.size());
  
  // 1) add the points from original zoltan load balance to send list
  for (int i=0; i<loadBalanceData.numExport; ++i) {
    point_partitioninfo.GlobalIds.push_back(loadBalanceData.exportGlobalGids[i]);
    point_partitioninfo.Procs.push_back(loadBalanceData.exportProcs[i]);
  }

  // 2) add the points from ghost tests just performed to send list
  for (std::vector<process_tuple>::iterator x=process_vector.begin(); x!=process_vector.end(); ++x) 
  {
    point_partitioninfo.GlobalIds.push_back(x->first);
    point_partitioninfo.Procs.push_back(x->second);
  }

  vtkDebugMacro(<<"FindPointsInHaloRegions "  << 
    " numImport : " << this->LoadBalanceData.numImport <<
    " numExport : " << point_partitioninfo.GlobalIds .size()
  );

}
//----------------------------------------------------------------------------
int vtkParticlePartitionFilter::ExchangeHaloPoints(vtkInformation* info,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  this->ExtentTranslator->InitWholeBounds();
  this->LocalBox     = &this->BoxList[this->UpdatePiece];
  this->LocalBoxHalo = &this->BoxListWithHalo[this->UpdatePiece];

  //
  // Find points which overlap other processes' ghost regions
  // note that we must use the 'new' migrated points which are not the same
  // as the original input points (might be bigger/smaller), so get the new IdArray 
  //
  //vtkIdTypeArray *newIds = vtkIdTypeArray::SafeDownCast(
  //  this->ZoltanCallbackData.Output->GetPointData()->GetArray(this->IdsName.c_str()));
  //if (!newIds || newIds->GetNumberOfTuples()!=this->ZoltanCallbackData.Output->GetPoints()->GetNumberOfPoints()) {
  //  vtkErrorMacro(<<"Fatal : Ids on migrated data corrupted");
  //  return 0;
  //}

  //PartitionInfo GhostIds;
//  this->FindPointsInHaloRegions(this->ZoltanCallbackData.Input->GetPoints(), GhostIds);

  int num_found = 0; //this->ManualPointMigrate(GhostIds, true, this->KeepInversePointLists);

  //
  // Ghost information : Paraview doesn't let us visualize an array called vtkGhostLevels
  // because it's an 'internal' array, so we make an extra one for debug purposes
  //
  vtkSmartPointer<vtkUnsignedCharArray> GhostArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  vtkSmartPointer<vtkIntArray> GhostArray2 = vtkSmartPointer<vtkIntArray>::New();
  vtkIdType N = this->ZoltanCallbackData.Output->GetNumberOfPoints();
  GhostArray->SetName("vtkGhostLevels");
  GhostArray->SetNumberOfComponents(1);
  GhostArray->SetNumberOfTuples(N);
  GhostArray2->SetName("GhostLevels");
  GhostArray2->SetNumberOfComponents(1);
  GhostArray2->SetNumberOfTuples(N);
  unsigned char *ghost = GhostArray->GetPointer(0);
  int          *ghost2 = GhostArray2->GetPointer(0);
  for (vtkIdType i=0; i<N; i++) {
    if (i<(N-num_found)) {
      ghost[i]  = 0;
      ghost2[i] = 0;
    }
    else {
      ghost[i]  = 1;
      ghost2[i] = 1;
    }
  }
  this->ZoltanCallbackData.Output->GetPointData()->AddArray(GhostArray);
  this->ZoltanCallbackData.Output->GetPointData()->AddArray(GhostArray2);

  return 1;
}
//----------------------------------------------------------------------------
