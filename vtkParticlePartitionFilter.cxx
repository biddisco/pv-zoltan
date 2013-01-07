/*=========================================================================

  Project                 : pv-meshless
  Module                  : vtkParticlePartitionFilter.h
  Revision of last commit : $Rev: 884 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2010-04-06 12:03:55 +0200 #$

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
#include "vtkParticlePartitionFilter.h"
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
#include "vtkParticlePartitionFilter.h"
#include "vtkDummyController.h"
//
#include "vtkBoundsExtentTranslator.h"
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
// Function Type: Pre migration callback for halo particle exchange
// the halo migration is controlled manually. We declare what we need to exchange
// and use this function to ensure all memory is allocated correctly for the receive
// to take place and unpack each object
//----------------------------------------------------------------------------
template <typename T>
void vtkParticlePartitionFilter::zoltan_pre_migrate_func_halo(void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  CallbackData *mesh = (CallbackData*)data;
  // resize points to accept ghost cell additions
  mesh->OutputNumberOfPointsWithHalo = mesh->OutputNumberOfLocalPoints + num_import;
  mesh->OutputPoints->GetData()->Resize(mesh->OutputNumberOfPointsWithHalo);
  mesh->OutputPoints->SetNumberOfPoints(mesh->OutputNumberOfPointsWithHalo);
  mesh->OutputPointsData = (T*)(mesh->OutputPoints->GetData()->GetVoidPointer(0));
  // copies are now being made from existing data, so use output as new input points
  mesh->InputPointsData = mesh->OutputPointsData;
  vtkPointData   *inPD  = mesh->Output->GetPointData();
  vtkPointData   *outPD = mesh->Output->GetPointData();
  //
  // we must resize all the scalar/vector fields for the point data
  // WARNING: because this routine takes places after the auto-migration of points during
  // the zoltan load-balance, when copying point data, we set the input pointer to the output pointer 
  // because the output pointer is in fact the 'new' input
  mesh->InputArrayPointers.clear();
  mesh->OutputArrayPointers.clear();
  for (int i=0; i<mesh->NumberOfFields; i++) {
    vtkDataArray *oarray = mesh->Output->GetPointData()->GetArray(i);
    oarray->Resize(mesh->OutputNumberOfPointsWithHalo);
    oarray->SetNumberOfTuples(mesh->OutputNumberOfPointsWithHalo);
    mesh->InputArrayPointers.push_back(oarray->GetVoidPointer(0));
    mesh->OutputArrayPointers.push_back(oarray->GetVoidPointer(0));
  }
}
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
void vtkParticlePartitionFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
//----------------------------------------------------------------------------
int vtkParticlePartitionFilter::RequestData(vtkInformation* info,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  //
  // Distribute points evenly across processes
  //
  this->PartitionPoints(info, inputVector, outputVector);

  //
  // clean up arrays that zoltan passed back to us
  //
  if (this->LoadBalanceData.importGlobalGids) {
    Zoltan_LB_Free_Part(&this->LoadBalanceData.importGlobalGids, &this->LoadBalanceData.importLocalGids, &this->LoadBalanceData.importProcs, &this->LoadBalanceData.importToPart);
    Zoltan_LB_Free_Part(&this->LoadBalanceData.exportGlobalGids, &this->LoadBalanceData.exportLocalGids, &this->LoadBalanceData.exportProcs, &this->LoadBalanceData.exportToPart);
  }
  
  this->ExchangeHaloPoints(info, inputVector, outputVector);


  this->Controller->Barrier();
  this->Timer->StopTimer();
  vtkDebugMacro(<<"Particle partitioning : " << this->Timer->GetElapsedTime() << " seconds");
  return 1;
}
//----------------------------------------------------------------------------
void vtkParticlePartitionFilter::InitBoundingBoxes(vtkDataSet *input, vtkBoundingBox &box) 
{
  this->vtkZoltanV1PartitionFilter::InitBoundingBoxes(input, box);
  // hardly need to bother since on one process we are not going to use the halo regions
  double bounds[6];
  input->GetBounds(bounds);
  this->BoxList.clear();
  this->BoxListWithHalo.clear();
  if (this->UpdateNumPieces==1) {
    vtkBoundingBox databox(bounds);
    this->BoxList.push_back(databox);
    // we add a ghost cell region to our boxes
    databox.Inflate(this->GhostCellOverlap);
    this->BoxListWithHalo.push_back(databox);
    this->ExtentTranslator->SetNumberOfPieces(1);
    // Copy the bounds to our piece to bounds translator
    this->ExtentTranslator->SetBoundsForPiece(0, bounds);
    this->ExtentTranslator->InitWholeBounds();
  } 
}
//-------------------------------------------------------------------------
void vtkParticlePartitionFilter::FindPointsInHaloRegions(vtkPoints *pts, vtkIdTypeArray *IdArray, PartitionInfo &ghostinfo)
{
  typedef std::pair<vtkBoundingBox*,int> boxproc;
  // we don't want to test against boxes that are far away, so first reject
  // any boxes that don't overlap with our local box.
  std::vector< boxproc > good_boxes;
  int proc = 0;
  for (std::vector<vtkBoundingBox>::iterator it=this->BoxListWithHalo.begin(); 
    it!=this->BoxListWithHalo.end(); ++it, proc++) 
  {
    vtkBoundingBox &b = *it;
#if 0 
    if (b!=*(this->LocalBoxHalo) && this->LocalBoxHalo->Intersects(b)) {
#else
    if (b!=*(this->LocalBoxHalo)) {
#endif
      good_boxes.push_back( boxproc(&b,proc));
    }
  }
  // now test all points against those boxes which are valid
  vtkIdType N = pts->GetNumberOfPoints();
  for (vtkIdType i=0; i<N; i++) {
    double *pt = pts->GetPoint(i);
    for (std::vector<boxproc>::iterator it=good_boxes.begin(); 
      it!=good_boxes.end(); ++it) 
    {
      vtkBoundingBox *b = (it->first);
      if (b->ContainsPoint(pt)) {
        ghostinfo.GlobalIds.push_back(IdArray->GetValue(i));
        ghostinfo.LocalIds.push_back(i);
        ghostinfo.Procs.push_back(it->second);
      }
    }
  }
//  std::cout << "There are " << ghostinfo.GlobalIds.size() << " ghosts on the list" <<std::endl;
//  std::ostream_iterator<vtkIdType> out_it (cout,", ");
//  copy ( ghostinfo.GlobalIds.begin(), ghostinfo.GlobalIds.end(), out_it );
}
//----------------------------------------------------------------------------
int vtkParticlePartitionFilter::ExchangeHaloPoints(vtkInformation* info,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  //
  // Set the halo/ghost regions we need around each process bounding box
  //  
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
    // do a calculation on all nodes
    std::vector<double> ghostOverlaps(this->UpdateNumPieces,this->GhostCellOverlap);
    for (int p=0; p<this->UpdateNumPieces; p++) {
      vtkBoundingBox box = this->BoxList[p];  
      box.Inflate(ghostOverlaps[p]);
      this->BoxListWithHalo.push_back(box);
    }
  }

  this->ExtentTranslator->InitWholeBounds();
  this->LocalBox     = &this->BoxList[this->UpdatePiece];
  this->LocalBoxHalo = &this->BoxListWithHalo[this->UpdatePiece];

  //
  // Find points which overlap other processes' ghost regions
  // note that we must use the 'new' migrated points which are not the same
  // as the original input points (might be bigger/smaller), so get the new IdArray 
  //
  vtkIdTypeArray *newIds = vtkIdTypeArray::SafeDownCast(
    this->ZoltanCallbackData.Output->GetPointData()->GetArray(this->IdsName.c_str()));
  if (!newIds || newIds->GetNumberOfTuples()!=this->ZoltanCallbackData.OutputPoints->GetNumberOfPoints()) {
    vtkErrorMacro(<<"Fatal : Ids on migrated data corrupted");
    return 0;
  }

  PartitionInfo GhostIds;
  this->FindPointsInHaloRegions(this->ZoltanCallbackData.OutputPoints, newIds, GhostIds);

  //
  // Pass the lists of ghost cells to zoltan so that it
  // can build a list of lists for exchanges between processes
  //
  size_t        num_known = GhostIds.GlobalIds.size(); 
  int           num_found = 0;
  ZOLTAN_ID_PTR found_global_ids = NULL;
  ZOLTAN_ID_PTR found_local_ids  = NULL;
  int          *found_procs      = NULL;
  int          *found_to_part    = NULL;
  //
  int zoltan_error = Zoltan_Invert_Lists(this->ZoltanData, 
    (int)num_known,
    num_known>0 ? &GhostIds.GlobalIds[0] : NULL,
    num_known>0 ? &GhostIds.LocalIds[0]  : NULL,
    num_known>0 ? &GhostIds.Procs[0]     : NULL,
    num_known>0 ? &GhostIds.Procs[0]     : NULL,
    &num_found,
    &found_global_ids,
    &found_local_ids,
    &found_procs,
    &found_to_part); 

  if (zoltan_error != ZOLTAN_OK){
    printf("Zoltan_LB_Partition NOT OK...\n");
    MPI_Finalize();
    Zoltan_Destroy(&this->ZoltanData);
    exit(0);
  }

  //
  // Before sending, we need to change the pre-migrate function as we are now adding
  // extra ghost cells and not starting our lists from a clean slate.
  //
  if (this->ZoltanCallbackData.PointType==VTK_FLOAT) {
    zprem_fn f4 = zoltan_pre_migrate_func_halo<float>;
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &this->ZoltanCallbackData);
  }
  else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
    zprem_fn f4 = zoltan_pre_migrate_func_halo<double>;
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &this->ZoltanCallbackData);
  }

  //
  // Now we can actually send ghost particles between processes
  //
  zoltan_error = Zoltan_Migrate (this->ZoltanData,
    (int)num_found,
    found_global_ids,
    found_local_ids,
    found_procs,
    found_to_part,
    (int)num_known,
    num_known>0 ? &GhostIds.GlobalIds[0] : NULL,
    num_known>0 ? &GhostIds.LocalIds[0]  : NULL,
    num_known>0 ? &GhostIds.Procs[0]     : NULL,
    num_known>0 ? &GhostIds.Procs[0]     : NULL
    );

  //
  // Release the arrays allocated during Zoltan_Invert_Lists
  //
  Zoltan_LB_Free_Part(&found_global_ids, &found_local_ids, 
    &found_procs, &found_to_part);

  //
  // Ghost information : Paraview doesn't let us visualize an array called vtkGhostLevels
  // because it's an 'internal' array, so we make an extra one for debug purposes
  //
  vtkSmartPointer<vtkUnsignedCharArray> GhostArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  vtkSmartPointer<vtkIntArray> GhostArray2 = vtkSmartPointer<vtkIntArray>::New();
  GhostArray->SetName("vtkGhostLevels");
  GhostArray->SetNumberOfComponents(1);
  GhostArray->SetNumberOfTuples(this->ZoltanCallbackData.OutputNumberOfLocalPoints + num_found);
  GhostArray2->SetName("GhostLevels");
  GhostArray2->SetNumberOfComponents(1);
  GhostArray2->SetNumberOfTuples(this->ZoltanCallbackData.OutputNumberOfLocalPoints + num_found);
  unsigned char *ghost = GhostArray->GetPointer(0);
  int          *ghost2 = GhostArray2->GetPointer(0);
  for (vtkIdType i=0; i<this->ZoltanCallbackData.OutputNumberOfLocalPoints + num_found; i++) {
    if (i<this->ZoltanCallbackData.OutputNumberOfLocalPoints) {
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
