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
void vtkParticlePartitionFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
//----------------------------------------------------------------------------
void vtkParticlePartitionFilter::SetupGlobalIds(vtkPointSet *ps) 
{
  /*
  //
  // Global Ids : always do them after other point arrays are setup 
  //
  if (this->IdChannelArray) {
    this->IdsName = this->IdChannelArray;
  }
  if (this->IdsName.empty() || this->IdsName==std::string("Not available")) {
    this->IdsName = "ZPF_PointIds";
  } 
*/
  vtkSmartPointer<vtkPointData> pd = ps->GetPointData();
  vtkSmartPointer<vtkIdTypeArray> Ids = NULL;
  //
  Ids = vtkIdTypeArray::SafeDownCast(pd->GetArray(this->IdsName.c_str()));
  if (!Ids) {
    // Try loading the user supplied global ids.
    Ids = vtkIdTypeArray::SafeDownCast(pd->GetGlobalIds());
  }
  if (!Ids) {
    // and increment the callbackdata field count
    this->ZoltanCallbackData.NumberOfFields++;
  }
  // Generate our own if none exist
  vtkDebugMacro(<<"About to Init Global Ids");
  vtkIdType numPoints = ps->GetNumberOfPoints();
  vtkIdType  numCells = ps->GetNumberOfCells();
  this->ComputeIdOffsets(numPoints, numCells);

  //
  // Global point IDs generated here
  //
  vtkIdType offset = this->ZoltanCallbackData.ProcessOffsetsPointId[this->UpdatePiece];
  //
  if (!Ids) {
    Ids = vtkSmartPointer<vtkIdTypeArray>::New();
    Ids->SetNumberOfValues(numPoints);
    for (vtkIdType id=0; id<numPoints; id++) {
      Ids->SetValue(id, id + offset);
    }
    Ids->SetName(this->IdsName.c_str());
    vtkDebugMacro(<< "Generated Ids with " << numPoints << " values");
  }

  ps->GetPointData()->AddArray(Ids);
  vtkDebugMacro(<<"Global Ids Initialized");
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
  
  //
  // Initialize the halo regions based on our redistribution
  //
  this->FillPartitionBoundingBoxWithHalo();

  //
  // Find points in halo regions and sent them to remote processes
  //
  if (this->UpdateNumPieces>1) {
    this->ExchangeHaloPoints(info, inputVector, outputVector);
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
void vtkParticlePartitionFilter::FillPartitionBoundingBoxWithHalo()
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
#if 1 
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
    for (std::vector<boxproc>::iterator it=good_boxes.begin(); it!=good_boxes.end(); ++it) 
    {
      vtkBoundingBox *b = (it->first);
      if (b->ContainsPoint(pt)) {
//        ghostinfo.LocalIds.push_back(i);
        ghostinfo.GlobalIds.push_back(IdArray->GetValue(i));
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
  if (!newIds || newIds->GetNumberOfTuples()!=this->ZoltanCallbackData.Output->GetPoints()->GetNumberOfPoints()) {
    vtkErrorMacro(<<"Fatal : Ids on migrated data corrupted");
    return 0;
  }

  PartitionInfo GhostIds;
  this->FindPointsInHaloRegions(this->ZoltanCallbackData.Output->GetPoints(), newIds, GhostIds);

  int num_found = this->ManualPointMigrate(GhostIds, true);

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
