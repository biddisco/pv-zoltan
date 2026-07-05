// Test for vtkPartitionOutline's ShowKdTreeBounds mode.
// Generates a random particle dataset, partitions it across MPI ranks,
// and then renders the k-d tree produced by vtkParticlePartitionFilter.

#ifdef _WIN32
  #include <windows.h>
#else
  #include <sys/time.h>
#endif

#define _USE_MATH_DEFINES
#include <math.h>
//
// For PARAVIEW_USE_MPI
#include "vtkPVVersion.h"
#ifdef PARAVIEW_USE_MPI
  #include "vtkMPI.h"
  #include "vtkMPIController.h"
  #include "vtkMPICommunicator.h"
#endif
#include "vtkDummyController.h"
//
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkDataSet.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkWindowToImageFilter.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkProperty.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkTimerLog.h"
#include "vtkBoundingBox.h"
#include "vtkOutlineSource.h"
#include "vtkTransform.h"
//
#include <vtksys/SystemTools.hxx>
#include <sstream>
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
//
#include "TestUtils.h"
//
#include "vtkParticlePartitionFilter.h"
#include "vtkPartitionOutline.h"

//----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int retVal = 1;
  bool ok = true;

  //--------------------------------------------------------------
  // Setup Test Params
  //--------------------------------------------------------------
  TestStruct test;
  initTest(argc, argv, test);

  double read_elapsed = 0.0;
  double partition_elapsed = 0.0;
  vtkIdType totalParticles = 0;

  //--------------------------------------------------------------
  // allocate points + arrays
  //--------------------------------------------------------------
  vtkSmartPointer<vtkPolyData>  Sprites = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints>     points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray>   verts = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkIdTypeArray>   Ids = vtkSmartPointer<vtkIdTypeArray>::New();
  vtkSmartPointer<vtkIntArray>    Ranks = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkFloatArray> Weights = vtkSmartPointer<vtkFloatArray>::New();

  points->SetNumberOfPoints(test.generateN);
  verts->AllocateExact(test.generateN, test.generateN);
  Sprites->SetPoints(points);
  Sprites->SetVerts(verts);

  Ids->SetNumberOfTuples(test.generateN);
  Ids->SetNumberOfComponents(1);
  Ids->SetName("PointIds");
  Sprites->GetPointData()->AddArray(Ids);

  Ranks->SetNumberOfTuples(test.generateN);
  Ranks->SetNumberOfComponents(1);
  Ranks->SetName("Rank");
  Sprites->GetPointData()->AddArray(Ranks);

  Weights->SetNumberOfTuples(test.generateN);
  Weights->SetNumberOfComponents(1);
  Weights->SetName("Weights");
  Sprites->GetPointData()->AddArray(Weights);

  //--------------------------------------------------------------
  // Create default scalar arrays
  //--------------------------------------------------------------
  double radius = 500.0;
  if (test.particleGenerator == 0) {
    SpherePoints(test.generateN, radius * (1.5 + test.myRank) / (test.numProcs + 0.5),
                 vtkFloatArray::SafeDownCast(points->GetData())->GetPointer(0));
  }
  else {
    CubePoints(test.generateN, radius,
               vtkFloatArray::SafeDownCast(points->GetData())->GetPointer(0),
               Weights->GetPointer(0));
  }
  for (vtkIdType Id = 0; Id < test.generateN; Id++) {
    Ids->SetTuple1(Id, Id + test.myRank * test.generateN);
    Ranks->SetTuple1(Id, test.myRank);
    verts->InsertNextCell(1, &Id);
  }

  //--------------------------------------------------------------
  // Parallel partition
  //--------------------------------------------------------------
  test.CreatePartitioner_Particles();
  test.partitioner->SetInputData(Sprites);
  if (test.useWeights) {
    test.partitioner->SetPointWeightsArrayName("Weights");
  }
  static_cast<vtkParticlePartitionFilter*>(test.partitioner.GetPointer())
      ->SetGhostHaloSize(test.ghostOverlap);
  partition_elapsed = test.UpdatePartitioner();

  //--------------------------------------------------------------
  // Generate the partition outline from the k-d tree
  //--------------------------------------------------------------
  vtkSmartPointer<vtkPartitionOutline> outline = vtkSmartPointer<vtkPartitionOutline>::New();
  outline->SetInputConnection(test.partitioner->GetOutputPort());
  outline->ShowKdTreeBoundsOn();

  vtkStreamingDemandDrivenPipeline* sddp =
      vtkStreamingDemandDrivenPipeline::SafeDownCast(outline->GetExecutive());
  sddp->UpdateInformation();
  {
    vtkInformation* outInfo = sddp->GetOutputInformation(0);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), test.myRank);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), test.numProcs);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), 0);
  }
  sddp->Update();

  test.controller->Barrier();
  if (test.myRank == 0) {
    testDebugMacro("Process Id : " << test.myRank << " Generated N Points : " << test.generateN);
  }

  if (test.doRender) {
    // Send all the outline data to process zero for display
    vtkSmartPointer<vtkPolyData> OutputData;
    OutputData.TakeReference(vtkPolyData::SafeDownCast(sddp->GetOutputData(0)->NewInstance()));
    OutputData->ShallowCopy(sddp->GetOutputData(0));
    if (test.myRank > 0) {
      test.controller->Send(OutputData, 0, DATA_SEND_TAG);
    }
    else if (test.myRank == 0) {
      retVal = test.RenderPieces(argc, argv, OutputData);
    }
  }

  if (ok && test.myRank == 0) {
    DisplayParameter<vtkIdType>("Total Particles", "", &totalParticles, 1, test.myRank);
    DisplayParameter<double>("Read Time", "", &read_elapsed, 1, test.myRank);
    DisplayParameter<double>("Partition Time", "", &partition_elapsed, 1, test.myRank);
  }

  FinalizeKokkos();
  test.controller->Finalize();

  return !retVal;
}
