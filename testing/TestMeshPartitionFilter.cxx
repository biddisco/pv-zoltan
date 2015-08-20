// mpishim : C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\Remote Debugger\x64\mpishim100.exe
// mpiexec : C:\Program Files\MPICH2\bin\mpiexec.exe 
// mpiargs : -localonly -n 2 -env PATH C:\cmakebuild\pv-meshless\bin\debug;c:\bin

#ifdef _WIN32
  #include <windows.h>
#else 
  #include <sys/time.h>
#endif

#define _USE_MATH_DEFINES
#include <math.h>
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
#include "vtkActor.h"
#include "vtkAppendPolyData.h"
#include "vtkCamera.h"
#include "vtkPointSource.h"
#include "vtkDataSet.h"
#include "vtkMath.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleSwitch.h"
#include "vtkRenderer.h"
#include "vtkWindowToImageFilter.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkDebugLeaks.h"
#include "vtkProperty.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkTimerLog.h"
#include "vtkBoundingBox.h"
#include "vtkOutlineSource.h"
#include "vtkProcessIdScalars.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPPolyDataReader.h"
#include "vtkTransform.h"
//
#include <vtksys/SystemTools.hxx>
#include <sstream>
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <algorithm>
//
#include "TestUtils.h"
//
#include "vtkMeshPartitionFilter.h"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int main (int argc, char* argv[])
{
  int retVal = 1;
  char *empty = "";
  bool ok = true;

  //--------------------------------------------------------------
  // Setup Test Params
  //--------------------------------------------------------------
  TestStruct test;
  initTest(argc, argv, test);

  // if testing partition from file
  double read_elapsed = 0.0;
  double partition_elapsed = 0.0;
  vtkSmartPointer<vtkAlgorithm> data_algorithm; 

  test.CreateXMLPolyDataReader();
  test.xmlreader->Update();
  test.ghostLevels = 0;

  //--------------------------------------------------------------
  // Parallel partition
  //--------------------------------------------------------------
  test.CreatePartitioner_Mesh();
  test.partitioner->SetInputConnection(test.xmlreader->GetOutputPort());
  test.partitioner->SetInputDisposable(1);
  test.partitioner->SetKeepInversePointLists(1);
  // setup ghost options
  static_cast<vtkMeshPartitionFilter*>(test.partitioner.GetPointer())->SetGhostMode(test.ghostMode);
  if (test.ghostMode==vtkMeshPartitionFilter::Neighbour) {
      static_cast<vtkMeshPartitionFilter*>(test.partitioner.GetPointer())->SetNumberOfGhostLevels(1);
  }
  else if (test.ghostMode==vtkMeshPartitionFilter::BoundingBox) {
      static_cast<vtkMeshPartitionFilter*>(test.partitioner.GetPointer())->SetGhostCellOverlap(test.ghostOverlap);
  }
  //
  partition_elapsed = test.UpdatePartitioner();

  //--------------------------------------------------------------
  // Add process Id's
  //--------------------------------------------------------------
  vtkSmartPointer<vtkProcessIdScalars> processId = vtkSmartPointer<vtkProcessIdScalars>::New();
  processId->SetInputConnection(test.partitioner->GetOutputPort());
  processId->SetController(test.controller);
  //
//  test.controller->Barrier();
  if (test.myRank==0) {
    testDebugMacro( "Process Id : " << test.myRank << " Generated N Points : " << test.generateN );
  }

  //--------------------------------------------------------------
  // Update in parallel
  //
  // To get parallel operation correct, we need to make sure that piece
  // information is passed upstream. first update information,
  // then set piece update extent,
  //--------------------------------------------------------------
  testDebugMacro( "Setting piece information " << test.myRank << " of " << test.numProcs );
  vtkStreamingDemandDrivenPipeline *sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(processId->GetExecutive());
  // no piece info set yet, assumes info is not piece dependent
  sddp->UpdateInformation();
  // now set piece info and update
  sddp->SetUpdateExtent(0, test.myRank, test.numProcs, 0);
  sddp->Update();
  testDebugMacro("Update completed . "<<test.myRank);

  if (test.doRender) {
    //
    // Send all the data to process zero for display
    //
    vtkSmartPointer<vtkPolyData> OutputData;
    OutputData.TakeReference(vtkPolyData::SafeDownCast(sddp->GetOutputData(0)->NewInstance()));
    OutputData->ShallowCopy(sddp->GetOutputData(0));

    if (test.myRank>0) {
      testDebugMacro("data sending from "<<test.myRank);
      test.controller->Send(OutputData, 0, DATA_SEND_TAG);
    }
    //
    // Rank 0 collect all data pieces from parallel processes
    //
    else if (test.myRank==0) {
        retVal = test.RenderPieces(argc, argv, OutputData);
    }
  }

  if (ok && test.myRank==0) {
//    DisplayParameter<vtkIdType>("Total Particles", "", &totalParticles, 1, test.myRank);
    DisplayParameter<double>("Read Time", "", &read_elapsed, 1, test.myRank);
    DisplayParameter<double>("Partition Time", "", &partition_elapsed, 1, test.myRank);
    DisplayParameter<char *>("====================", "", &empty, 1, test.myRank);
  }

  // manually free partitioner so Zoltan structures are freed before MPI finalize

  processId->SetInputConnection(NULL);
  processId = NULL;
  //
  test.controller->Barrier();
  //
  test.DeleteXMLPolyDataReader();
  test.DeletePartitioner();
  //
  test.controller->Finalize();
  //
  return !retVal;
}
//----------------------------------------------------------------------------
