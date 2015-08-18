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
#define DATA_SEND_TAG 301
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
  vtkIdType totalParticles = 0;

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
  static_cast<vtkMeshPartitionFilter*>(test.partitioner.GetPointer())->SetNumberOfGhostLevels(test.ghostLevels);
//  test.partitioner->SetGhostCellOverlap(test.ghostOverlap); // no delta required in mesh partition
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
    
//    std::stringstream ss;
//    #define temp_debug(a) \
//    ss<<a;
//
//    for (vtkIdType i=0; i<OutputData->GetNumberOfCells(); i++) {
//      vtkIdType *pts;
//      vtkIdType npts;
//      OutputData->GetCellPoints(i, npts, pts);
//      temp_debug("cell: "<<i);
//      for (int j = 0; j<npts; j++){
//        temp_debug(" ,"<<pts[j]);
//      }
//      temp_debug(" >>> "<<test.myRank<<std::endl);
//    }
//    temp_debug(" >>> "<<test.myRank<<std::endl);
//    
//    for (vtkIdType i=0; i<OutputData->GetNumberOfPoints(); i++) {
//      vtkIdType *pts;
//      unsigned short npts;
//      OutputData->GetPointCells(i, npts, pts);
//      temp_debug("point: "<<i);
//      for (int j=0; j<npts; j++) {
//        temp_debug(" ,"<<pts[j]);
//      }
//      temp_debug(" >>> "<<test.myRank<<std::endl);
//    }
    
    if (test.myRank>0) {
      testDebugMacro("data sending from "<<test.myRank);
      test.controller->Send(OutputData, 0, DATA_SEND_TAG);
    }
    //
    // Rank 0 collect all data pieces from parallel processes
    //
    else if (test.myRank==0) {
      //
      vtkSmartPointer<vtkRenderer>                ren = vtkSmartPointer<vtkRenderer>::New();
      vtkSmartPointer<vtkRenderWindow>      renWindow = vtkSmartPointer<vtkRenderWindow>::New();
      vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
      vtkSmartPointer<vtkInteractorStyleSwitch> style = vtkSmartPointer<vtkInteractorStyleSwitch>::New();
      iren->SetRenderWindow(renWindow);
      iren->SetInteractorStyle(style);
      style->SetCurrentStyleToTrackballCamera();
      ren->SetBackground(0.1, 0.1, 0.1);
      renWindow->SetSize(test.windowSize);
      renWindow->AddRenderer(ren);
      //
      // To make display of ghost cells and boundary regions better, find the centre of
      // all the pieces and use that to apply a transform to the actors to shift them
      // away from the centre so the edges don't touch and overlapping boundary cells are visible
      //
      double centre[3], midpoint[3]={0,0,0};
      for (int i=0; i<test.numProcs; i++) {
          vtkBoundingBox *box = test.partitioner->GetPartitionBoundingBox(i);
          box->GetCenter(centre);
          for (int d=0; d<3; d++) midpoint[d] += centre[d]/test.numProcs;
      }
      //
      for (int i=0; i<test.numProcs; i++) {
        vtkSmartPointer<vtkPolyData> pd;
        if (i==0) {
          pd = OutputData;
        }
        else {
          pd = vtkSmartPointer<vtkPolyData>::New();
          testDebugMacro("data receiving at "<<test.myRank);
          test.controller->Receive(pd, i, DATA_SEND_TAG);
          testDebugMacro("data received at "<<test.myRank<<" from "<< i);
        }
        pd->PrintSelf(std::cout, vtkIndent(0));
        vtkSmartPointer<vtkPolyDataMapper>       mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSmartPointer<vtkActor>                 actor = vtkSmartPointer<vtkActor>::New();
        mapper->SetInputData(pd);
        mapper->SetImmediateModeRendering(1);
        mapper->SetColorModeToMapScalars();
        if (test.scalarmode==0) {
            mapper->SetScalarModeToUsePointFieldData();
            mapper->SetUseLookupTableScalarRange(0);
            mapper->SetScalarRange(0,test.numProcs-1);
            mapper->SetInterpolateScalarsBeforeMapping(0);
        }
        else {
            mapper->SetScalarModeToUseCellFieldData();
            mapper->SetUseLookupTableScalarRange(0);
            mapper->SetInterpolateScalarsBeforeMapping(0);
            mapper->SetScalarRange(0,6);
        }
        mapper->SelectColorArray(test.scalarname.c_str());
        actor->SetMapper(mapper);
        actor->GetProperty()->SetPointSize(2);
        ren->AddActor(actor);
        // move each actor away from the midpoint so we can see ghost cells better
        vtkBoundingBox *box = test.partitioner->GetPartitionBoundingBox(i);
        box->GetCenter(centre);
        vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        transform->PostMultiply();
        transform->Translate(-0.25*(midpoint[0]-centre[0]), -0.25*(midpoint[1]-centre[1]), -0.25*(midpoint[2]-centre[2]));
        actor->SetUserTransform(transform);
        //
        if (test.cameraSet) {
          ren->GetActiveCamera()->SetPosition(test.cameraPosition);
          ren->GetActiveCamera()->SetFocalPoint(test.cameraFocus);
          ren->GetActiveCamera()->SetViewUp(test.cameraViewUp);
          ren->ResetCameraClippingRange();
        }
        else {
          ren->ResetCamera();
        }
      }
      //
      // Display boxes for each partition
      //
      for (int i=0; i<test.numProcs; i++) {
        vtkBoundingBox *box = test.partitioner->GetPartitionBoundingBox(i);
        double bounds[6];
        box->GetBounds(bounds);
        vtkSmartPointer<vtkOutlineSource> boxsource = vtkSmartPointer<vtkOutlineSource>::New();
        boxsource->SetBounds(bounds);
        vtkSmartPointer<vtkPolyDataMapper> bmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSmartPointer<vtkActor>          bactor = vtkSmartPointer<vtkActor>::New();
        bmapper->SetInputConnection(boxsource->GetOutputPort());
        bactor->SetMapper(bmapper);
        ren->AddActor(bactor);
        // move each box away from the midpoint so we can see ghost cells better
        box->GetCenter(centre);
        vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        transform->PostMultiply();
        transform->Translate(-0.25*(midpoint[0]-centre[0]), -0.25*(midpoint[1]-centre[1]), -0.25*(midpoint[2]-centre[2]));
        bactor->SetUserTransform(transform);
      }
      
      testDebugMacro( "Process Id : " << test.myRank << " About to Render" );
      renWindow->Render();

      retVal = vtkRegressionTester::Test(argc, argv, renWindow, 10);

      if ( retVal == vtkRegressionTester::DO_INTERACTOR) {
        iren->Start();
      }
      ok = (retVal==vtkRegressionTester::PASSED);
      testDebugMacro( "Process Id : " << test.myRank << " Rendered " << (ok?"Pass":"Fail"));
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
