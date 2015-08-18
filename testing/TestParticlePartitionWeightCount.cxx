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
#include "vtkParticlePartitionFilter.h"

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
    char *empty = " ";
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

    //--------------------------------------------------------------
    // allocate points + arrays
    //--------------------------------------------------------------
    vtkSmartPointer<vtkPolyData>  Sprites = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints>     points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray>   verts = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkIdTypeArray>   Ids = vtkSmartPointer<vtkIdTypeArray>::New();
    vtkSmartPointer<vtkIntArray>    Ranks = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkFloatArray> Weights = vtkSmartPointer<vtkFloatArray>::New();
    //
    points->SetNumberOfPoints(test.generateN);
    //
    verts->Allocate(test.generateN,test.generateN);
    Sprites->SetPoints(points);
    Sprites->SetVerts(verts);
    //
    Ids->SetNumberOfTuples(test.generateN);
    Ids->SetNumberOfComponents(1);
    Ids->SetName("PointIds");
    Sprites->GetPointData()->AddArray(Ids);
    //
    Ranks->SetNumberOfTuples(test.generateN);
    Ranks->SetNumberOfComponents(1);
    Ranks->SetName("Rank");
    Sprites->GetPointData()->AddArray(Ranks);
    //
    Weights->SetNumberOfTuples(test.generateN);
    Weights->SetNumberOfComponents(1);
    Weights->SetName("Weights");
    Sprites->GetPointData()->AddArray(Weights);
    //
    //--------------------------------------------------------------
    // Create default scalar arrays
    //--------------------------------------------------------------
    double radius  = 500.0;
    test.ghostOverlap = radius*0.1; // ghost_region
    test.ghostLevels = 0;

    known_seed();
    if (test.particleGenerator==0) {
        SpherePoints(test.generateN, radius*(1.5+test.myRank)/(test.numProcs+0.5), vtkFloatArray::SafeDownCast(points->GetData())->GetPointer(0));
    }
    else {
        CubePoints(test.generateN, radius,
                vtkFloatArray::SafeDownCast(points->GetData())->GetPointer(0), Weights->GetPointer(0));
    }
    for (vtkIdType Id=0; Id<test.generateN; Id++) {
        Ids->SetTuple1(Id, Id + test.myRank*test.generateN);
        Ranks->SetTuple1(Id, test.myRank);
        verts->InsertNextCell(1,&Id);
    }
    /*
  // Randomly give some processes zero points to improve test coverage
  random_seed();
  if (test.numProcs>1 && rand()%2==1) {
    test.generateN = 0;
    Sprites = vtkSmartPointer<vtkPolyData>::New();
  }
     */
    //--------------------------------------------------------------
    // Parallel partition
    //--------------------------------------------------------------
    test.CreatePartitioner_Particles();
    test.partitioner->SetInputData(Sprites);
    if (test.useWeights) {
        test.partitioner->SetPointWeightsArrayName("Weights");
    }
    //  test.partitioner->SetIdChannelArray("PointIds");
    static_cast<vtkParticlePartitionFilter*>(test.partitioner.GetPointer())->SetGhostCellOverlap(test.ghostOverlap);
    static_cast<vtkParticlePartitionFilter*>(test.partitioner.GetPointer())->SetNumberOfGhostLevels(test.ghostLevels);
    partition_elapsed = test.UpdatePartitioner();

    //--------------------------------------------------------------
    // Add process Id's
    //--------------------------------------------------------------
    vtkSmartPointer<vtkProcessIdScalars> processId = vtkSmartPointer<vtkProcessIdScalars>::New();
    processId->SetInputConnection(test.partitioner->GetOutputPort());
    processId->SetController(test.controller);
    //
    test.controller->Barrier();
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

    //
    // Send all the data to process zero for display
    //
    vtkSmartPointer<vtkPolyData> OutputData;
    OutputData.TakeReference(vtkPolyData::SafeDownCast(sddp->GetOutputData(0)->NewInstance()));
    OutputData->ShallowCopy(sddp->GetOutputData(0));
    if (test.myRank>0) {
        test.controller->Send(OutputData, 0, DATA_SEND_TAG);
    }
    //
    // Rank 0 collect all data pieces from parallel processes
    //
    else if (test.myRank==0) {
        //
        std::vector<int> pointsCounts(test.numProcs, 0);
        std::vector<double> weightCounts(test.numProcs, 0);
        //
        for (int i=0; i<test.numProcs; i++) {
            vtkSmartPointer<vtkPolyData> pd;
            if (i==0) {
                pd = OutputData;
            }
            else {
                pd = vtkSmartPointer<vtkPolyData>::New();
                testDebugMacro("data receiving at " << test.myRank);
                test.controller->Receive(pd, i, DATA_SEND_TAG);
                testDebugMacro("data received at " << test.myRank << " from " << i);
            }
            //        pd->PrintSelf(std::cout, vtkIndent(0));
            vtkFloatArray *weights = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetArray("Weights"));
            vtkUnsignedCharArray *ghosts = vtkUnsignedCharArray::SafeDownCast(pd->GetPointData()->GetArray("vtkGhostLevels"));

            for (int n = 0; n < weights->GetNumberOfTuples(); ++n)
            {
                // ghost cells are trnasferred in and are were not used for the weighting
                int ghost = ghosts->GetValue(n);
                if (ghost==0) {
                    weightCounts[i] += weights->GetValue(n);
                    pointsCounts[i] ++;
                }
            }
        }
        for (int i=0; i<test.numProcs; i++){
            std::cout << "Rank " << i << " num points : " << pointsCounts[i] << "\t, Weight sum : " << weightCounts[i] << "\n";
        }
        double sum = std::accumulate(weightCounts.begin(), weightCounts.end(), 0.0);
        double mean = sum / weightCounts.size();

        double sq_sum = std::inner_product(weightCounts.begin(), weightCounts.end(), weightCounts.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / weightCounts.size() - mean * mean);
        std::cout << "standard deviation : " << stdev << ")\n";
        ok = (stdev<0.2);
    }

    if (ok && test.myRank==0) {
        DisplayParameter<vtkIdType>("Total Particles", "", &totalParticles, 1, test.myRank);
        DisplayParameter<double>("Read Time", "", &read_elapsed, 1, test.myRank);
        DisplayParameter<double>("Partition Time", "", &partition_elapsed, 1, test.myRank);
        DisplayParameter<char *>("====================", "", &empty, 1, test.myRank);
    }
    retVal = (ok==true);

    test.controller->Finalize();

    return !retVal;
}
//----------------------------------------------------------------------------
