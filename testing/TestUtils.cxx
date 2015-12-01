#define _USE_MATH_DEFINES
#include <cmath>
//
#include "TestUtils.h"
//
#ifdef _WIN32
  #include <windows.h>
  #undef min
  #undef max
#else 
  #include <sys/time.h>
#endif
//
#include <algorithm>
#include <random>
#include <cmath>
#include <vtksys/SystemTools.hxx>
//
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPPolyDataReader.h"
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
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#ifdef _WIN32
void known_seed() {
  srand(1);
}
#else
void known_seed() {
  srandom(12345);
}
#endif
//----------------------------------------------------------------------------
// Template specialization for string types to handle spaces in names
//
template <>
std::string GetParameter(const char *argstr, const char *message, int argc, char **argv, std::string defaultvalue, int rank, bool &valueset)
{
  char *tempChar = vtkTestUtilities::GetArgOrEnvOrDefault(argstr, argc, argv, "", "");
  std::string newValue = defaultvalue;
  valueset = false;
  if (std::string(tempChar).size()) {
    std::stringstream temp(tempChar);
    while (temp.good()) {
        std::string tempval;
        temp >> tempval;
        newValue += tempval;
        if (temp.good()) newValue += " ";
    }
    if (rank==0) {
      DisplayParameter<std::string>(message, "", &newValue, 1, rank);
    }
    valueset = true;
  }
  delete []tempChar;
  return newValue;
}
//----------------------------------------------------------------------------
void SpherePoints(int n, float radius, float X[]) {
  double x, y, z, w, t;
  std::default_random_engine generator(12345);
  std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
  for (int i=0; i<n; i++ ) {
    double r1 = uniform_dist(generator);
    double r2 = uniform_dist(generator);
    z = 2.0 * r1 - 1.0;
    t = 2.0 * M_PI * r2;
    w = radius * sqrt( 1 - z*z );
    x = w * cos( t );
    y = w * sin( t );
    X[3*i+0] = static_cast<float>(x);
    X[3*i+1] = static_cast<float>(y);
    X[3*i+2] = static_cast<float>(z*radius);
  }
}

//----------------------------------------------------------------------------
void CubePoints(int n, float radius, float X[], float W[]) {
  // generate uniformly distributed data to evenly fill a cube
  std::default_random_engine generator(12345);
  std::uniform_real_distribution<double> uniform_dist(0.0, radius);

  for (int i=0; i<n; i++ ) {
    double r1 = uniform_dist(generator);
    double r2 = uniform_dist(generator);
    double r3 = uniform_dist(generator);
    X[3*i+0] = static_cast<float>(r1) - radius/2.0;
    X[3*i+1] = static_cast<float>(r2) - radius/2.0;
    X[3*i+2] = static_cast<float>(r3) - radius/2.0;
    W[i] = r1*r2*r3 / (radius*radius*radius);
  }
}
//----------------------------------------------------------------------------
int initTest(int argc, char* argv[], TestStruct &test)
{
  MPI_Init(&argc,&argv);
  test.controller = vtkSmartPointer<vtkMPIController>::New();
  test.controller->Initialize(&argc, &argv, 1);
  
  bool unused;
  // Obtain the id of the running process and the total number of processes
  test.myRank = test.controller->GetLocalProcessId();
  test.numProcs = test.controller->GetNumberOfProcesses();
  //
  test.gridSpacing[0] = test.gridSpacing[1] = test.gridSpacing[2] = 0.0;
  test.gridResolution[0] = test.gridResolution[1] = test.gridResolution[2] = -1;
  test.vminmax[0] = test.vminmax[1] = 0.0;
  test.vpos[0] = test.vpos[1] = test.vpos[2] = 0.0;
  test.cameraPosition[0] = test.cameraPosition[1] = test.cameraPosition[2] = 0.0;
  test.cameraFocus[0] = test.cameraFocus[1] = test.cameraFocus[2] = 0.0;
  test.cameraViewUp[0] = 0.0;
  test.cameraViewUp[1] = 0.0;
  test.cameraViewUp[2] = 1.0;
  test.windowSize[0] = test.windowSize[1] = 400; // +8;

  test.scalarMode = 0;
  test.actor_shift = 0.0;

  // uncomment this to wait for debugger attach
  // DEBUG_WAIT
  //
  test.controller->Barrier();

  //--------------------------------------------------------------
  // command line params : Setup testing utilities/args etc
  //--------------------------------------------------------------
  vtkSmartPointer<vtkTesting> vtktest = vtkSmartPointer<vtkTesting>::New();
  for (int c=1; c<argc; c++ ) {
    vtktest->AddArgument(argv[c]);
  }
  char *empty = " ";

  //
  // Force the creation of our output window object
  //
//  vtkSmartPointer<vtkStreamOutputWindow> outwin = vtkSmartPointer<vtkStreamOutputWindow>::New();
//  vtkOutputWindow::SetInstance(outwin);
//  outwin->SetOutputStream(&std::cout);

  //
  // General test flags/info
  //
  DisplayParameter<char *>("====================", "", &empty, 1, (test.myRank==0)?0:-1);
  test.testName = GetParameter<std::string>("-testName", "Test name", argc, argv, "", test.myRank, unused);
  test.doRender = GetParameter<bool>("-doRender", "Enable Render", argc, argv, 0, test.myRank, unused);
  test.keepTempFiles = GetParameter<bool>("-X", "Keep Temporary Files", argc, argv, 0, test.myRank, unused);

  //
  // ParticleGenerate info
  //
  test.generateN = GetParameter<vtkIdType>("-generateParticles", "Generated Particles", argc, argv, 0, test.myRank, unused);
  test.particleGenerator = GetParameter<int>("-particleGenerator", "Generator for particles (sphere=0, cube=1)", argc, argv, 0, test.myRank, unused);
  test.useWeights = GetParameter<bool>("-useWeights", "Enable weights in partitioning", argc, argv, 0, test.myRank, unused);

  //
  // File load / H5Part info
  //
  std::string filename = GetParameter<std::string>("-F", "Filename", argc, argv, "", test.myRank, unused);
  std::string filepath = GetParameter<std::string>("-D", "Filepath", argc, argv, "", test.myRank, unused);
  if (filename.size() && filepath.size()) {
    test.fullName = vtksys::SystemTools::ConvertToOutputPath(std::string(filepath+"/"+filename).c_str());
    DisplayParameter<std::string>("FullName", "", &test.fullName, 1, (test.myRank==0)?0:-1);
  }
  test.Xarray = GetParameter<std::string>("-Xarray", "Xarray name", argc, argv, "", test.myRank, unused);
  test.Yarray = GetParameter<std::string>("-Yarray", "Yarray name", argc, argv, "", test.myRank, unused);
  test.Zarray = GetParameter<std::string>("-Zarray", "Zarray name", argc, argv, "", test.myRank, unused);
  test.ignorePartitions = GetParameter<bool>("-ignorePartitions", "Ignore Partitions", argc, argv, 0, test.myRank, unused);
  test.randomizeExtents = GetParameter<bool>("-randomizeExtents", "Randomize Extents", argc, argv, 0, test.myRank, unused);

  //
  // SPH kernel or neighbour info
  //
  test.particleSize = GetParameter<double>("-particlesize", "Particle Size", argc, argv, 0, test.myRank, test.fixRadius);
  test.ghostOverlap = GetParameter<double>("-ghostOverlap", "Ghost Region size", argc, argv, 0.0, test.myRank, unused);
  GetArrayParameter<double>("-gridSpacing", "Grid Spacing", test.gridSpacing, 3, argc, argv, test.myRank);
  GetArrayParameter<int>("-gridResolution", "Grid Resolution", test.gridResolution, 3, argc, argv, test.myRank);
  test.maxN = GetParameter<int>("-neighbours", "Fixed Neighbours", argc, argv, 0, test.myRank, test.fixNeighbours);
  test.massScalars = GetParameter<std::string>("-massScalars", "Mass Scalar Array", argc, argv, "", test.myRank, unused);
  test.densityScalars = GetParameter<std::string>("-densityScalars", "Density Scalar Array", argc, argv, "", test.myRank, unused);
  test.expectedN = GetParameter<vtkIdType>("-expectedparticles", "Simulation Particles", argc, argv, 0, test.myRank, unused);

  //
  // Test/Display of results
  //
  test.scalarName = GetParameter<std::string>("-scalarName", "Testing Scalar Array", argc, argv, "", test.myRank, unused);
  test.scalarMode = GetParameter<bool>("-scalarMode", "Point(0) or Cell(1) data", argc, argv, "", test.myRank, unused);
  test.contourVal = GetParameter<double>("-contour", "Contour Value", argc, argv, 0.0, test.myRank, unused);
  GetArrayParameter<double>("-scalarRange", "Scalar Range", test.scalarRange, 2, argc, argv, test.myRank);
  test.actor_shift = GetParameter<double>("-actorShift", "Shift pieces by amount", argc, argv, 0.0, test.myRank, unused);

  test.imageResample = GetParameter<bool>("-imageResample", "imageResample", argc, argv, 0, test.myRank, unused);
  test.skipImageTest = GetParameter<bool>("-skipImageTest", "skipImageTest", argc, argv, 0, test.myRank, unused);
  test.imageScalars = GetParameter<std::string>("-imageScalars", "Image Scalar Array", argc, argv, "", test.myRank, unused);
  GetArrayParameter<double>("-value_range", "Expected Value Range", test.vminmax, 2, argc, argv, test.myRank);
  GetArrayParameter<double>("-peak_position", "Expected Peak Position", test.vpos, 3, argc, argv, test.myRank);
  test.imageThreshold = GetParameter<int>("-imageThreshold", "Image Threshold Pass/Fail", argc, argv, 1, test.myRank, unused);
  test.benchmarkPartition = GetParameter<bool>("-benchmarkPartition", "benchmarkPartition", argc, argv, 0, test.myRank, unused);

  //
  // Window/Camera
  //
  test.cameraSet = GetArrayParameter<double>("-cameraPosition", "Camera Position", test.cameraPosition, 3, argc, argv, test.myRank);
  unused = GetArrayParameter<double>("-cameraFocus", "Camera Focus", test.cameraFocus, 3, argc, argv, test.myRank);
  unused = GetArrayParameter<double>("-cameraViewUp", "Camera ViewUp", test.cameraViewUp, 3, argc, argv, test.myRank);
  unused = GetArrayParameter<int>("-windowSize", "Window Size", test.windowSize, 2, argc, argv, test.myRank);
//#ifdef WIN32
  if (0 && unused) { // why have window sizes changed?
 //   test.windowSize[0] += 8;
 //   test.windowSize[1] += 8;
  }
//#endif  
  // bug fix for cmd line params on windows with debugger (only first read properly)
  test.gridSpacing[2] = test.gridSpacing[1] = test.gridSpacing[0];
  test.gridResolution[2] = test.gridResolution[1] = test.gridResolution[0];
  //
  DisplayParameter<vtkTypeInt64>("No. of Processes", "", &test.numProcs, 1, (test.myRank==0)?0:-1);
  DisplayParameter<char *>("--------------------", "", &empty, 1, (test.myRank==0)?0:-1);
  //
  return 1;
}  
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void TestStruct::CreateXMLPolyDataReader()
{
  this->xmlreader = vtkSmartPointer<vtkXMLPPolyDataReader>::New();
  this->xmlreader->SetFileName(this->fullName.c_str());
}
//----------------------------------------------------------------------------
void TestStruct::DeleteXMLPolyDataReader()
{
  this->xmlreader->SetFileName(NULL);
  this->xmlreader = NULL;
}
//----------------------------------------------------------------------------
void TestStruct::CreatePartitioner_Particles()
{
  testDebugMacro( "Creating Partitioner " << this->myRank << " of " << this->numProcs );
  this->partitioner = vtkSmartPointer<vtkParticlePartitionFilter>::New();
  this->partitioner->SetController(this->controller);
}
//----------------------------------------------------------------------------
void TestStruct::CreatePartitioner_Mesh()
{
  testDebugMacro( "Creating Partitioner " << this->myRank << " of " << this->numProcs );
  this->partitioner = vtkSmartPointer<vtkMeshPartitionFilter>::New();
  this->partitioner->SetController(this->controller);
}
//----------------------------------------------------------------------------
double TestStruct::UpdatePartitioner()
{
  vtkSmartPointer<vtkTimerLog> partitiontimer = vtkSmartPointer<vtkTimerLog>::New();
  partitiontimer->StartTimer();
  //
  vtkStreamingDemandDrivenPipeline *partition_sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(this->partitioner->GetExecutive());
  partition_sddp->UpdateDataObject();
  testDebugMacro( "Partition DataObject Updated " << this->myRank << " of " << this->numProcs );
  partition_sddp->SetUpdateExtent(0, this->myRank, this->numProcs, 0);
  partition_sddp->UpdateInformation();
  testDebugMacro( "Partition Information Updated " << this->myRank << " of " << this->numProcs );
  testDebugMacro( "Partition Update coming " << this->myRank << " of " << this->numProcs );
  partition_sddp->Update();
  testDebugMacro( "Partition Updated " << this->myRank << " of " << this->numProcs );
  this->controller->Barrier();
  partitiontimer->StopTimer();
  double partition_elapsed = partitiontimer->GetElapsedTime();
  testDebugMacro( "Partition completed in " << partition_elapsed << " seconds" );
  return partition_elapsed;
}
//----------------------------------------------------------------------------
void TestStruct::DeletePartitioner()
{
  this->partitioner->SetInputConnection(NULL);
  this->partitioner = NULL;
}
//----------------------------------------------------------------------------
int TestStruct::RenderPieces(int argc, char **argv, vtkPolyData *OutputData)
{
    //
    vtkSmartPointer<vtkRenderer>                ren = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow>      renWindow = vtkSmartPointer<vtkRenderWindow>::New();
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    vtkSmartPointer<vtkInteractorStyleSwitch> style = vtkSmartPointer<vtkInteractorStyleSwitch>::New();
    iren->SetRenderWindow(renWindow);
    iren->SetInteractorStyle(style);
    style->SetCurrentStyleToTrackballCamera();
    ren->SetBackground(0.1, 0.1, 0.1);
    renWindow->SetSize(windowSize);
    renWindow->AddRenderer(ren);
    //
    // To make display of ghost cells and boundary regions better, find the centre of
    // all the pieces and use that to apply a transform to the actors to shift them
    // away from the centre so the edges don't touch and overlapping boundary cells are visible
    //
    double centre[3], midpoint[3]={0,0,0};
    for (int i=0; i<numProcs; i++) {
        vtkBoundingBox *box = partitioner->GetPartitionBoundingBox(i);
        box->GetCenter(centre);
        for (int d=0; d<3; d++) midpoint[d] += centre[d]/numProcs;
    }
    //
    for (int i=0; i<numProcs; i++) {
        vtkSmartPointer<vtkPolyData> pd;
        if (i==0) {
            pd = OutputData;
        }
        else {
            pd = vtkSmartPointer<vtkPolyData>::New();
            testDebugMacro("data receiving at "<<myRank);
            controller->Receive(pd, i, DATA_SEND_TAG);
            testDebugMacro("data received at "<<myRank<<" from "<< i);
        }
//        pd->PrintSelf(std::cout, vtkIndent(0));
        vtkSmartPointer<vtkPolyDataMapper>       mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        vtkSmartPointer<vtkActor>                 actor = vtkSmartPointer<vtkActor>::New();
        mapper->SetInputData(pd);
        mapper->SetImmediateModeRendering(1);
        // scalars
        mapper->SetColorModeToMapScalars();
        if (scalarMode==0) {
            mapper->SetScalarModeToUsePointFieldData();
        }
        else {
            mapper->SetScalarModeToUseCellFieldData();
        }
        testDebugMacro("setting scalar colours to " << scalarName.c_str());
        mapper->SelectColorArray(scalarName.c_str());
        mapper->SetUseLookupTableScalarRange(0);
        mapper->SetScalarRange(scalarRange[0], scalarRange[1]);
        mapper->SetInterpolateScalarsBeforeMapping(0);
        //
        actor->SetMapper(mapper);
        actor->GetProperty()->SetPointSize(2);
        ren->AddActor(actor);
        // move each actor away from the midpoint so we can see ghost cells better
        vtkBoundingBox *box = partitioner->GetPartitionBoundingBox(i);
        box->GetCenter(centre);
        vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        transform->PostMultiply();
        transform->Translate(
                -actor_shift*(midpoint[0]-centre[0]),
                -actor_shift*(midpoint[1]-centre[1]),
                -actor_shift*(midpoint[2]-centre[2]));
        actor->SetUserTransform(transform);
        //
        if (cameraSet) {
            ren->GetActiveCamera()->SetPosition(cameraPosition);
            ren->GetActiveCamera()->SetFocalPoint(cameraFocus);
            ren->GetActiveCamera()->SetViewUp(cameraViewUp);
            ren->ResetCameraClippingRange();
        }
        else {
            ren->ResetCamera();
        }
    }
    //
    // Display boxes for each partition
    //
    for (int i=0; i<numProcs; i++) {
        vtkBoundingBox *box = partitioner->GetPartitionBoundingBox(i);
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
        transform->Translate(
                -actor_shift*(midpoint[0]-centre[0]),
                -actor_shift*(midpoint[1]-centre[1]),
                -actor_shift*(midpoint[2]-centre[2]));
        bactor->SetUserTransform(transform);
    }

    testDebugMacro( "Process Id : " << myRank << " About to Render" );
    renWindow->Render();

    int retVal = vtkRegressionTester::Test(argc, argv, renWindow, 10);

    if ( retVal == vtkRegressionTester::DO_INTERACTOR) {
        iren->Start();
    }
    bool ok = (retVal==vtkRegressionTester::PASSED);
    testDebugMacro( "Process Id : " << myRank << " Rendered " << (ok?"Pass":"Fail"));
    return retVal;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
