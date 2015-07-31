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
//
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
#ifdef _WIN32
unsigned long int random_seed()
{
  LARGE_INTEGER lpPerformanceCount;
  QueryPerformanceCounter(&lpPerformanceCount);
  long int seed = lpPerformanceCount.LowPart + lpPerformanceCount.HighPart;
  srand(seed);
  return seed;
}
#else
unsigned long int random_seed()
{
  unsigned int seed;
  struct timeval tv;
  FILE *devrandom;
  if ((devrandom = fopen("/dev/random","r")) == NULL) {
    gettimeofday(&tv,0);
    seed = tv.tv_sec + tv.tv_usec;
  } 
  else {
    if (fread(&seed,sizeof(seed),1,devrandom) == 1) {
      fclose(devrandom);
    } 
    else {
      gettimeofday(&tv,0);
      seed = tv.tv_sec + tv.tv_usec;
    }
  }
  srandom(seed);
  return seed;
}
#endif
//----------------------------------------------------------------------------
void SpherePoints(int n, float radius, float X[]) {
  double x, y, z, w, t;
  Random r(12345);
  double rmin=1E6, rmax=-1E6;
  for(int i=0; i<n; i++ ) {
    double r1 = r.nextNumber(); // double(rand())/RAND_MAX;
    double r2 = r.nextNumber(); // double(rand())/RAND_MAX;
    rmin = std::min(rmin, r1);
    rmin = std::min(rmin, r2);
    rmax = std::max(rmax, r1);
    rmax = std::max(rmax, r2);
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

  double rmin=1E6, rmax=-1E6;
  for(int i=0; i<n; i++ ) {
    double r1 = uniform_dist(generator);
    double r2 = uniform_dist(generator);
    double r3 = uniform_dist(generator);
    X[3*i+0] = static_cast<float>(r1);
    X[3*i+1] = static_cast<float>(r2);
    X[3*i+2] = static_cast<float>(r3);
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

  // uncomment this to wait for debugger attach
//   DEBUG_WAIT
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
  test.ghostOverlap = GetParameter<double>("-ghost_region", "Ghost Region", argc, argv, 0.0, test.myRank, unused);
  GetArrayParameter<double>("-gridSpacing", "Grid Spacing", test.gridSpacing, 3, argc, argv, test.myRank);
  GetArrayParameter<int>("-gridResolution", "Grid Resolution", test.gridResolution, 3, argc, argv, test.myRank);
  test.maxN = GetParameter<int>("-neighbours", "Fixed Neighbours", argc, argv, 0, test.myRank, test.fixNeighbours);
  test.massScalars = GetParameter<std::string>("-massScalars", "Mass Scalar Array", argc, argv, "", test.myRank, unused);
  test.densityScalars = GetParameter<std::string>("-densityScalars", "Density Scalar Array", argc, argv, "", test.myRank, unused);
  test.expectedN = GetParameter<vtkIdType>("-expectedparticles", "Simulation Particles", argc, argv, 0, test.myRank, unused);

  //
  // Test/Display of results
  //
  test.scalarname = GetParameter<std::string>("-scalar", "Testing Scalar Array", argc, argv, "", test.myRank, unused);
  test.contourVal = GetParameter<double>("-contour", "Contour Value", argc, argv, 0.0, test.myRank, unused);
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
  this->xmlreader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
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
        cout<<"Here "<<1<<"by >>> "<<this->myRank<<endl;
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
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
