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
typedef vtkZoltanV1PartitionFilter::ProcessExchangeVariables Exchange;
//----------------------------------------------------------------------------
// #define JB_DEBUG__
#if defined JB_DEBUG__
#define OUTPUTTEXT(a) std::cout <<(a); std::cout.flush();

  #undef vtkDebugMacro
  #define vtkDebugMacro(a)  \
  { \
    if (this->UpdatePiece>=0) { \
      vtkOStreamWrapper::EndlType endl; \
      vtkOStreamWrapper::UseEndl(endl); \
      vtkOStrStreamWrapper vtkmsg; \
      vtkmsg << "P(" << this->UpdatePiece << "): " a << "\n"; \
      OUTPUTTEXT(vtkmsg.str()); \
      vtkmsg.rdbuf()->freeze(0); \
    } \
  }

  #undef  vtkErrorMacro
  #define vtkErrorMacro(a) vtkDebugMacro(a)  
#endif
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkParticlePartitionFilter);
//----------------------------------------------------------------------------
// vtkParticlePartitionFilter :: implementation 
//----------------------------------------------------------------------------
vtkParticlePartitionFilter::vtkParticlePartitionFilter()
{
  this->GhostCellOverlap          = 0.0;
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
int vtkParticlePartitionFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkPointSet     *output = vtkPointSet::GetData(outputVector,0);
  vtkInformation  *inInfo = inputVector[0]->GetInformationObject(0);
  vtkPointSet      *input = vtkPointSet::GetData(inputVector[0]);
  //
  this->UpdatePiece     = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  int ghostLevel        = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
  vtkDebugMacro(<<"Partition filter " << this->UpdatePiece << " Ghost level " << ghostLevel);
  //
#ifdef VTK_USE_MPI
  vtkMPICommunicator *communicator = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  MPI_Comm mpiComm = MPI_COMM_NULL;
  if (communicator) {
    mpiComm = *(communicator->GetMPIComm()->GetHandle());
  }
#else
  int mpiComm = 0;
#endif

  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  // Get input
  vtkIdType       numPoints = input->GetNumberOfPoints();
  vtkIdType        numCells = input->GetNumberOfCells();
  vtkDataArray    *inPoints = numPoints>0 ? input->GetPoints()->GetData() : NULL;
  vtkDebugMacro(<<"Partitioning on " << this->UpdatePiece << " Points Input : " << numPoints);

  // collect bounding boxes for all MPI ranks, we'll use it later to clamp the BSP limits
  // do this even if just one process as it also sets up the extent translator
  vtkBoundingBox globalBounds;
  this->InitBoundingBoxes(input,globalBounds);

  // if only one process, we can exit quietly just passing data through
  if (this->UpdateNumPieces==1) {
    output->ShallowCopy(input);
    return 1;
  }

  //
  // we make a temp copy of the input so we can add Ids if necessary
  //
  vtkSmartPointer<vtkPointSet> inputCopy = input->NewInstance();
  inputCopy->ShallowCopy(input);

  // Setup output points/cells
  vtkSmartPointer<vtkPoints>   outPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray>    cells = vtkSmartPointer<vtkCellArray>::New();

  // if input had 0 points, make sure output is still setup correctly (float/double?)
  // collective exchanges will break if this is wrong as we may still receive data
  int pointstype = this->GatherDataTypeInfo(input->GetPoints());
  outPoints->SetDataType(pointstype);
  output->SetPoints(outPoints);

  //--------------------------------------------------------------
  // Use Zoltan library to re-partition data in parallel
  // declare pointers which hold array info
  //--------------------------------------------------------------
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport=0, numExport=0;
  ZOLTAN_ID_PTR importGlobalGids = NULL;
  ZOLTAN_ID_PTR importLocalGids  = NULL; 
  ZOLTAN_ID_PTR exportGlobalGids = NULL;
  ZOLTAN_ID_PTR exportLocalGids  = NULL;
  int *importProcs  = NULL;
  int *importToPart = NULL;
  int *exportProcs  = NULL;
  int *exportToPart = NULL;

  //
  // A structure we can pass in/out of zoltan callbacks holding the info we
  // need for our allocation and exchange of particles
  //
  Exchange mesh;

  vtkDebugMacro(<<"Initializing Zoltan on " << this->UpdatePiece);
  float ver;
  int rc = Zoltan_Initialize(0, NULL, &ver);
  if (rc != ZOLTAN_OK){
    printf("Zoltan initialization failed ...\n");
    return 0;
  }
  vtkDebugMacro(<<"Zoltan Initialized on " << this->UpdatePiece);

  //
  // Setup mesh structure as a user parameter to zoltan 
  // with partition info (H5Part reader can generate this)
  //
  mesh.ProcessRank              = this->UpdatePiece;
  mesh.Input                    = inputCopy;
  mesh.Output                   = output;
  mesh.InputNumberOfLocalPoints = numPoints;
  mesh.InputNumberOfLocalCells  = numCells;
  mesh.InputPointsData          = inPoints ? inPoints->GetVoidPointer(0) : NULL;
  mesh.OutputPoints             = outPoints;
  mesh.OutPointCount            = 0;

  //
  // if a process has zero points, we need to make dummy point data arrays to allow 
  // space for when data gets sent in from other processes in the zoltan unpack function 
  // This also stops hangs during collective operations by ensuring all ranks participate
  //
  vtkSmartPointer<vtkPointData> PointDataCopy = inputCopy->GetPointData();
  this->SetupFieldArrayPointers(PointDataCopy, mesh, mesh.NumberOfPointFields);

  //
  // Global Ids : always do them after other point arrays 
  //
  std::string IdsName;
  if (this->IdChannelArray) {
    IdsName = this->IdChannelArray;
  }
  if (IdsName.empty() || IdsName==std::string("Not available")) {
    IdsName = "PPF_PointIds";
  } 

  vtkSmartPointer<vtkDataArray> Ids = NULL;
  Ids = PointDataCopy->GetArray(IdsName.c_str());
  if (!Ids) {
    // Try loading the user supplied global ids.
    Ids = PointDataCopy->GetGlobalIds();
  }
  if (!Ids) {
    // Generate our own since none exist
    Ids = this->GenerateGlobalIds(numPoints, numCells, IdsName.c_str(), &mesh);
    inputCopy->GetPointData()->AddArray(Ids);
    // and increment the mesh field count
    mesh.NumberOfPointFields++;
  }

  //***************************************************************
  //* Create a Zoltan library structure for this instance of load
  //* balancing.  Set the parameters and query functions that will
  //* govern the library's calculation.  See the Zoltan User's
  //* Guide for the definition of these and many other parameters.
  //***************************************************************

  zz = Zoltan_Create(mpiComm); 

  // we don't need any debug info
  Zoltan_Set_Param(zz, "RCB_OUTPUT_LEVEL", "0");
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");

  // Method for subdivision
  Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");
  Zoltan_Set_Param(zz, "LB_METHOD",   "RCB");
  //  Zoltan_Set_Param(zz, "LB_METHOD", "PARMETIS");

  // Global and local Ids are a single integer
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");

  // divide into N global and M local partitions
  std::stringstream global;
  global << this->UpdateNumPieces << ends;
  std::stringstream local;
  local << 1 << ends;

  Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", global.str().c_str());
//  Zoltan_Set_Param(zz, "NUM_LOCAL_PARTS",  local.str().c_str());

  // All points have the same weight
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  // RCB parameters
  // Zoltan_Set_Param(zz, "PARMETIS_METHOD", "PARTKWAY");
  Zoltan_Set_Param(zz, "RCB_RECOMPUTE_BOX", "0");
  Zoltan_Set_Param(zz, "REDUCE_DIMENSIONS", "0");

  // Don't allow very extended regions
  std::stringstream aspect;
  aspect << this->MaxAspectRatio << std::ends;  
  Zoltan_Set_Param(zz, "RCB_MAX_ASPECT_RATIO", aspect.str().c_str());

  // we need the cuts to get BBoxes for partitions later
  Zoltan_Set_Param(zz, "KEEP_CUTS", "1");

  // don't allow points on cut to be in different partitions
  // not likely/useful for particle data anyway
  Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1"); 

  // Let Zoltan do the load balance step automatically
  // particles will be transferred as required between processes
  Zoltan_Set_Param(zz, "AUTO_MIGRATE", "1");  

  //
  // Query functions, to provide geometry to Zoltan 
  //
  Zoltan_Set_Num_Obj_Fn(zz,    get_number_of_objects_points, &mesh);
  Zoltan_Set_Obj_List_Fn(zz,   get_object_list_points,       &mesh);
  Zoltan_Set_Num_Geom_Fn(zz,   get_num_geometry,             &mesh);
  if (pointstype==VTK_FLOAT) {
    vtkDebugMacro(<<"Using float data pointers ");
    Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list<float>, &mesh);
  }
  else if (pointstype==VTK_DOUBLE) {
    vtkDebugMacro(<<"Using double data pointers ");
    Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list<double>, &mesh);
  }

  //
  // Register functions for packing and unpacking data
  // by migration tools.  
  // GCC has trouble resolving the templated function pointers, so we explicitly
  // declare the types and then cast them as args
  //    
  
  typedef int  (*zsize_fn) (void *, int , int , ZOLTAN_ID_PTR , ZOLTAN_ID_PTR , int *);
  typedef void (*zpack_fn) (void *, int , int , ZOLTAN_ID_PTR , ZOLTAN_ID_PTR , int , int , char *, int *);
  typedef void (*zupack_fn)(void *, int , ZOLTAN_ID_PTR , int , char *, int *);
  typedef void (*zprem_fn) (void *, int , int , int , ZOLTAN_ID_PTR , ZOLTAN_ID_PTR , int *, int *, int , ZOLTAN_ID_PTR ,
    ZOLTAN_ID_PTR , int *, int *, int *);

  if (pointstype==VTK_FLOAT) {
    zsize_fn  f1 = zoltan_obj_size_func_points<float>;
    zpack_fn  f2 = zoltan_pack_obj_func_points<float>;
    zupack_fn f3 = zoltan_unpack_obj_func_points<float>;
    zprem_fn  f4 = zoltan_pre_migrate_func_points<float>; 
    Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &mesh); 
    Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &mesh); 
    Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &mesh); 
    Zoltan_Set_Fn(zz, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &mesh); 
  }
  else if (pointstype==VTK_DOUBLE) {
    zsize_fn  f1 = zoltan_obj_size_func_points<double>;
    zpack_fn  f2 = zoltan_pack_obj_func_points<double>;
    zupack_fn f3 = zoltan_unpack_obj_func_points<double>;
    zprem_fn  f4 = zoltan_pre_migrate_func_points<double>;
    Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &mesh);
    Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &mesh);
    Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &mesh);
    Zoltan_Set_Fn(zz, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &mesh);
  }

  //
  // Check the input to see if it has a bounds translator already initialized
  // with partition info (H5Part reader can generate this)
  //
  vtkExtentTranslator *translator = inInfo ? vtkExtentTranslator::SafeDownCast(
    inInfo->Get(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR())) : NULL;
  vtkBoundsExtentTranslator *input_bet = vtkBoundsExtentTranslator::SafeDownCast(translator);
  // if the extent translator has not been initialized well - don't use it
  if (input_bet && input_bet->GetNumberOfPieces()==0) {
    input_bet = NULL;
  }

  //
  // if the input had a BoundsExtentTranslator, then we don't need to load-balance particles
  // we only need to exchange halo particles.
  //
  if (!input_bet) {
    //
    // Zoltan can now partition our particles. 
    // After this returns, we have redistributed particles and the Output holds
    // the list of correct points/fields etc for each process
    //
    rc = Zoltan_LB_Partition(zz, // input (all remaining fields are output)
          &changes,              // 1 if partitioning was changed, 0 otherwise 
          &numGidEntries,        // Number of integers used for a global ID
          &numLidEntries,        // Number of integers used for a local ID
          &numImport,            // Number of vertices to be sent to me
          &importGlobalGids,     // Global IDs of vertices to be sent to me
          &importLocalGids,      // Local IDs of vertices to be sent to me
          &importProcs,          // Process rank for source of each incoming vertex
          &importToPart,         // New partition for each incoming vertex
          &numExport,            // Number of vertices I must send to other processes
          &exportGlobalGids,     // Global IDs of the vertices I must send
          &exportLocalGids,      // Local IDs of the vertices I must send
          &exportProcs,          // Process to which I send each of the vertices
          &exportToPart);        // Partition to which each vertex will belong

    if (rc != ZOLTAN_OK){
      printf("Zoltan_LB_Partition NOT OK...\n");
      MPI_Finalize();
      Zoltan_Destroy(&zz);
      exit(0);
    }

#ifdef EXTRA_ZOLTAN_DEBUG
    vtkDebugMacro(<<"Partitioning complete on " << this->UpdatePiece << 
      " pack_count : " << pack_count <<
      " size_count : " << size_count <<
      " unpack_count : " << unpack_count 
     );
#endif

    if (!this->ExchangeCells) {
      Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                          &importProcs, &importToPart);
      Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                          &exportProcs, &exportToPart);
    }
  }
  else {
    //
    // If we skipped the zoltan repartitioning, then do a copy (setup everything)
    // just like would have taken place at the start of the load balancing
    //
    if (pointstype==VTK_FLOAT) {
      zprem_fn  f4 = zoltan_pre_migrate_func_points<float>; 
      f4(&mesh, 0, 0, 0, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL); 
    }
    else if (pointstype==VTK_DOUBLE) {
      zprem_fn  f4 = zoltan_pre_migrate_func_points<double>;
      f4(&mesh, 0, 0, 0, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL); 
    }
  }

  //
  // Get the partition bounding boxes from Zoltan if we don't have them
  //  
  this->ExtentTranslator->SetNumberOfPieces(this->UpdateNumPieces);
  if (!input_bet) {
    for (int p=0; p<this->UpdateNumPieces; p++) {
      double bounds[6];
      int ndim;
      if (ZOLTAN_OK==Zoltan_RCB_Box(zz, p, &ndim, &bounds[0], &bounds[2], &bounds[4], &bounds[1], &bounds[3], &bounds[5])) {
        if (bounds[0]==-DBL_MAX) { bounds[0] = globalBounds.GetMinPoint()[0]; }
        if (bounds[1]== DBL_MAX) { bounds[1] = globalBounds.GetMaxPoint()[0]; }
        if (bounds[2]==-DBL_MAX) { bounds[2] = globalBounds.GetMinPoint()[1]; }
        if (bounds[3]== DBL_MAX) { bounds[3] = globalBounds.GetMaxPoint()[1]; }
        if (bounds[4]==-DBL_MAX) { bounds[4] = globalBounds.GetMinPoint()[2]; }
        if (bounds[5]== DBL_MAX) { bounds[5] = globalBounds.GetMaxPoint()[2]; }
        vtkBoundingBox box(bounds);
        this->BoxList.push_back(box);
        this->ExtentTranslator->SetBoundsForPiece(p, bounds);
      }
    }
  }
  else {
    for (int p=0; p<this->UpdateNumPieces; p++) {
      vtkBoundingBox box;  
      box.SetBounds(input_bet->GetBoundsForPiece(p));
      this->BoxList.push_back(box);
      this->ExtentTranslator->SetBoundsForPiece(p, input_bet->GetBoundsForPiece(p));
    }
  }

  vtkSmartPointer<vtkIntArray> boundary;
  PartitionInfo partitioninfo;
  if (this->ExchangeCells) {
    boundary = this->BuildCellToProcessList(mesh.Input, 
      partitioninfo,
      mesh.ProcessOffsetsCellId,
      numExport,            // Number of vertices I must send to other processes
      exportGlobalGids,     // Global IDs of the vertices I must send
      exportLocalGids,      // Local IDs of the vertices I must send
      exportProcs           // Process to which I send each of the vertices
      );

    this->SetupFieldArrayPointers(mesh.Output->GetCellData(), mesh, mesh.NumberOfCellFields);
    //
    // now we have a map of cells to processId, so do a collective 'invert lists' 
    // operation to compute the global exchange map of who sends cells to who 
    //
    size_t        num_known = partitioninfo.GlobalIds.size(); 
    int           num_found = 0;
    ZOLTAN_ID_PTR found_global_ids = NULL;
    ZOLTAN_ID_PTR found_local_ids  = NULL;
    int          *found_procs      = NULL;
    int          *found_to_part    = NULL;
    //
    rc = Zoltan_Invert_Lists(zz, 
      (int)num_known,
      num_known>0 ? &partitioninfo.GlobalIds[0] : NULL,
      num_known>0 ? &partitioninfo.LocalIds[0]  : NULL,
      num_known>0 ? &partitioninfo.Procs[0]     : NULL,
      num_known>0 ? &partitioninfo.Procs[0]     : NULL,
      &num_found,
      &found_global_ids,
      &found_local_ids,
      &found_procs,
      &found_to_part); 
    //
    if (rc != ZOLTAN_OK){
      printf("Zoltan_LB_Partition NOT OK...\n");
      MPI_Finalize();
      Zoltan_Destroy(&zz);
      exit(0);
    }

    //
    // we need to change the callback functions for migration to use cell based
    // functions we declare instead of the point exchange ones
    //
    zsize_fn  f1 = zoltan_obj_size_func_cell;
    zpack_fn  f2 = zoltan_pack_obj_func_cell;
    zupack_fn f3 = zoltan_unpack_obj_func_cell;
    zprem_fn  f4 = zoltan_pre_migrate_func_cell<float>;
    Zoltan_Set_Fn(zz, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &mesh); 
    Zoltan_Set_Fn(zz, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &mesh); 
    Zoltan_Set_Fn(zz, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &mesh); 
    Zoltan_Set_Fn(zz, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &mesh); 

    //
    // Perform the cell exchange
    //
    rc = Zoltan_Migrate (zz,
      (int)num_found,
      found_global_ids,
      found_local_ids,
      found_procs,
      found_to_part,
      (int)num_known,
      num_known>0 ? &partitioninfo.GlobalIds[0] : NULL,
      num_known>0 ? &partitioninfo.LocalIds[0]  : NULL,
      num_known>0 ? &partitioninfo.Procs[0]     : NULL,
      num_known>0 ? &partitioninfo.Procs[0]     : NULL
      );

  }

  if (this->ExchangeHaloPoints) {
    //
    // Set the halo/ghost regions we need around each process bounding box
    //  
    if (input_bet && input_bet->GetBoundsHalosPresent()) {
      this->ExtentTranslator->SetBoundsHalosPresent(1);
      for (int p=0; p<this->UpdateNumPieces; p++) {
        vtkBoundingBox box;  
        box.SetBounds(input_bet->GetBoundsHaloForPiece(p));
        this->BoxListWithHalo.push_back(box);
        this->ExtentTranslator->SetBoundsHaloForPiece(p,input_bet->GetBoundsHaloForPiece(p));
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
      mesh.Output->GetPointData()->GetArray(IdsName.c_str()));
    if (!newIds || newIds->GetNumberOfTuples()!=mesh.OutputPoints->GetNumberOfPoints()) {
      vtkErrorMacro(<<"Fatal : Ids on migrated data corrupted");
      return 0;
    }

    PartitionInfo GhostIds;
    this->FindPointsInHaloRegions(mesh.OutputPoints, newIds, GhostIds);

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
    rc = Zoltan_Invert_Lists(zz, 
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

    if (rc != ZOLTAN_OK){
      printf("Zoltan_LB_Partition NOT OK...\n");
      MPI_Finalize();
      Zoltan_Destroy(&zz);
      exit(0);
    }

    //
    // Before sending, we need to change the pre-migrate function as we are now adding
    // extra ghost cells and not starting our lists from a clean slate.
    //
    if (pointstype==VTK_FLOAT) {
      zprem_fn f4 = zoltan_pre_migrate_func_halo<float>;
      Zoltan_Set_Fn(zz, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &mesh);
    }
    else if (pointstype==VTK_DOUBLE) {
      zprem_fn f4 = zoltan_pre_migrate_func_halo<double>;
      Zoltan_Set_Fn(zz, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &mesh);
    }

    //
    // Now we can actually send ghost particles between processes
    //
    rc = Zoltan_Migrate (zz,
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
    GhostArray->SetNumberOfTuples(mesh.OutputNumberOfLocalPoints + num_found);
    GhostArray2->SetName("GhostLevels");
    GhostArray2->SetNumberOfComponents(1);
    GhostArray2->SetNumberOfTuples(mesh.OutputNumberOfLocalPoints + num_found);
    unsigned char *ghost = GhostArray->GetPointer(0);
    int          *ghost2 = GhostArray2->GetPointer(0);
    for (vtkIdType i=0; i<mesh.OutputNumberOfLocalPoints + num_found; i++) {
      if (i<mesh.OutputNumberOfLocalPoints) {
        ghost[i]  = 0;
        ghost2[i] = 0;
      }
      else {
        ghost[i]  = 1;
        ghost2[i] = 1;
      }
    }
    output->GetPointData()->AddArray(GhostArray);
    output->GetPointData()->AddArray(GhostArray2);

    //
    //
    //
    vtkDebugMacro(<<"Process " << this->UpdatePiece << " Points Output : " << mesh.OutPointCount);
    //  for (int i=0; i<mesh.NumberOfPointFields; i++) {
    //    vtkDataArray *darray = output->GetPointData()->GetArray(i);
    //    vtkDebugMacro(<<"Process " << this->UpdatePiece << " Array Output : " << darray->GetNumberOfTuples());
    //  }

    //
    // If polydata create Vertices for each point
    //
    if (vtkPolyData::SafeDownCast(output)) {
      vtkIdType *arraydata = cells->WritePointer(mesh.OutPointCount, 2*mesh.OutPointCount);
      for (int i=0; i<mesh.OutPointCount; i++) {
        arraydata[i*2]   = 1;
        arraydata[i*2+1] = i;
      }
      vtkPolyData::SafeDownCast(output)->SetVerts(cells);
    }
    //
    //*****************************************************************
    // Free the arrays allocated by Zoltan_LB_Partition, and free
    // the storage allocated for the Zoltan structure.
    //*****************************************************************
    //
    Zoltan_Destroy(&zz);

  }

  this->Controller->Barrier();
  timer->StopTimer();
  vtkDebugMacro(<<"Particle partitioning : " << timer->GetElapsedTime() << " seconds");
  return 1;
}
//----------------------------------------------------------------------------
