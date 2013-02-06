/*=========================================================================

  Module                  : vtkMeshPartitionFilter.cxx

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
#include "vtkPKdTree.h"
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
#include "vtkMeshPartitionFilter.h"
//
#include <sstream>
//
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <numeric>
#include <algorithm>
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkMeshPartitionFilter);
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Function Type: Pre migration callback for cell exchange
// the cell exchange is controlled manually. We declare what we need to exchange
// and use this function to ensure all memory is allocated correctly for the receive
// to take place and unpack each object
//----------------------------------------------------------------------------
template <typename T>
void vtkMeshPartitionFilter::zoltan_pre_migrate_function_cell(void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  // How many cells will we have at the end
  vtkIdType OutputNumberOfLocalCells = callbackdata->Input->GetNumberOfCells();
  vtkIdType OutputNumberOfFinalCells = OutputNumberOfLocalCells + num_import - num_export;
  vtkIdType       ReservedPointCount = callbackdata->Output->GetNumberOfPoints();
  callbackdata->OutputPointsData     = callbackdata->Output->GetPoints()->GetData()->GetVoidPointer(0);
  // we'll need enough space to handle the largest cells
  callbackdata->MaxCellSize     = callbackdata->Input->GetMaxCellSize();
  callbackdata->OutputCellArray = vtkSmartPointer<vtkCellArray>::New();
  callbackdata->OutputCellArray->Allocate(OutputNumberOfLocalCells*(callbackdata->MaxCellSize+1));
  //
  vtkCellData *inCD  = callbackdata->Input->GetCellData();
  vtkCellData *outCD = callbackdata->Output->GetCellData();
  outCD->CopyAllocate(inCD, OutputNumberOfFinalCells);
  //
  callbackdata->self->InitializeFieldDataArrayPointers(callbackdata, inCD, outCD, OutputNumberOfFinalCells);

  std::vector<bool> local(OutputNumberOfLocalCells, true);
  for (vtkIdType i=0; i<num_export; i++) {
    local[export_global_ids[i]-callbackdata->ProcessOffsetsCellId[callbackdata->ProcessRank]] = false;    
  }
  //
  vtkIdType npts, *pts, newpts[32];
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(callbackdata->Input);
//  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(callbackdata->Input);
  vtkPolyData         *pdata2 = vtkPolyData::SafeDownCast(callbackdata->Output);
//  vtkUnstructuredGrid *udata2 = vtkUnstructuredGrid::SafeDownCast(callbackdata->Output);
  if (pdata2) {
    pdata2->SetPolys(callbackdata->OutputCellArray);
  }
  //
  callbackdata->OutCellCount = 0;
  for (vtkIdType cellId=0; cellId<OutputNumberOfLocalCells; cellId++) {
    if (local[cellId]) {
      // copy cell data from old to new datasets
      outCD->CopyData(inCD, cellId, callbackdata->OutCellCount);
      // copy cell point Ids to new dataset, 
      if (pdata) { 
        int ctype = pdata->GetCellType(cellId);
        pdata->GetCellPoints(cellId, npts, pts); 
        for (int i=0; i<npts; i++) {
          if (callbackdata->LocalToLocalIdMap[pts[i]]!=-1) newpts[i] = callbackdata->LocalToLocalIdMap[pts[i]];
          else {
/*
            // a point which has been migrated has been referenced. 
            // It can only happen on a boundary cell which has had some points sent away
            // option 1. Make a local copy of the point we've already sent, 
            // option 2. Send the other points too and make the cell remote. 
            // easiest solution is 1 currently.
            if (callbackdata->OutPointCount>=ReservedPointCount) {
              // output points should have been resized to fit all existing and migrated points, 
              // but we need to allow a small amount (try 5%) extra for duplicated boundary cell points
              ReservedPointCount = callbackdata->OutPointCount * 1.05;
              callbackdata->Output->GetPoints()->GetData()->Resize(ReservedPointCount);
              callbackdata->OutputPointsData = callbackdata->Output->GetPoints()->GetData()->GetVoidPointer(0);
            }
            callbackdata->Output->GetPointData()->CopyData(callbackdata->Input->GetPointData(), pts[i], callbackdata->OutPointCount);
            newpts[i] = callbackdata->LocalToLocalIdMap[pts[i]] = callbackdata->OutPointCount;
            memcpy(&((T*)(callbackdata->OutputPointsData))[callbackdata->OutPointCount*3], &((T*)(callbackdata->InputPointsData))[pts[i]*3], sizeof(T)*3);
            callbackdata->OutPointCount++;
*/
          }
        }
        pdata2->InsertNextCell(ctype, npts, newpts);
      }
      //else if (udata) { 
      //  int ctype = udata->GetCellType(cellId);
      //  udata->GetCellPoints(cellId, npts, pts); 
      //  udata2->InsertNextCell(ctype, npts, pts);
      //}
      callbackdata->OutCellCount++;
    }
  }
  // make sure final point count is actual point count, not the reserved amount
//  callbackdata->Output->GetPoints()->SetNumberOfPoints(callbackdata->OutPointCount);
}
//----------------------------------------------------------------------------
// Zoltan callback function : returns size of each cell and all its data
//----------------------------------------------------------------------------
int vtkMeshPartitionFilter::zoltan_obj_size_function_cell(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
{
  INC_SIZE_COUNT
  vtkZoltanV1PartitionFilter::CallbackData *callbackdata = (vtkZoltanV1PartitionFilter::CallbackData*)data;
  *ierr = ZOLTAN_OK;
  // return the size of the cell data + number of points in cell + point Ids
  vtkIdType GID = *global_id;
  vtkIdType LID = GID - callbackdata->ProcessOffsetsCellId[callbackdata->ProcessRank];
  //
  vtkIdType npts, *pts;
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(callbackdata->Input);
//  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(callbackdata->Input);
  //
  int ctype = pdata->GetCellType(LID);
  pdata->GetCellPoints(LID,npts,pts);
  int size = 2 + npts; // npts + ctype + pts
  return callbackdata->TotalSizePerId + size*sizeof(vtkIdType);
}
//----------------------------------------------------------------------------
// Zoltan callback function : pack one cell and all its data
//----------------------------------------------------------------------------
void vtkMeshPartitionFilter::zoltan_pack_obj_function_cell(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size, char *buf, int *ierr)
{
  INC_PACK_COUNT
  vtkZoltanV1PartitionFilter::CallbackData *callbackdata = (vtkZoltanV1PartitionFilter::CallbackData*)data;
  vtkIdType GID = *global_id;
  vtkIdType LID = GID - callbackdata->ProcessOffsetsCellId[callbackdata->ProcessRank];
  //
  //
  vtkIdType npts, *pts, newPts[32];
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(callbackdata->Input);
//  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(callbackdata->Input);
  //
  // Copy all Cell data arrays into the buffer provided by zoltan
  //
  for (int i=0; i<callbackdata->NumberOfFields; i++) {
    int asize = callbackdata->MemoryPerTuple[i];
    char *dataptr = (char*)(callbackdata->InputArrayPointers[i]) + asize*LID;
    memcpy(buf, dataptr, asize);
    buf += asize;
  }
  //
  // copy the cell point Ids into the buffer provided by zoltan.
  // before sending, we will convert local to global Ids because when the cell is unpacked on the 
  // remote process, the local Ids are meaningless - there is no way to find where they really came from 
  //
  if (pdata) {
    pdata->GetCellPoints(LID, npts, pts);
    // copy the number of points and cell type so we know what's been sent when we unpack
    newPts[0] = npts;
    newPts[1] = pdata->GetCellType(LID);
    // and the points Ids converted to global Ids
    for (int i=0; i<npts; i++) {
      newPts[i+2] = pts[i] + callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];
    }
    memcpy(buf, newPts, sizeof(vtkIdType)*(npts+2));  
  }
  //else if (udata) { 
  //  throw std::string("Implement this");
  //}

  *ierr = ZOLTAN_OK;
  return;
}
//----------------------------------------------------------------------------
// Zoltan callback function : unpack one cell and all its data
//----------------------------------------------------------------------------
void vtkMeshPartitionFilter::zoltan_unpack_obj_function_cell(void *data, int num_gid_entries,
  ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr)
{
  INC_UNPACK_COUNT
  vtkZoltanV1PartitionFilter::CallbackData *callbackdata = (vtkZoltanV1PartitionFilter::CallbackData*)data;
  //
  vtkIdType GID = *global_id;
//  vtkIdType LID = GID - callbackdata->ProcessOffsetsCellId[callbackdata->ProcessRank];
  //
  vtkIdType npts, *pts, ctype, newPts[32];
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(callbackdata->Input);
//  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(callbackdata->Input);
  vtkPolyData         *pdata2 = vtkPolyData::SafeDownCast(callbackdata->Output);
//  vtkUnstructuredGrid *udata2 = vtkUnstructuredGrid::SafeDownCast(callbackdata->Output);
  //
  for (int i=0; i<callbackdata->NumberOfFields; i++) {
    int asize = callbackdata->MemoryPerTuple[i];
    char *dataptr = (char*)(callbackdata->InputArrayPointers[i]) + asize*(callbackdata->OutCellCount);
    memcpy(dataptr, buf, asize);
    buf += asize;
  }
  // we have received a cell with a list of Global point Ids, 
  // We need to convert each of the global Ids to the local Ids they became 
  // when they were copied into the local point list
  npts  =  reinterpret_cast<vtkIdType*>(buf)[0];
  ctype =  reinterpret_cast<vtkIdType*>(buf)[1];
  pts   = &reinterpret_cast<vtkIdType*>(buf)[2];
  for (int i=0; i<npts; i++) {
    newPts[i] = callbackdata->self->global_to_local_Id(pts[i]);
  }
  if (pdata) {
    pdata2->InsertNextCell(ctype, npts, newPts);
  }
  //else if (udata) { 
  //  throw std::string("Implement this");
  //}

  callbackdata->OutCellCount++;
  *ierr = ZOLTAN_OK;
  return;
}
//----------------------------------------------------------------------------
// vtkMeshPartitionFilter :: implementation 
//----------------------------------------------------------------------------
vtkMeshPartitionFilter::vtkMeshPartitionFilter()
{
}
//----------------------------------------------------------------------------
vtkMeshPartitionFilter::~vtkMeshPartitionFilter()
{
}
//----------------------------------------------------------------------------
void vtkMeshPartitionFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
//----------------------------------------------------------------------------
int vtkMeshPartitionFilter::RequestData(vtkInformation* info,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  //
  // Distribute points evenly across processes
  //
  this->PartitionPoints(info, inputVector, outputVector);

  if (this->UpdateNumPieces==1) {
    // input has been copied to output 
    return 1;
  }
  //
  // Distribute cells based on the usage of the points already distributed
  //
  this->PartitionCells(info, inputVector, outputVector);

  this->CreatePkdTree();
  this->ExtentTranslator->SetKdTree(this->GetKdtree());

  //*****************************************************************
  // Free the storage allocated for the Zoltan structure.
  //*****************************************************************
  Zoltan_Destroy(&this->ZoltanData);

//  this->Controller->Barrier();
  this->Timer->StopTimer();
  vtkDebugMacro(<<"Particle partitioning : " << this->Timer->GetElapsedTime() << " seconds");
  return 1;
}
//----------------------------------------------------------------------------
int vtkMeshPartitionFilter::PartitionCells(vtkInformation* info,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  PartitionInfo cell_partitioninfo, point_partitioninfo;

  vtkDebugMacro(<<"Entering BuildCellToProcessList");
  if (this->ZoltanCallbackData.PointType==VTK_FLOAT) {
    this->BuildCellToProcessList<float>(this->ZoltanCallbackData.Input, 
      cell_partitioninfo,  // lists of which cells to send to which process
      point_partitioninfo, // list of which points to send to which process
      this->LoadBalanceData.numExport,            // Number of vertices I must send to other processes
      this->LoadBalanceData.exportGlobalGids,     // Global IDs of the vertices I must send
      this->LoadBalanceData.exportLocalGids,      // Local IDs of the vertices I must send
      this->LoadBalanceData.exportProcs           // Process to which I send each of the vertices
      );
  }
  else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
    this->BuildCellToProcessList<double>(this->ZoltanCallbackData.Input, 
      cell_partitioninfo,  // lists of which cells to send to which process
      point_partitioninfo, // list of which points to send to which process
      this->LoadBalanceData.numExport,            // Number of vertices I must send to other processes
      this->LoadBalanceData.exportGlobalGids,     // Global IDs of the vertices I must send
      this->LoadBalanceData.exportLocalGids,      // Local IDs of the vertices I must send
      this->LoadBalanceData.exportProcs           // Process to which I send each of the vertices
      );
  }

  vtkDebugMacro(<<"Freeing Zoltan LB arrays");
  //*****************************************************************
  // Free the arrays allocated by Zoltan_LB_Partition
  // before we do a manual migration.
  //*****************************************************************
  if (this->LoadBalanceData.importGlobalGids) {
    Zoltan_LB_Free_Part(&this->LoadBalanceData.importGlobalGids, &this->LoadBalanceData.importLocalGids, &this->LoadBalanceData.importProcs, &this->LoadBalanceData.importToPart);
    Zoltan_LB_Free_Part(&this->LoadBalanceData.exportGlobalGids, &this->LoadBalanceData.exportLocalGids, &this->LoadBalanceData.exportProcs, &this->LoadBalanceData.exportToPart);
  }
  // after deleting memory, add a barrier to let ranks free as much as possible before the next big allocation
  this->Controller->Barrier();

  //
  // cells that were split over remote processes require another point migration
  // to ensure all the points are available/complete
  //
  vtkDebugMacro(<<"Entering ManualPointMigrate");
  int num_found = this->ManualPointMigrate(point_partitioninfo, false);

  //
  // now we have a map of cells to processId, so do a collective 'invert lists' 
  // operation to compute the global exchange map of who sends cells to who 
  //
  size_t        num_known = cell_partitioninfo.GlobalIds.size(); 
  num_found = 0;
  ZOLTAN_ID_PTR found_global_ids = NULL;
  ZOLTAN_ID_PTR found_local_ids  = NULL;
  int          *found_procs      = NULL;
  int          *found_to_part    = NULL;
  //

  vtkDebugMacro(<<"About to invert lists (cell migration)");
  int zoltan_error = Zoltan_Invert_Lists(this->ZoltanData, 
    (int)num_known,
    num_known>0 ? &cell_partitioninfo.GlobalIds[0] : NULL,
    /*num_known>0 ? &cell_partitioninfo.LocalIds[0]  : */NULL,
    num_known>0 ? &cell_partitioninfo.Procs[0]     : NULL,
    /*num_known>0 ? &cell_partitioninfo.Procs[0]     : */NULL,
    &num_found,
    &found_global_ids,
    &found_local_ids,
    &found_procs,
    &found_to_part); 
  //
  if (zoltan_error != ZOLTAN_OK){
    printf("Zoltan_LB_Partition NOT OK...\n");
    MPI_Finalize();
    Zoltan_Destroy(&this->ZoltanData);
    exit(0);
  }

  this->AllocateFieldArrays(this->ZoltanCallbackData.Output->GetCellData());

  //
  // we need to change the callback functions for migration to use cell based
  // functions we declare instead of the point exchange ones
  //
  zsize_fn  f1 = zoltan_obj_size_function_cell;
  zpack_fn  f2 = zoltan_pack_obj_function_cell;
  zupack_fn f3 = zoltan_unpack_obj_function_cell;
  zprem_fn  f4 = zoltan_pre_migrate_function_cell<float>;
  Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &this->ZoltanCallbackData); 
  Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &this->ZoltanCallbackData); 
  Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &this->ZoltanCallbackData); 
  Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &this->ZoltanCallbackData); 

  //
  // Perform the cell exchange
  //
  vtkDebugMacro(<<"About to Zoltan_Migrate (cells)");
  zoltan_error = Zoltan_Migrate (this->ZoltanData,
    (int)num_found,
    found_global_ids,
    found_local_ids,
    found_procs,
    found_to_part,
    (int)num_known,
    num_known>0 ? &cell_partitioninfo.GlobalIds[0] : NULL,
    /*num_known>0 ? &cell_partitioninfo.LocalIds[0]  : */NULL,
    num_known>0 ? &cell_partitioninfo.Procs[0]     : NULL,
    /*num_known>0 ? &cell_partitioninfo.Procs[0]     : */NULL
    );

  //
  // Release the arrays allocated during Zoltan_Invert_Lists
  //
  vtkDebugMacro(<<"About to Free Zoltan_Migrate (cells)");
  Zoltan_LB_Free_Part(
    &found_global_ids, 
    &found_local_ids, 
    &found_procs, 
    &found_to_part);
  vtkDebugMacro(<<"Done Migration (cells)");

  return 1;
}
//----------------------------------------------------------------------------
// 
// Build a list of cell to process Ids based on the already performed 
// point redistibution. Points which are being sent away will be taking 
// their cells with them. Caution, some points are shared between cells and may 
// go to different processes. For these cells, we must send another copy 
// of the points to the relevant process.
//
template <typename T>
void vtkMeshPartitionFilter::BuildCellToProcessList(     
  vtkDataSet *data, 
  PartitionInfo &cell_partitioninfo, 
  PartitionInfo &point_partitioninfo, 
  int numExport,
  ZOLTAN_ID_PTR exportGlobalGids,
  ZOLTAN_ID_PTR exportLocalGids,
  int *exportProcs)
{
  vtkIdType numPts = data->GetNumberOfPoints();
  vtkIdType numCells = data->GetNumberOfCells();
  vtkIdType j, cellId;

  // we will put the results of the cell tests in these arrays
  cell_partitioninfo.Procs.reserve(numCells/this->UpdateNumPieces);
//  cell_partitioninfo.LocalIds.reserve(numCells);
  cell_partitioninfo.GlobalIds.reserve(numCells/this->UpdateNumPieces);

  // we know that some points on this process are being exported to remote processes
  // so build a point to process map to quickly lookup the process Id from the point Id
  // 1) set all points as local to this process
  std::vector<int> localId_to_process_map(numPts, this->UpdatePiece); 
  // 2) loop over all to be exported and note the destination
  for (vtkIdType i=0; i<numExport; i++) {
    vtkIdType id               = exportGlobalGids[i] - this->ZoltanCallbackData.ProcessOffsetsPointId[this->ZoltanCallbackData.ProcessRank];
    localId_to_process_map[id] = exportProcs[i];
  }

  // when we find cells split across multiple remote processes, we must resend some points
  // to one or other of the remote processes. This vector will store {Id,process} pairs
  // wilst the list is being scanned
  typedef std::pair<vtkIdType, int> process_tuple;
  std::vector<process_tuple> process_vector;

  vtkIdType npts, *pts;
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(data);
  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(data);
  if (!pdata && !udata) {
    vtkErrorMacro(<<"Only Polydata and UnstructuredGrid supported so far");  
  }
  // polydata requires a cell map (verts/lines/polys/strips) to be present before we traverse cells
  if (pdata) pdata->BuildCells();

  // a simple bitmask with one entry per process, used for each cell to count processes for the cell
  std::vector<unsigned char> process_flag(this->UpdateNumPieces,0);;

  //
  // Some points have been migrated already, but are still needed locally, so we
  // must monitor the points copied locally
  //
  vtkIdType   OutputPointCount = this->ZoltanCallbackData.Output->GetNumberOfPoints();
  vtkIdType ReservedPointCount = OutputPointCount;
  
  //
  // for each cell, find if all points are designated as remote and cell needs to be sent away
  //
  for (cellId=0; cellId<numCells; cellId++) {
    // get a pointer to the cell points
    if (pdata) { pdata->GetCellPoints(cellId, npts, pts); }
    else if (udata) { udata->GetCellPoints(cellId, npts, pts); }

    //
    // we need to examine all points of the cell and classify it : there are several possibile actions
    //
    // 1) all points are local                     : keep cell
    // 2) some points are local, some remote       : keep cell, make a local copy of the points already migrated
    // 3) all points on same remote process        : send cell to remote process
    // 4) all points on different remote processes : send cell to remote process, also send previously unsent points
    //
    // use a bit index to mark/mask processes receiving points
    process_flag.assign(this->UpdateNumPieces,0);
    int points_remote = 0;
    for (j=0; j<npts; j++) {
      int process = localId_to_process_map[pts[j]];
      process_flag[process] = 1;
      if (process!=this->UpdatePiece) {
        points_remote++;
      }
    }
    int cellstatus = 0;
    if (points_remote==0) {
      cellstatus = 1;
    }
    else {
      int process_count = std::count(process_flag.begin(), process_flag.end(), 1);
      // all points are on the same remote process
      if (process_count==1) {
        cellstatus = 3;
      }
      // some local, others all on remote process(es)
      else if (process_count>1 && points_remote<npts) {
        cellstatus = 2;
      }
      // all on remote processes, but not all the same process
      else if (process_count>1 && points_remote==npts) {
        cellstatus = 4;
      }
      else {
        throw std::string("This should not be possible");
      }
    }
    //
    // We will use the process receiving the first point as the desired final location :
    // ideally we'd count them and choose the process with the most, but the gain would be tiny
    //
    vtkIdType destProcess = localId_to_process_map[pts[0]];
    if (cellstatus>=3) {
      // The cell is going to be sent away, so add it to our send list
      cell_partitioninfo.Procs.push_back(destProcess); 
      cell_partitioninfo.GlobalIds.push_back(cellId + this->ZoltanCallbackData.ProcessOffsetsCellId[this->UpdatePiece]);
    }

    //
    // cells of type 2 and 4 need special treatment to handle the missing points
    //
    if (cellstatus==2 || cellstatus==4) {
      for (int i=0; i<npts; i++) {
        // if the point has been sent away - we need it back again
        if (cellstatus==2 && this->ZoltanCallbackData.LocalToLocalIdMap[pts[i]]==-1) {
          // add the point to the new output, but check for space first
          if (OutputPointCount>=ReservedPointCount) {
            // allow a small amount (try 5%) extra for duplicated boundary cell points
            // but if N is very small, double it
            vtkIdType increment = OutputPointCount/20;
            if (increment==0) increment = OutputPointCount;
            ReservedPointCount = OutputPointCount + increment;
            this->ZoltanCallbackData.Output->GetPoints()->GetData()->Resize(ReservedPointCount);
            this->ZoltanCallbackData.OutputPointsData = this->ZoltanCallbackData.Output->GetPoints()->GetData()->GetVoidPointer(0);
          }
          this->ZoltanCallbackData.Output->GetPointData()->CopyData(this->ZoltanCallbackData.Input->GetPointData(), pts[i], OutputPointCount);
          this->ZoltanCallbackData.LocalToLocalIdMap[pts[i]] = OutputPointCount;
          memcpy(&((T*)(this->ZoltanCallbackData.OutputPointsData))[OutputPointCount*3], &((T*)(this->ZoltanCallbackData.InputPointsData))[pts[i]*3], sizeof(T)*3);
          OutputPointCount++;
        }
        // if the cell was split over multiple remote processes, we must send the points
        // needed to complete the cell to the correct process
        if (cellstatus==4 && localId_to_process_map[pts[i]] != destProcess) {
//          std::cout << ", " << pts[i] << " will be sent to " << destProcess << " instead of " << localId_to_process_map[pts[i]] << std::endl;
          // The point is going to be sent away, so add it to our send list
          process_vector.push_back( process_tuple(pts[i], destProcess) );
        }
      }
    }
  }
  // make sure final point count is actual point count, not the reserved amount
  this->ZoltanCallbackData.Output->GetPoints()->SetNumberOfPoints(OutputPointCount);
    
  std::stringstream temp1, temp2;

  std::sort(process_vector.begin(), process_vector.end());
  process_vector.resize(std::unique(process_vector.begin(), process_vector.end()) - process_vector.begin());

  point_partitioninfo.GlobalIds.reserve(process_vector.size());
  point_partitioninfo.Procs.reserve(process_vector.size());
  for (std::vector<process_tuple>::iterator x=process_vector.begin(); x!=process_vector.end(); ++x) 
  {
    point_partitioninfo.GlobalIds.push_back(x->first + this->ZoltanCallbackData.ProcessOffsetsPointId[this->UpdatePiece]);
    point_partitioninfo.Procs.push_back(x->second);
  }
//  copy(point_partitioninfo.LocalIds.begin(), point_partitioninfo.LocalIds.end(), std::ostream_iterator<vtkIdType>(temp2,", ") );
//  std::cout << "Sorted " << point_partitioninfo.LocalIds.size() << std::endl;
//  std::cout << temp2.str() << std::endl;

  return;
}
//----------------------------------------------------------------------------
