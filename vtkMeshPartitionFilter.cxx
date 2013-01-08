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
void vtkMeshPartitionFilter::zoltan_pre_migrate_func_cell(void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  CallbackData *mesh = (CallbackData*)data;
  // How many cells will we have at the end
  mesh->OutputNumberOfLocalCells = mesh->InputNumberOfLocalCells + num_import - num_export;
  vtkIdType ReservedPointCount = mesh->OutputNumberOfLocalPoints;
  mesh->OutputPointsData = mesh->OutputPoints->GetData()->GetVoidPointer(0);
  // we will assume all cells are the same type (triangle/quad), simplify allocation/access
  mesh->MaxCellSize = mesh->Input->GetMaxCellSize();
  mesh->OutputCellArray = vtkSmartPointer<vtkCellArray>::New();
  mesh->OutputCellArray->Allocate(mesh->OutputNumberOfLocalCells*(mesh->MaxCellSize+1));

  vtkCellData *inCD  = mesh->Input->GetCellData();
  vtkCellData *outCD = mesh->Output->GetCellData();
  outCD->CopyAllocate(inCD, mesh->OutputNumberOfLocalCells);

  // make sure pointers for copying cell data are set to each data array start
  // and the size of each data array is recorded so we can quickly do a memcpy
  // for each ID that is transferred to/from this process
  mesh->TotalSizePerId = 0;
  mesh->InputArrayPointers.clear();
  mesh->OutputArrayPointers.clear();
  for (int i=0; i<mesh->NumberOfFields; i++) {
    vtkDataArray *iarray = inCD->GetArray(i);
    vtkDataArray *oarray = outCD->GetArray(i);
    oarray->Resize(mesh->OutputNumberOfLocalCells);
    oarray->SetNumberOfTuples(mesh->OutputNumberOfLocalCells);
    mesh->InputArrayPointers.push_back(iarray->GetVoidPointer(0));
    mesh->OutputArrayPointers.push_back(oarray->GetVoidPointer(0));
    // we need to know the amount of data to copy for each array tuple
    int Nc = iarray->GetNumberOfComponents();
    int Ns = iarray->GetDataTypeSize();
    mesh->MemoryPerTuple.push_back(Nc*Ns);
    mesh->TotalSizePerId += Nc*Ns;
  }

  std::vector<bool> local(mesh->InputNumberOfLocalCells, true);
  for (vtkIdType i=0; i<num_export; i++) {
    local[export_local_ids[i]] = false;    
  }
  //
  vtkIdType npts, *pts, newpts[32];
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(mesh->Input);
  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(mesh->Input);
  vtkPolyData         *pdata2 = vtkPolyData::SafeDownCast(mesh->Output);
  vtkUnstructuredGrid *udata2 = vtkUnstructuredGrid::SafeDownCast(mesh->Output);
  if (pdata2) {
    pdata2->SetPolys(mesh->OutputCellArray);
  }
  //
  mesh->OutCellCount = 0;
  for (vtkIdType cellId=0; cellId<mesh->InputNumberOfLocalCells; cellId++) {
    if (local[cellId]) {
      // copy cell data from old to new datasets
      outCD->CopyData(inCD, cellId, mesh->OutCellCount);
      // copy cell point Ids to new dataset, 
      if (pdata) { 
        int ctype = pdata->GetCellType(cellId);
        pdata->GetCellPoints(cellId, npts, pts); 
        for (int i=0; i<npts; i++) {
          if (mesh->LocalToLocalIdMap[pts[i]]!=-1) newpts[i] = mesh->LocalToLocalIdMap[pts[i]];
          else {
            // a point which has been migrated has been referenced. 
            // It can only happen on a boundary cell which has had some points sent away
            // option 1. Make a local copy of the point we've already sent, 
            // option 2. Send the other points too and make the cell remote. 
            // easiest solution is 1 currently.
            if (mesh->OutPointCount>=ReservedPointCount) {
              // output points should have been resized to fit all existing and migrated points, 
              // but we need to allow a small amount (try 5%) extra for duplicated boundary cell points
              ReservedPointCount = mesh->OutPointCount * 1.05;
              mesh->OutputPoints->GetData()->Resize(ReservedPointCount);
              mesh->OutputPointsData = mesh->OutputPoints->GetData()->GetVoidPointer(0);
            }
            mesh->Output->GetPointData()->CopyData(mesh->Input->GetPointData(), pts[i], mesh->OutPointCount);
            newpts[i] = mesh->LocalToLocalIdMap[pts[i]] = mesh->OutPointCount;
            memcpy(&((T*)(mesh->OutputPointsData))[mesh->OutPointCount*3], &((T*)(mesh->InputPointsData))[pts[i]*3], sizeof(T)*3);
            mesh->OutPointCount++;
          }
        }
        pdata2->InsertNextCell(ctype, npts, newpts);
      }
      else if (udata) { 
        int ctype = udata->GetCellType(cellId);
        udata->GetCellPoints(cellId, npts, pts); 
        udata2->InsertNextCell(ctype, npts, pts);
      }
      mesh->OutCellCount++;
    }
  }
  // make sure final point count is actual point count, not the reserved amount
  mesh->OutputPoints->SetNumberOfPoints(mesh->OutPointCount);
  mesh->OutputNumberOfLocalPoints = mesh->OutPointCount;
}
//----------------------------------------------------------------------------
// Zoltan callback function : returns size of each cell and all its data
//----------------------------------------------------------------------------
int vtkMeshPartitionFilter::zoltan_obj_size_func_cell(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
{
  INC_SIZE_COUNT
  vtkZoltanV1PartitionFilter::CallbackData *mesh = (vtkZoltanV1PartitionFilter::CallbackData*)data;
  *ierr = ZOLTAN_OK;
  // return the size of the cell data + number of points in cell + point Ids
  vtkIdType LID = *local_id;
  //
  vtkIdType npts, *pts;
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(mesh->Input);
  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(mesh->Input);
  //
  int ctype = pdata->GetCellType(LID);
  pdata->GetCellPoints(LID,npts,pts);
  int size = 2 + npts; // npts + ctype + pts
  return mesh->TotalSizePerId + size*sizeof(vtkIdType);
}
//----------------------------------------------------------------------------
// Zoltan callback function : pack one cell and all its data
//----------------------------------------------------------------------------
void vtkMeshPartitionFilter::zoltan_pack_obj_func_cell(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size, char *buf, int *ierr)
{
  INC_PACK_COUNT
  vtkZoltanV1PartitionFilter::CallbackData *mesh = (vtkZoltanV1PartitionFilter::CallbackData*)data;
  vtkIdType GID = *global_id;
  vtkIdType LID = *local_id;
  //
  //
  vtkIdType npts, *pts, newPts[32];
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(mesh->Input);
  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(mesh->Input);
  //
  // Copy all Cell data arrays into the buffer provided by zoltan
  //
  for (int i=0; i<mesh->NumberOfFields; i++) {
    int asize = mesh->MemoryPerTuple[i];
    char *dataptr = (char*)(mesh->InputArrayPointers[i]) + asize*LID;
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
      newPts[i+2] = pts[i] + mesh->ProcessOffsetsPointId[mesh->ProcessRank];
    }
    memcpy(buf, newPts, sizeof(vtkIdType)*(npts+2));  
  }
  else if (udata) { 
    throw std::exception("Implement this");
  }

  *ierr = ZOLTAN_OK;
  return;
}
//----------------------------------------------------------------------------
// Zoltan callback function : unpack one cell and all its data
//----------------------------------------------------------------------------
void vtkMeshPartitionFilter::zoltan_unpack_obj_func_cell(void *data, int num_gid_entries,
  ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr)
{
  INC_UNPACK_COUNT
  vtkZoltanV1PartitionFilter::CallbackData *mesh = (vtkZoltanV1PartitionFilter::CallbackData*)data;
  //
  vtkIdType GID = *global_id;
  //
  vtkIdType npts, *pts, ctype, newPts[32];
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(mesh->Input);
  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(mesh->Input);
  vtkPolyData         *pdata2 = vtkPolyData::SafeDownCast(mesh->Output);
  vtkUnstructuredGrid *udata2 = vtkUnstructuredGrid::SafeDownCast(mesh->Output);
  //
  for (int i=0; i<mesh->NumberOfFields; i++) {
    int asize = mesh->MemoryPerTuple[i];
    char *dataptr = (char*)(mesh->InputArrayPointers[i]) + asize*(mesh->OutCellCount);
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
    newPts[i] = mesh->self->global_to_local_Id(pts[i]);
  }
  if (pdata) {
    pdata2->InsertNextCell(ctype, npts, newPts);
  }
  else if (udata) { 
    throw std::exception("Implement this");
  }

  mesh->OutCellCount++;
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

  //
  // Distribute cells based on the usage of the points already distributed
  //
  this->PartitionCells(info, inputVector, outputVector);

  //*****************************************************************
  // Free the arrays allocated by Zoltan_LB_Partition, and free
  // the storage allocated for the Zoltan structure.
  //*****************************************************************
  if (this->LoadBalanceData.importGlobalGids) {
    Zoltan_LB_Free_Part(&this->LoadBalanceData.importGlobalGids, &this->LoadBalanceData.importLocalGids, &this->LoadBalanceData.importProcs, &this->LoadBalanceData.importToPart);
    Zoltan_LB_Free_Part(&this->LoadBalanceData.exportGlobalGids, &this->LoadBalanceData.exportLocalGids, &this->LoadBalanceData.exportProcs, &this->LoadBalanceData.exportToPart);
  }
  Zoltan_Destroy(&this->ZoltanData);

  this->Controller->Barrier();
  this->Timer->StopTimer();
  vtkDebugMacro(<<"Particle partitioning : " << this->Timer->GetElapsedTime() << " seconds");
  return 1;
}
//----------------------------------------------------------------------------
int vtkMeshPartitionFilter::PartitionCells(vtkInformation* info,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  vtkSmartPointer<vtkIntArray> boundary;
  PartitionInfo partitioninfo;

  boundary = this->BuildCellToProcessList(this->ZoltanCallbackData.Input, 
    partitioninfo,
    this->ZoltanCallbackData.ProcessOffsetsCellId,
    this->LoadBalanceData.numExport,            // Number of vertices I must send to other processes
    this->LoadBalanceData.exportGlobalGids,     // Global IDs of the vertices I must send
    this->LoadBalanceData.exportLocalGids,      // Local IDs of the vertices I must send
    this->LoadBalanceData.exportProcs           // Process to which I send each of the vertices
    );

  this->SetupFieldArrayPointers(this->ZoltanCallbackData.Output->GetCellData());
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
  int zoltan_error = Zoltan_Invert_Lists(this->ZoltanData, 
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
  if (zoltan_error != ZOLTAN_OK){
    printf("Zoltan_LB_Partition NOT OK...\n");
    MPI_Finalize();
    Zoltan_Destroy(&this->ZoltanData);
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
  Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &this->ZoltanCallbackData); 
  Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &this->ZoltanCallbackData); 
  Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &this->ZoltanCallbackData); 
  Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &this->ZoltanCallbackData); 

  //
  // Perform the cell exchange
  //
  zoltan_error = Zoltan_Migrate (this->ZoltanData,
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

  return 1;
}

//----------------------------------------------------------------------------
