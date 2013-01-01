/*=========================================================================

  Module                  : vtkZoltanV1PartitionFilter.h
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
#include "vtkZoltanV1PartitionFilter.h"
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
#include "vtkZoltanV1PartitionFilter.h"
#include "vtkDummyController.h"
//
#include "vtkBoundsExtentTranslator.h"
//
#include <sstream>
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>
#include <algorithm>
#include <map>
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
vtkStandardNewMacro(vtkZoltanV1PartitionFilter);
vtkCxxSetObjectMacro(vtkZoltanV1PartitionFilter, Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
#ifdef EXTRA_ZOLTAN_DEBUG
  static int pack_count   = 0;
  static int unpack_count = 0;
  static int size_count   = 0;
  #define INC_PACK_COUNT pack_count++;
  #define INC_UNPACK_COUNT unpack_count++;
  #define INC_SIZE_COUNT size_count++;
  #define CLEAR_ZOLTAN_DEBUG pack_count = 0; size_count = 0; unpack_count = 0;
#else
  #define INC_PACK_COUNT 
  #define INC_UNPACK_COUNT 
  #define INC_SIZE_COUNT 
  #define CLEAR_ZOLTAN_DEBUG 
#endif
//
typedef vtkZoltanV1PartitionFilter::ProcessExchangeVariables Exchange;
//----------------------------------------------------------------------------
// Zoltan callback interface
//----------------------------------------------------------------------------
void add_Id_to_interval_map(Exchange *mesh, vtkIdType GID, vtkIdType LID) {
  vtkIdType diff = GID-LID;
  std::map<vtkIdType,vtkIdType>::reverse_iterator prev = mesh->ReceivedGlobalToLocalIdMap.rbegin();
  if (prev!=mesh->ReceivedGlobalToLocalIdMap.rend()) {
    vtkIdType key = prev->first;
    vtkIdType delta = prev->second;
    if (diff==delta) {
      // no need to add this value to the map as it has the same delta as the previous entry
    }
    else {
      // add a new entry to the 'interval' map
      mesh->ReceivedGlobalToLocalIdMap[GID] = diff;
    }
  }
  else {
    // add a new entry to the 'interval' map
    mesh->ReceivedGlobalToLocalIdMap[GID] = diff;
  }
}
//----------------------------------------------------------------------------
vtkIdType global_to_local_Id(Exchange *mesh, vtkIdType GID) {
  std::map<vtkIdType,vtkIdType>::iterator ub = mesh->ReceivedGlobalToLocalIdMap.upper_bound(GID);
    ub--;
    vtkIdType key = ub->first;
    vtkIdType delta = ub->second;
    return GID - delta;
}
//----------------------------------------------------------------------------
// Zoltan callback which returns number of objects participating in exchange
//----------------------------------------------------------------------------
static int get_number_of_objects_points(void *data, int *ierr)
{
  int res = static_cast<Exchange*>(data)->InputNumberOfLocalPoints; 
  *ierr = (res < 0) ? ZOLTAN_FATAL : ZOLTAN_OK;
  return res;
}
//----------------------------------------------------------------------------
// Zoltan callback which fills the Ids for each object in the exchange
//----------------------------------------------------------------------------
static void get_object_list_points(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr)
{
  Exchange *mesh = (Exchange*)data;
  //
  // Return the IDs of our objects, but no weights.
  // Zoltan will assume equally weighted objects.
  //
  for (int i=0; i<mesh->InputNumberOfLocalPoints; i++){
    globalID[i] = i + mesh->ProcessOffsetsPointId[mesh->ProcessRank];
    localID[i] = i;
  }
  *ierr = ZOLTAN_OK;
}
//----------------------------------------------------------------------------
// Zoltan callback which returns the dimension of geometry (3D for us)
//----------------------------------------------------------------------------
static int get_num_geometry(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 3;
}
//----------------------------------------------------------------------------
// Zoltan callback which returns coordinate geometry data (points)
// templated here to alow float/double instances in our implementation
//----------------------------------------------------------------------------
template<typename T>
void get_geometry_list(
  void *data, int sizeGID, int sizeLID, int num_obj, 
  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
  int num_dim, double *geom_vec, int *ierr)
{
  Exchange *mesh = (Exchange*)data;
  for (int i=0;  i<num_obj; i++){
    geom_vec[3*i]   = ((T*)(mesh->InputPointsData))[3*i+0];
    geom_vec[3*i+1] = ((T*)(mesh->InputPointsData))[3*i+1];
    geom_vec[3*i+2] = ((T*)(mesh->InputPointsData))[3*i+2];
  }
  *ierr = ZOLTAN_OK;
  return;
}
//----------------------------------------------------------------------------
// A ZOLTAN_OBJ_SIZE_FN query function returns the size (in bytes) of the data buffer 
// that is needed to pack all of a single object's data.
//
// Here we add up the size of all the field arrays for points + the geometry itself
//  
// Function Type:   ZOLTAN_OBJ_SIZE_FN_TYPE
// Arguments:   
//  data             Pointer to user-defined data.
//  num_gid_entries  The number of array entries used to describe a single global ID.  
//  num_lid_entries  The number of array entries used to describe a single local ID.  
//  global_id        Pointer to the global ID of the object.
//  local_id         Pointer to the local ID of the object.
//  ierr             Error code to be set by function.
// Returned Value:   
//  int              The size (in bytes) of the required data buffer (per object).
//----------------------------------------------------------------------------
template<typename T>
int zoltan_obj_size_func_points(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
{
  INC_SIZE_COUNT
  Exchange *mesh = (Exchange*)data;
  *ierr = ZOLTAN_OK;
  return mesh->TotalSizePerId + sizeof(T)*3;
}
//----------------------------------------------------------------------------
// Zoltan callback to pack all the data for one point into a buffer
//----------------------------------------------------------------------------
template<typename T>
void zoltan_pack_obj_func_points(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size, char *buf, int *ierr)
{
  INC_PACK_COUNT
  Exchange *mesh = (Exchange*)data;
  vtkIdType GID = *global_id;
  vtkIdType LID = *local_id;
  //
  for (int i=0; i<mesh->NumberOfPointFields; i++) {
    int asize = mesh->MemoryPerTuple[i];
    char *dataptr = (char*)(mesh->InputArrayPointers[i]) + asize*LID;
    memcpy(buf, dataptr, asize);
    buf += asize;
  }
  memcpy(buf, &((T*)(mesh->InputPointsData))[(*local_id)*3], sizeof(T)*3);  
  *ierr = ZOLTAN_OK;
  return;
}
//----------------------------------------------------------------------------
// Zoltan callback to unpack all the data for one point from a buffer
//----------------------------------------------------------------------------
template<typename T>
void zoltan_unpack_obj_func_points(void *data, int num_gid_entries,
  ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr)
{
  INC_UNPACK_COUNT
  if (num_gid_entries != 1) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  Exchange *mesh = (Exchange*)data;
  vtkIdType GID = *global_id;
  //
  vtkPointData *inPD  = mesh->Input->GetPointData();
  vtkPointData *outPD = mesh->Output->GetPointData();
  //
  for (int i=0; i<mesh->NumberOfPointFields; i++) {
    int asize = mesh->MemoryPerTuple[i];
    char *dataptr = (char*)(mesh->OutputArrayPointers[i]) + asize*(mesh->OutPointCount);
    memcpy(dataptr, buf, asize);
    buf += asize;
  }
  add_Id_to_interval_map(mesh, GID, mesh->OutPointCount);
  memcpy(&((T*)(mesh->OutputPointsData))[mesh->OutPointCount*3], buf, sizeof(T)*3);  
  mesh->OutPointCount++;
  *ierr = ZOLTAN_OK;
  return;
}
//----------------------------------------------------------------------------
// Function Type: Pre migration callback
// Arguments:   
//  data              Pointer to user-defined data.
//  num_gid_entries   The number of array entries used to describe a single global ID.  
//  num_lid_entries   The number of array entries used to describe a single local ID.  
//  num_import        The number of objects that will be received by this processor.
//  import_global_ids An array of num_import global IDs of objects to be received by this processor. 
//                    may be NULL, as the processor does not necessarily need to know which objects it will receive.
//  import_local_ids  An array of num_import local IDs of objects to be received by this processor. 
//                    may be NULL, as the processor does not necessarily need to know which objects it will receive.
//  import_procs      An array of size num_import listing the processor IDs of the source processors. 
//                    may be NULL, as the processor does not necessarily need to know which objects is will receive.
//  import_to_part    An array of size num_import listing the parts to which objects will be imported. 
//                    may be NULL, as the processor does not necessarily need to know which objects it will receive.
//  num_export        The number of objects that will be sent from this processor to other processors.
//  export_global_ids An array of num_export global IDs of objects to be sent from this processor.
//  export_local_ids  An array of num_export local IDs of objects to be sent from this processor.
//  export_procs      An array of size num_export listing the processor IDs of the destination processors.
//  export_to_part    An array of size num_export listing the parts to which objects will be sent.
//  ierr              Error code to be set by function.
template<typename T>
void zoltan_pre_migrate_func_points(void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  Exchange *mesh = (Exchange*)data;
  // newTotal = original points - sent away + received
  mesh->OutputNumberOfLocalPoints = mesh->InputNumberOfLocalPoints + num_import - num_export;
  mesh->OutputPoints->SetNumberOfPoints(mesh->OutputNumberOfLocalPoints);
  mesh->OutputPointsData = mesh->OutputPoints->GetData()->GetVoidPointer(0);
  vtkPointData    *inPD  = mesh->Input->GetPointData();
  vtkPointData    *outPD = mesh->Output->GetPointData();
  outPD->CopyAllocate(inPD, mesh->OutputNumberOfLocalPoints);

  // make sure pointer for copying point data are set to each data array start
  // and the size of each data array is recorded so we can quickly do a memcpy
  // for each ID that is transferred to/from this process
  mesh->TotalSizePerId = 0;
  for (int i=0; i<mesh->NumberOfPointFields; i++) {
    vtkDataArray *iarray = mesh->Input->GetPointData()->GetArray(i);
    vtkDataArray *oarray = mesh->Output->GetPointData()->GetArray(i);
    oarray->SetNumberOfTuples(mesh->OutputNumberOfLocalPoints);
    mesh->InputArrayPointers.push_back(iarray->GetVoidPointer(0));
    mesh->OutputArrayPointers.push_back(oarray->GetVoidPointer(0));
    // we need to know the amount of data to copy for each array tuple
    int Nc = iarray->GetNumberOfComponents();
    int Ns = iarray->GetDataTypeSize();
    mesh->MemoryPerTuple.push_back(Nc*Ns);
    mesh->TotalSizePerId += Nc*Ns;
  }
  // some points are being sent away, some will be received, we must copy
  // the ones that are not moving from the input to the output.
  // Mark points so we know which ones are (will be) local after sending away
  mesh->LocalToLocalIdMap.assign(mesh->InputNumberOfLocalPoints, 0);
  for (vtkIdType i=0; i<num_export; i++) mesh->LocalToLocalIdMap[export_local_ids[i]] = -1;    
  // Loop over each local point and copy it to the output.
  // WARNING: point Ids are changing so any cells referencing the points
  // must have their Ids updated to the new index - create an IdMap to hold this info.
  mesh->OutPointCount = 0;
  for (vtkIdType i=0; i<mesh->InputNumberOfLocalPoints; i++) {
    if (mesh->LocalToLocalIdMap[i]==0) {
      outPD->CopyData(inPD, i, mesh->OutPointCount);
      memcpy(&((T*)(mesh->OutputPointsData))[mesh->OutPointCount*3], &((T*)(mesh->InputPointsData))[i*3], sizeof(T)*3);
      mesh->LocalToLocalIdMap[i] = mesh->OutPointCount;
      mesh->OutPointCount++;
    }
  }
}
//----------------------------------------------------------------------------
// Function Type: Pre migration callback for halo particle exchange
// the halo migration is controlled manually. We declare what we need to exchange
// and use this function to ensure all memory is allocated correctly for the receive
// to take place and unpack each object
//----------------------------------------------------------------------------
template <typename T>
void zoltan_pre_migrate_func_halo(void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  Exchange *mesh = (Exchange*)data;
  // resize points to accept ghost cell additions
  mesh->OutputNumberOfPointsWithHalo = mesh->OutputNumberOfLocalPoints + num_import;
  mesh->OutputPoints->GetData()->Resize(mesh->OutputNumberOfPointsWithHalo);
  mesh->OutputPoints->SetNumberOfPoints(mesh->OutputNumberOfPointsWithHalo);
  mesh->OutputPointsData = (T*)(mesh->OutputPoints->GetData()->GetVoidPointer(0));
  // copies are now being made from existing data, so use output as new input points
  mesh->InputPointsData = mesh->OutputPointsData;
  vtkPointData   *inPD  = mesh->Output->GetPointData();
  vtkPointData   *outPD = mesh->Output->GetPointData();
  //
  // we must resize all the scalar/vector fields for the point data
  // WARNING: because this routine takes places after the auto-migration of points during
  // the zoltan load-balance, when copying point data, we set the input pointer to the output pointer 
  // because the output pointer is in fact the 'new' input
  mesh->InputArrayPointers.clear();
  mesh->OutputArrayPointers.clear();
  for (int i=0; i<mesh->NumberOfPointFields; i++) {
    vtkDataArray *oarray = mesh->Output->GetPointData()->GetArray(i);
    oarray->Resize(mesh->OutputNumberOfPointsWithHalo);
    oarray->SetNumberOfTuples(mesh->OutputNumberOfPointsWithHalo);
    mesh->InputArrayPointers.push_back(oarray->GetVoidPointer(0));
    mesh->OutputArrayPointers.push_back(oarray->GetVoidPointer(0));
  }
}
//----------------------------------------------------------------------------
// Function Type: Pre migration callback for cell exchange
// the cell exchange is controlled manually. We declare what we need to exchange
// and use this function to ensure all memory is allocated correctly for the receive
// to take place and unpack each object
//----------------------------------------------------------------------------
template <typename T>
void zoltan_pre_migrate_func_cell(void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  Exchange *mesh = (Exchange*)data;
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
  for (int i=0; i<mesh->NumberOfCellFields; i++) {
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
int zoltan_obj_size_func_cell(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
{
  INC_SIZE_COUNT
  Exchange *mesh = (Exchange*)data;
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
void zoltan_pack_obj_func_cell(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size, char *buf, int *ierr)
{
  INC_PACK_COUNT
  Exchange *mesh = (Exchange*)data;
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
  for (int i=0; i<mesh->NumberOfCellFields; i++) {
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
void zoltan_unpack_obj_func_cell(void *data, int num_gid_entries,
  ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr)
{
  INC_UNPACK_COUNT
  Exchange *mesh = (Exchange*)data;
  vtkIdType GID = *global_id;
  //
  vtkIdType npts, *pts, ctype, newPts[32];
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(mesh->Input);
  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(mesh->Input);
  vtkPolyData         *pdata2 = vtkPolyData::SafeDownCast(mesh->Output);
  vtkUnstructuredGrid *udata2 = vtkUnstructuredGrid::SafeDownCast(mesh->Output);
  //
  for (int i=0; i<mesh->NumberOfCellFields; i++) {
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
    newPts[i] = global_to_local_Id(mesh, pts[i]);
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
// vtkZoltanV1PartitionFilter :: implementation 
//----------------------------------------------------------------------------
vtkZoltanV1PartitionFilter::vtkZoltanV1PartitionFilter()
{
  this->UpdatePiece               = 0;
  this->UpdateNumPieces           = 1;
  this->IdChannelArray            = NULL;
  this->GhostCellOverlap          = 0.0;
  this->MaxAspectRatio            = 5.0;
  this->ExtentTranslator          = vtkBoundsExtentTranslator::New();
  this->ExchangePoints            = 1;
  this->ExchangeHaloPoints        = 0;
  this->ExchangeCells             = 1;
  this->Controller                = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
}

//----------------------------------------------------------------------------
vtkZoltanV1PartitionFilter::~vtkZoltanV1PartitionFilter()
{
  this->SetController(NULL);
  this->SetIdChannelArray(NULL);
  this->ExtentTranslator->Delete();
}

//----------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // supports any vtkPointSet type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
  return 1;
}
//----------------------------------------------------------------------------
vtkBoundingBox *vtkZoltanV1PartitionFilter::GetPartitionBoundingBox(int partition)
{
  if (partition<this->BoxList.size()) {
    return &this->BoxList[partition];
  }
  vtkErrorMacro(<<"Partition not found in Bounding Box list");
  return NULL;
}
//----------------------------------------------------------------------------
vtkBoundingBox *vtkZoltanV1PartitionFilter::GetPartitionBoundingBoxWithHalo(int partition)
{
  if (partition<this->BoxListWithHalo.size()) {
    return &this->BoxListWithHalo[partition];
  }
  vtkErrorMacro(<<"Partition not found in Bounding Box list");
  return NULL;
}
//----------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::FindPointsInHaloRegions(vtkPoints *pts, vtkIdTypeArray *IdArray, PartitionInfo &ghostinfo)
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
#if 0 
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
    for (std::vector<boxproc>::iterator it=good_boxes.begin(); 
      it!=good_boxes.end(); ++it) 
    {
      vtkBoundingBox *b = (it->first);
      if (b->ContainsPoint(pt)) {
        ghostinfo.GlobalIds.push_back(IdArray->GetValue(i));
        ghostinfo.LocalIds.push_back(i);
        ghostinfo.Procs.push_back(it->second);
      }
    }
  }
//  std::cout << "There are " << ghostinfo.GlobalIds.size() << " ghosts on the list" <<std::endl;
//  std::ostream_iterator<vtkIdType> out_it (cout,", ");
//  copy ( ghostinfo.GlobalIds.begin(), ghostinfo.GlobalIds.end(), out_it );
}
//----------------------------------------------------------------------------
vtkSmartPointer<vtkIdTypeArray> vtkZoltanV1PartitionFilter::GenerateGlobalIds(vtkIdType Npoints, vtkIdType Ncells, const char *ptidname, Exchange *mesh)
{
  vtkSmartPointer<vtkIdTypeArray> ptIds = vtkSmartPointer<vtkIdTypeArray>::New();

  std::vector<vtkIdType> PointsPerProcess(this->UpdateNumPieces);
  std::vector<vtkIdType> CellsPerProcess(this->UpdateNumPieces);
  mesh->ProcessOffsetsPointId.assign(this->UpdateNumPieces+1, 0);
  mesh->ProcessOffsetsCellId.assign(this->UpdateNumPieces+1, 0);

  // create array of Id offsets increasing for each process rank
#ifdef VTK_USE_MPI
  vtkMPICommunicator* com = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  com->AllGather(&Npoints, &PointsPerProcess[0], 1);
  com->AllGather(&Ncells,  &CellsPerProcess[0], 1);
  std::partial_sum(PointsPerProcess.begin(), PointsPerProcess.end(), mesh->ProcessOffsetsPointId.begin()+1);
  std::partial_sum(CellsPerProcess.begin(), CellsPerProcess.end(), mesh->ProcessOffsetsCellId.begin()+1);
#endif

  //
  // Global point IDs generated here
  //
  vtkIdType offset = mesh->ProcessOffsetsPointId[this->UpdatePiece];
  ptIds->SetNumberOfValues(Npoints);
  for (vtkIdType id=0; id<Npoints; id++) {
    ptIds->SetValue(id, id + offset);
  }
  ptIds->SetName(ptidname);
  return ptIds;
}
//----------------------------------------------------------------------------
struct vtkPPF_datainfo {
  int  datatype;
  int  numC;
  char name[64];
  vtkPPF_datainfo() : datatype(-1), numC(-1) {};
};
//----------------------------------------------------------------------------
bool vtkZoltanV1PartitionFilter::GatherDataArrayInfo(vtkDataArray *data, 
  int &datatype, std::string &dataname, int &numComponents)
{
#ifdef VTK_USE_MPI
  std::vector< vtkPPF_datainfo > datatypes(this->UpdateNumPieces);
  if (data) {
    ((vtkPPF_datainfo*)&datatypes[this->UpdatePiece])->datatype = data->GetDataType();
    ((vtkPPF_datainfo*)&datatypes[this->UpdatePiece])->numC     = data->GetNumberOfComponents();
    strncpy(((vtkPPF_datainfo*)&datatypes[this->UpdatePiece])->name, data->GetName(), 64);
  }
  vtkMPICommunicator* com = vtkMPICommunicator::SafeDownCast(
    this->Controller->GetCommunicator()); 
  int result = com->AllGather((char*)MPI_IN_PLACE, (char*)&datatypes[0], sizeof(vtkPPF_datainfo));
  for (int i=0; i<this->UpdateNumPieces; i++) {
    vtkPPF_datainfo &newdata = datatypes[i];
    if (newdata.datatype!=-1) {
      datatype = newdata.datatype;
      numComponents = newdata.numC;
      dataname = newdata.name;
    }
  }
  return (result == 1) ;
#else
  return 1;
#endif
}
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::GatherDataTypeInfo(vtkPoints *points)
{
#ifdef VTK_USE_MPI
  if (this->UpdateNumPieces==1) {
      return points->GetDataType();
  }
  std::vector< int > datatypes(this->UpdateNumPieces, -1);
  int datatype = -1;
  datatypes[this->UpdatePiece] = points->GetDataType();
  vtkMPICommunicator* com = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator()); 
  int result = com->AllGather((int*)MPI_IN_PLACE, (int*)&datatypes[0], 1);
  for (int i=0; i<this->UpdateNumPieces; i++) {
    int &newdatatype = datatypes[i];
    if (datatype==-1 && newdatatype!=-1) {
      datatype = newdatatype;
    }
    else if (datatype!=-1 && newdatatype!=-1 && newdatatype!=datatype) {
      vtkErrorMacro(<<"Fatal datatype error in Point DataType Gather");
    }
  }
  return datatype;
#else
  return pInput->GetPoints()->GetDataType();
#endif
}
//-------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::InitBoundingBoxes(vtkDataSet *input, vtkBoundingBox &box) 
{
  double bounds[6];
  input->GetBounds(bounds);
  this->BoxList.clear();
  this->BoxListWithHalo.clear();
  if (this->UpdateNumPieces==1) {
    vtkBoundingBox databox(bounds);
    this->BoxList.push_back(databox);
    // we add a ghost cell region to our boxes
    databox.Inflate(this->GhostCellOverlap);
    this->BoxListWithHalo.push_back(databox);
    this->ExtentTranslator->SetNumberOfPieces(1);
    // Copy the bounds to our piece to bounds translator
    this->ExtentTranslator->SetBoundsForPiece(0, bounds);
    this->ExtentTranslator->InitWholeBounds();
  } 
  else {
    double bmin[3] = {bounds[0], bounds[2], bounds[4]};
    double bmax[3] = {bounds[1], bounds[3], bounds[5]};
    if (!vtkMath::AreBoundsInitialized(bounds)) {
      bmin[0] = bmin[1] = bmin[2] = VTK_DOUBLE_MAX;
      bmax[0] = bmax[1] = bmax[2] = VTK_DOUBLE_MIN;
    }
    if (this->Controller) {
      double globalMins[3], globalMaxes[3];
      this->Controller->AllReduce(bmin, globalMins,  3, vtkCommunicator::MIN_OP);
      this->Controller->AllReduce(bmax, globalMaxes, 3, vtkCommunicator::MAX_OP);
      bmin[0] = globalMins[0];  bmax[0] = globalMaxes[0];
      bmin[1] = globalMins[1];  bmax[1] = globalMaxes[1];
      bmin[2] = globalMins[2];  bmax[2] = globalMaxes[2];
    }
    box.SetMinPoint(bmin);
    box.SetMaxPoint(bmin);
  }
}
//-------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::SetupFieldArrayPointers(vtkDataSetAttributes *fields, Exchange &mesh, int &NumberOfFields) 
{
  int NumberOfFieldArrays = fields->GetNumberOfArrays();
  this->Controller->AllReduce(&NumberOfFieldArrays, &NumberOfFields, 1, vtkCommunicator::MAX_OP);
  for (int i=0; i<NumberOfFields; i++) {
    vtkSmartPointer<vtkDataArray> darray = fields->GetArray(i);
    //
    int correctType = -1, numComponents = -1;
    std::string correctName;
    this->GatherDataArrayInfo(darray, correctType, correctName, numComponents);
    if (!darray) {
      vtkDebugMacro(<<"NULL data found, used MPI_Gather to find :" 
        << " DataType " << correctType
        << " Name " << correctName.c_str()
        << " NumComponents " << numComponents);
      darray.TakeReference(vtkDataArray::CreateDataArray(correctType));
      darray->SetNumberOfComponents(numComponents);
      darray->SetName(correctName.c_str());
      fields->AddArray(darray);
    }
  }
}
//-------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  int piece, numPieces;
  piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), piece);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), numPieces);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);

  return 1;
}
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::RequestInformation(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
#ifdef VTK_USE_MPI
  vtkMPICommunicator *communicator = 
    vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  int maxpieces = communicator ? communicator->GetNumberOfProcesses() : 1;
#else
  int maxpieces = 1;
#endif

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  outInfo->Set(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR(), this->ExtentTranslator);
  //
//  outInfo->Set(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR(),
//               inInfo->Get(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR()));
  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  //
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()),
               6);
  return 1;
}
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  CLEAR_ZOLTAN_DEBUG
  //
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
// 
// Build a list of cell to process Ids. We will use this to send cells to their
// destination. We must also build a list of cells which are split across processes
// this is because some of the cell vertices will already have been redistributed
// during the initial partition, but we must also send the points which are on this process
// and part of the cell which ...hold one, we can keep partial cells ...
//
vtkSmartPointer<vtkIntArray> vtkZoltanV1PartitionFilter::BuildCellToProcessList(     
  vtkDataSet *data, PartitionInfo &partitioninfo, std::vector<int> &ProcessOffsetsCellId,
  int numExport,
  ZOLTAN_ID_PTR exportGlobalGids,
  ZOLTAN_ID_PTR exportLocalGids,
  int *exportProcs)
{
  vtkIdType numPts = data->GetNumberOfPoints();
  vtkIdType numCells = data->GetNumberOfCells();
  vtkIdType j, cellId;
  // we will put the results of the cell tests in these arrays
  partitioninfo.Procs.reserve(numCells);
  partitioninfo.LocalIds.reserve(numCells);
  partitioninfo.GlobalIds.reserve(numCells);

  // an optional flag for boundary cells we might use in the future
  vtkSmartPointer<vtkIntArray> boundaryFlag = vtkSmartPointer<vtkIntArray>::New();
  boundaryFlag->SetNumberOfTuples(numCells);
  boundaryFlag->SetName("Boundary");
  int *boundaryData = boundaryFlag->GetPointer(0);

  // initialize point to process array with this process Id
  std::vector<int> processdata(numPts, this->UpdatePiece); 
  // for each (local) point that we know is being exported, flag the point in the array
  for (vtkIdType i=0; i<numExport; i++) {
    vtkIdType id = exportLocalGids[i];
    processdata[id] = exportProcs[i];
  }

  // identify cells we need to export
  vtkIdType npts, *pts;
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(data);
  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(data);
  if (pdata || udata) {
    // polydata requires a cell map (verts/lines/polys/strips) to be present before we traverse cells
    if (pdata) pdata->BuildCells();
    // for each cell, find if all points are designated as remote and cell needs to be sent away
    for (cellId=0; cellId<numCells; cellId++) {
      // get a pointer to the cell points
      if (pdata) { pdata->GetCellPoints(cellId, npts, pts); }
      else if (udata) { udata->GetCellPoints(cellId, npts, pts); }
      // if all the points are remote, we will send this cell away
      vtkIdType destProcess = processdata[pts[0]];
      // Optional : if any of the points are on a different process, it's a boundary cell
      for (j=1; j<npts; j++) {
        if (processdata[pts[j]]!=destProcess) {
          boundaryData[cellId] = 1;
        }
      }
      // Check if all points are remote, partial cells will be retained locally
      bool remote = (destProcess!=this->UpdatePiece);
      for (j=1; remote && j<npts; j++) {
        destProcess = processdata[pts[j]];
        remote = remote && (destProcess!=this->UpdatePiece);
      }
      if (remote) {
        partitioninfo.Procs.push_back(destProcess); 
        partitioninfo.LocalIds.push_back(cellId);
        partitioninfo.GlobalIds.push_back(cellId + ProcessOffsetsCellId[this->UpdatePiece]);
      }
    }
  }
  else {
    vtkErrorMacro(<<"Only Polydata and UnstructuredGrid supported so far");  
  }
  return boundaryFlag;
}
//----------------------------------------------------------------------------
