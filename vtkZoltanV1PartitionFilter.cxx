/*=========================================================================

  Module                  : vtkPartitionOutline.h

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
#include "vtkCubeSource.h"
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
#include "vtkNew.h"
#include "vtkHexahedron.h"
#include "vtkPKdTree.h"
#include "vtkBSPCuts.h"
#include "vtkKdTreeGenerator.h"
//
#include "vtkPKdTree2.h"
#include "vtkBoundsExtentTranslator.h"
#include "vtkZoltanV1PartitionFilter.h"
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>
#include <algorithm>
#include <map>
#include <iostream>
#include <ostream>
#include <sstream>
#include <iterator>
//
#include <stack>
#include "zz_const.h"
#include "rcb.h"
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkZoltanV1PartitionFilter);
vtkCxxSetObjectMacro(vtkZoltanV1PartitionFilter, Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
// Zoltan callback which returns number of objects participating in exchange
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::get_number_of_objects_points(void *data, int *ierr)
{
  int res = static_cast<CallbackData*>(data)->Input->GetNumberOfPoints(); 
  *ierr = (res < 0) ? ZOLTAN_FATAL : ZOLTAN_OK;
  return res;
}
//----------------------------------------------------------------------------
// Zoltan callback which fills the Ids for each object in the exchange
//----------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::get_object_list_points(void *data, int sizeGID, int sizeLID,
  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  //
  // Return the IDs of our objects, but no weights.
  // Zoltan will assume equally weighted objects.
  //
  vtkIdType N = callbackdata->Input->GetNumberOfPoints();
  for (int i=0; i<N; i++){
    globalID[i] = i + callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];
//    localID[i] = i;
  }
  *ierr = ZOLTAN_OK;
}
//----------------------------------------------------------------------------
// Zoltan callback which fills the Ids for each object in the exchange
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::get_first_object_points(
  void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
  int wdim, float *wgt, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);

  *ierr = ZOLTAN_OK;
  vtkIdType N = callbackdata->Input->GetNumberOfPoints();
  if (N>0) {
    global_id[0] = 0 + callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];
//    local_id[0] = 0;
    return 1;
  }
  return 0;
}

int vtkZoltanV1PartitionFilter::get_next_object_points(
  void * data, int num_gid_entries, int num_lid_entries, 
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, ZOLTAN_ID_PTR next_global_id, ZOLTAN_ID_PTR next_local_id, 
  int wgt_dim, float *next_obj_wgt, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  //
  vtkIdType GID = *global_id;
  vtkIdType LID = GID - callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];

  *ierr = ZOLTAN_OK;
  vtkIdType N = callbackdata->Input->GetNumberOfPoints();
  if (LID<N) {
    next_global_id[0] = LID + 1  + callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];
//    next_local_id[0] = LID + 1;
    return 1;
  }
  return 0;
}


//----------------------------------------------------------------------------
// Zoltan callback which returns the dimension of geometry (3D for us)
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::get_num_geometry(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 3;
}
//----------------------------------------------------------------------------
// Zoltan callback which returns coordinate geometry data (points)
// templated here to alow float/double instances in our implementation
//----------------------------------------------------------------------------
template<typename T>
void vtkZoltanV1PartitionFilter::get_geometry_list(
  void *data, int sizeGID, int sizeLID, int num_obj, 
  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
  int num_dim, double *geom_vec, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  for (int i=0;  i<num_obj; i++){
    geom_vec[3*i]   = ((T*)(callbackdata->InputPointsData))[3*i+0];
    geom_vec[3*i+1] = ((T*)(callbackdata->InputPointsData))[3*i+1];
    geom_vec[3*i+2] = ((T*)(callbackdata->InputPointsData))[3*i+2];
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
int vtkZoltanV1PartitionFilter::zoltan_obj_size_function_points(void *data, 
  int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, 
  ZOLTAN_ID_PTR local_id, int *ierr)
{
  INC_SIZE_COUNT
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  *ierr = ZOLTAN_OK;
  return callbackdata->TotalSizePerId + sizeof(T)*3;
}
//----------------------------------------------------------------------------
// Zoltan callback to pack all the data for one point into a buffer
//----------------------------------------------------------------------------
template<typename T>
void vtkZoltanV1PartitionFilter::zoltan_pack_obj_function_points(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size, char *buf, int *ierr)
{
  INC_PACK_COUNT
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  vtkIdType GID = *global_id;
  vtkIdType LID = GID - callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];
  //
  for (int i=0; i<callbackdata->NumberOfFields; i++) {
    int asize = callbackdata->MemoryPerTuple[i];
    char *dataptr = (char*)(callbackdata->InputArrayPointers[i]) + asize*LID;
    memcpy(buf, dataptr, asize);
    buf += asize;
  }
  memcpy(buf, &((T*)(callbackdata->InputPointsData))[LID*3], sizeof(T)*3);  
  *ierr = ZOLTAN_OK;
  return;
}
//----------------------------------------------------------------------------
// Zoltan callback to unpack all the data for one point from a buffer
//----------------------------------------------------------------------------
template<typename T>
void vtkZoltanV1PartitionFilter::zoltan_unpack_obj_function_points(void *data, int num_gid_entries,
  ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr)
{
  INC_UNPACK_COUNT
  if (num_gid_entries != 1) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  vtkIdType GID = *global_id;
  //
  vtkPointData *inPD  = callbackdata->Input->GetPointData();
  vtkPointData *outPD = callbackdata->Output->GetPointData();
  //
  for (int i=0; i<callbackdata->NumberOfFields; i++) {
    int asize = callbackdata->MemoryPerTuple[i];
    char *dataptr = (char*)(callbackdata->OutputArrayPointers[i]) + asize*(callbackdata->OutPointCount);
    memcpy(dataptr, buf, asize);
    buf += asize;
  }
//  if (callbackdata->self->UpdatePiece==2 && GID <100) { std::cout <<"Received " << GID << std::endl; }
  add_Id_to_interval_map(callbackdata, GID, callbackdata->OutPointCount);
  memcpy(&((T*)(callbackdata->OutputPointsData))[callbackdata->OutPointCount*3], buf, sizeof(T)*3);  
  callbackdata->OutPointCount++;
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
void vtkZoltanV1PartitionFilter::zoltan_pre_migrate_function_points(
  void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  // newTotal = original points - sent away + received
  vtkIdType N  = callbackdata->Input->GetNumberOfPoints();
  vtkIdType N2 = N + num_import - num_export;
  callbackdata->Output->GetPoints()->SetNumberOfPoints(N2);
  callbackdata->OutputPointsData = callbackdata->Output->GetPoints()->GetData()->GetVoidPointer(0);
  vtkPointData    *inPD  = callbackdata->Input->GetPointData();
  vtkPointData    *outPD = callbackdata->Output->GetPointData();
  outPD->CopyAllocate(inPD, N2);
  //
  callbackdata->self->InitializeFieldDataArrayPointers(callbackdata, inPD, outPD, N2);

  // some points are being sent away, some will be received, we must copy
  // the ones that are not moving from the input to the output.
  // Mark points so we know which local points will still be local after the exchange
  callbackdata->LocalToLocalIdMap.assign(N, 0);
  for (vtkIdType i=0; i<num_export; i++) {
    vtkIdType GID = export_global_ids[i];
    vtkIdType LID = GID - callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];    
    callbackdata->LocalToLocalIdMap[LID] = -1;    
  }

  // Loop over each local point and copy it to the output.
  // WARNING: point Ids are changing so any cells referencing the points
  // must have their Ids updated to the new index - create an IdMap to hold this info.
  callbackdata->OutPointCount = 0;
  for (vtkIdType i=0; i<N; i++) {
    if (callbackdata->LocalToLocalIdMap[i]==0) {
      outPD->CopyData(inPD, i, callbackdata->OutPointCount);
      memcpy(&((T*)(callbackdata->OutputPointsData))[callbackdata->OutPointCount*3], &((T*)(callbackdata->InputPointsData))[i*3], sizeof(T)*3);
      callbackdata->LocalToLocalIdMap[i] = callbackdata->OutPointCount;
      callbackdata->OutPointCount++;
    }
  }
}
//----------------------------------------------------------------------------
// Function Type: Pre migration callback for halo/other particle exchange
// The difference between this and the standard pre_migrate function for points 
// is that we only add new points and do not remove any.
//----------------------------------------------------------------------------
template <typename T>
void vtkZoltanV1PartitionFilter::zoltan_pre_migrate_function_points_add(
  void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);

  // this is the current actual size of the output point list
  vtkIdType OutputNumberOfPoints = callbackdata->Output->GetNumberOfPoints();
  callbackdata->OutPointCount = OutputNumberOfPoints;
  // resize points to accept ghost cell additions
  OutputNumberOfPoints = OutputNumberOfPoints + num_import;
  callbackdata->Output->GetPoints()->GetData()->Resize(OutputNumberOfPoints);
  callbackdata->Output->GetPoints()->SetNumberOfPoints(OutputNumberOfPoints);
  callbackdata->OutputPointsData = (T*)(callbackdata->Output->GetPoints()->GetData()->GetVoidPointer(0));
  // 
  // The migration taking place might be using the original points (standard send/receive)
  // or the already migrated points (like a halo exchange). When doing a halo exchange, we might
  // be resending points we have just received, which are in our output list and not our
  // input list, so change the pointers accordingly
  //
  vtkPointData *inPD;
  if (callbackdata->CopyFromOutput) {
    callbackdata->InputPointsData = callbackdata->OutputPointsData;
    inPD  = callbackdata->Output->GetPointData();
  }
  else {
    callbackdata->InputPointsData = callbackdata->Input->GetPoints() ? (T*)(callbackdata->Input->GetPoints()->GetData()->GetVoidPointer(0)) : NULL;
    inPD  = callbackdata->Input->GetPointData();
  }
  vtkPointData *outPD = callbackdata->Output->GetPointData();
  callbackdata->self->InitializeFieldDataArrayPointers(callbackdata, inPD, outPD, OutputNumberOfPoints);
}
//----------------------------------------------------------------------------
// vtkZoltanV1PartitionFilter :: implementation 
//----------------------------------------------------------------------------
vtkZoltanV1PartitionFilter::vtkZoltanV1PartitionFilter()
{
  this->UpdatePiece               = 0;
  this->UpdateNumPieces           = 1;
  this->IdChannelArray            = NULL;
  this->MaxAspectRatio            = 5.0;
  this->ExtentTranslator          = vtkSmartPointer<vtkBoundsExtentTranslator>::New();
  this->InputExtentTranslator     = NULL;
  this->ZoltanData                = NULL;
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
}
//----------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter supports any vtkPointSet type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // output should be the same as the input (a subclass of vtkPointSet)
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
//template<typename T>
//std::ostream& PrintVector(std::ostream& out, int width, std::vector<T> &vec) {
//  out << "[" << std::setprecision(1) << std::fixed;
//  for (auto & x : vec) out << std::setw(width) << x << ", ";
//  return out << "]";
//}
//----------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::SetupGlobalIds(vtkPointSet *ps)
{
  this->ComputeIdOffsets(ps->GetNumberOfPoints(), ps->GetNumberOfCells());
}
//----------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::ComputeIdOffsets(vtkIdType Npoints, vtkIdType Ncells)
{
  // offset arrays
  std::vector<vtkIdType> PointsPerProcess(this->UpdateNumPieces);
  std::vector<vtkIdType> CellsPerProcess(this->UpdateNumPieces);
  this->ZoltanCallbackData.ProcessOffsetsPointId.assign(this->UpdateNumPieces+1, 0);
  this->ZoltanCallbackData.ProcessOffsetsCellId.assign(this->UpdateNumPieces+1, 0);
  // create array of Id offsets increasing for each process rank
#ifdef VTK_USE_MPI
  vtkMPICommunicator* com = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  com->AllGather(&Npoints, &PointsPerProcess[0], 1);
  com->AllGather(&Ncells,  &CellsPerProcess[0], 1);
  std::partial_sum(PointsPerProcess.begin(), PointsPerProcess.end(), this->ZoltanCallbackData.ProcessOffsetsPointId.begin()+1);
  std::partial_sum(CellsPerProcess.begin(), CellsPerProcess.end(), this->ZoltanCallbackData.ProcessOffsetsCellId.begin()+1);
#endif
//  std::stringstream temp1, temp2;
//  copy(this->ZoltanCallbackData.ProcessOffsetsPointId.begin(), this->ZoltanCallbackData.ProcessOffsetsPointId.end(), std::ostream_iterator<vtkIdType>(temp1,", ") );
//  copy(this->ZoltanCallbackData.ProcessOffsetsCellId.begin(), this->ZoltanCallbackData.ProcessOffsetsCellId.end(), std::ostream_iterator<vtkIdType>(temp2,", ") );
//  PrintVector<int>(temp1, 6, this->ZoltanCallbackData.ProcessOffsetsPointId);
//  PrintVector<int>(temp2, 6, this->ZoltanCallbackData.ProcessOffsetsCellId);
//  vtkDebugMacro(<< "Offsets generated { pts : " << temp1.str() << "} { cells : " << temp2.str() << "}" );
}
//----------------------------------------------------------------------------
struct vtkZPF_datainfo {
  int  datatype;
  int  numC;
  char name[64];
  vtkZPF_datainfo() : datatype(-1), numC(-1) {};
};
//----------------------------------------------------------------------------
bool vtkZoltanV1PartitionFilter::GatherDataArrayInfo(vtkDataArray *data, 
  int &datatype, std::string &dataname, int &numComponents)
{
#ifdef VTK_USE_MPI
  std::vector< vtkZPF_datainfo > datatypes(this->UpdateNumPieces);
  if (data) {
    ((vtkZPF_datainfo*)&datatypes[this->UpdatePiece])->datatype = data->GetDataType();
    ((vtkZPF_datainfo*)&datatypes[this->UpdatePiece])->numC     = data->GetNumberOfComponents();
    strncpy(((vtkZPF_datainfo*)&datatypes[this->UpdatePiece])->name, data->GetName(), 64);
  }
  vtkMPICommunicator* com = vtkMPICommunicator::SafeDownCast(
    this->Controller->GetCommunicator()); 
  int result = com->AllGather((char*)MPI_IN_PLACE, (char*)&datatypes[0], sizeof(vtkZPF_datainfo));
  for (int i=0; i<this->UpdateNumPieces; i++) {
    vtkZPF_datainfo &newdata = datatypes[i];
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
  datatypes[this->UpdatePiece] = points ? points->GetDataType() : -1;
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
vtkBoundingBox vtkZoltanV1PartitionFilter::GetGlobalBounds(vtkDataSet *input) 
{
  double bounds[6];
  input->GetBounds(bounds);
  vtkBoundingBox globalBounds(bounds);
  if (this->UpdateNumPieces>1) {
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
    globalBounds.SetMinPoint(bmin);
    globalBounds.SetMaxPoint(bmax);
  }
  return globalBounds;
}
//-------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::AllocateFieldArrays(vtkDataSetAttributes *fields) 
{
  int NumberOfFieldArrays = fields->GetNumberOfArrays();
  this->Controller->AllReduce(&NumberOfFieldArrays, &this->ZoltanCallbackData.NumberOfFields, 1, vtkCommunicator::MAX_OP);
  for (int i=0; i<this->ZoltanCallbackData.NumberOfFields; i++) {
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
void vtkZoltanV1PartitionFilter::InitializeZoltanLoadBalance()
{
  //***************************************************************
  //* Create a Zoltan library structure for this instance of load
  //* balancing.  Set the parameters and query functions that will
  //* govern the library's calculation.  See the Zoltan User's
  //* Guide for the definition of these and many other parameters.
  //***************************************************************

  this->ZoltanData = Zoltan_Create(this->GetMPIComm()); 

  // we don't need any debug info
  Zoltan_Set_Param(this->ZoltanData, "RCB_OUTPUT_LEVEL", "0");
  Zoltan_Set_Param(this->ZoltanData, "DEBUG_LEVEL", "0");

  // Method for subdivision
  Zoltan_Set_Param(this->ZoltanData, "LB_APPROACH", "REPARTITION");
  Zoltan_Set_Param(this->ZoltanData, "LB_METHOD",   "RCB");
  //  Zoltan_Set_Param(this->ZoltanData, "LB_METHOD", "PARMETIS");

  // Global and local Ids are a single integer
  Zoltan_Set_Param(this->ZoltanData, "NUM_GID_ENTRIES", "1"); 
  Zoltan_Set_Param(this->ZoltanData, "NUM_LID_ENTRIES", "0");

  // divide into N global and M local partitions
  std::stringstream global;
  global << this->UpdateNumPieces << ends;
  // just one local piece for now
  std::stringstream local;
  local << 1 << ends;
  //
  Zoltan_Set_Param(this->ZoltanData, "NUM_GLOBAL_PARTS", global.str().c_str());
  Zoltan_Set_Param(this->ZoltanData, "NUM_LOCAL_PARTS",  local.str().c_str());

  // All points have the same weight
  Zoltan_Set_Param(this->ZoltanData, "OBJ_WEIGHT_DIM", "0");
  Zoltan_Set_Param(this->ZoltanData, "RETURN_LISTS", "ALL");

  // RCB parameters
  // Zoltan_Set_Param(this->ZoltanData, "PARMETIS_METHOD", "PARTKWAY");
  Zoltan_Set_Param(this->ZoltanData, "RCB_RECOMPUTE_BOX", "0");
  Zoltan_Set_Param(this->ZoltanData, "REDUCE_DIMENSIONS", "0");

  // Don't allow very extended regions
  std::stringstream aspect;
  aspect << this->MaxAspectRatio << std::ends;  
  Zoltan_Set_Param(this->ZoltanData, "RCB_MAX_ASPECT_RATIO", aspect.str().c_str());

  // we need the cuts to get BBoxes for partitions later
  Zoltan_Set_Param(this->ZoltanData, "KEEP_CUTS", "1");

  // don't allow points on cut to be in different partitions
  // not likely/useful for particle data anyway
  Zoltan_Set_Param(this->ZoltanData, "RCB_RECTILINEAR_BLOCKS", "1"); 

  // Let Zoltan do the load balance step automatically
  // particles will be transferred as required between processes
  Zoltan_Set_Param(this->ZoltanData, "AUTO_MIGRATE", "1");
  Zoltan_Set_Param(this->ZoltanData, "MIGRATE_ONLY_PROC_CHANGES", "1"); 

  //
  // Query functions, to provide geometry to Zoltan 
  //
  Zoltan_Set_Num_Obj_Fn(this->ZoltanData,    get_number_of_objects_points, &this->ZoltanCallbackData);
//  Zoltan_Set_Obj_List_Fn(this->ZoltanData,   get_object_list_points,       &this->ZoltanCallbackData);
  Zoltan_Set_First_Obj_Fn(this->ZoltanData,  get_first_object_points,      &this->ZoltanCallbackData);
  Zoltan_Set_Next_Obj_Fn(this->ZoltanData,   get_next_object_points,       &this->ZoltanCallbackData);

  Zoltan_Set_Num_Geom_Fn(this->ZoltanData,   get_num_geometry,             &this->ZoltanCallbackData);
  if (this->ZoltanCallbackData.PointType==VTK_FLOAT) {
//    vtkDebugMacro(<<"Using float data pointers ");
    Zoltan_Set_Geom_Multi_Fn(this->ZoltanData, get_geometry_list<float>, &this->ZoltanCallbackData);
  }
  else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
//    vtkDebugMacro(<<"Using double data pointers ");
    Zoltan_Set_Geom_Multi_Fn(this->ZoltanData, get_geometry_list<double>, &this->ZoltanCallbackData);
  }

  //
  // Register functions for packing and unpacking data
  // by migration tools.  
  if (this->ZoltanCallbackData.PointType==VTK_FLOAT) {
    zsize_fn  f1 = zoltan_obj_size_function_points<float>;
    zpack_fn  f2 = zoltan_pack_obj_function_points<float>;
    zupack_fn f3 = zoltan_unpack_obj_function_points<float>;
    zprem_fn  f4 = zoltan_pre_migrate_function_points<float>; 
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &this->ZoltanCallbackData); 
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &this->ZoltanCallbackData); 
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &this->ZoltanCallbackData); 
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &this->ZoltanCallbackData); 
  }
  else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
    zsize_fn  f1 = zoltan_obj_size_function_points<double>;
    zpack_fn  f2 = zoltan_pack_obj_function_points<double>;
    zupack_fn f3 = zoltan_unpack_obj_function_points<double>;
    zprem_fn  f4 = zoltan_pre_migrate_function_points<double>;
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &this->ZoltanCallbackData);
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &this->ZoltanCallbackData);
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &this->ZoltanCallbackData);
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &this->ZoltanCallbackData);
  }
}
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::PartitionPoints(vtkInformation*,
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
//  int ghostLevel        = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
//  vtkDebugMacro(<<"Partition filter " << this->UpdatePiece << " Ghost level " << ghostLevel);
  //
  Timer = vtkSmartPointer<vtkTimerLog>::New();
  Timer->StartTimer();

  // Get input
  vtkIdType       numPoints = input->GetNumberOfPoints();
  vtkIdType        numCells = input->GetNumberOfCells();
  vtkDataArray    *inPoints = numPoints>0 ? input->GetPoints()->GetData() : NULL;
  vtkDebugMacro(<<"Partitioning { Points : " << numPoints << " } | { Cells " << numCells << " }");

  //--------------------------------------------------------------
  // Use Zoltan library to re-partition data in parallel
  // declare pointers which hold returned array info
  //--------------------------------------------------------------
  this->LoadBalanceData.numImport=0;
  this->LoadBalanceData.numExport=0;
  this->LoadBalanceData.importGlobalGids = NULL;
  this->LoadBalanceData.importLocalGids  = NULL; 
  this->LoadBalanceData.exportGlobalGids = NULL;
  this->LoadBalanceData.exportLocalGids  = NULL;
  this->LoadBalanceData.importProcs  = NULL;
  this->LoadBalanceData.importToPart = NULL;
  this->LoadBalanceData.exportProcs  = NULL;
  this->LoadBalanceData.exportToPart = NULL;

  // collect bounding boxes for all MPI ranks, we'll use it later to clamp the BSP limits
  vtkBoundingBox globalBounds = this->GetGlobalBounds(input);

  //
  // we make a temp copy of the input so we can add Ids if necessary
  //
  vtkSmartPointer<vtkPointSet> inputCopy;
  vtkSmartPointer<vtkPoints>   outPoints;
  if (this->UpdateNumPieces>1) {
    inputCopy.TakeReference(input->NewInstance());
    inputCopy->ShallowCopy(input);
    // create a new output points array to be filled
    outPoints = vtkSmartPointer<vtkPoints>::New();
    // if input had 0 points, make sure output is still setup correctly (float/double?)
    // collective exchanges will break if this is wrong as we may still receive data from another process
    // even though we are not sending any
    this->ZoltanCallbackData.PointType = this->GatherDataTypeInfo(input->GetPoints());
    outPoints->SetDataType(this->ZoltanCallbackData.PointType);
    output->SetPoints(outPoints);
  }
  else {
    // if only one process, we can just pass data through
    output->ShallowCopy(input);
    // vertex generation will fail if we don't set certain values
//    this->ZoltanCallbackData.OutPointCount = output->GetNumberOfPoints();
    this->ZoltanCallbackData.Output        = output;
    // make sure bounding boxes are set so that objects querying the output get the correct results
    this->BoxList.clear();
    this->BoxList.push_back(globalBounds);
    this->ExtentTranslator->SetNumberOfPieces(1);
    this->ExtentTranslator->SetBoundsForPiece(0, globalBounds);
    this->ExtentTranslator->InitWholeBounds();
    return 1;
  }

  //
  // Setup callbackdata structure as a user parameter to zoltan 
  // with partition info (H5Part reader can generate this)
  //
  this->ZoltanCallbackData.ProcessRank              = this->UpdatePiece;
  this->ZoltanCallbackData.Input                    = inputCopy;
  this->ZoltanCallbackData.Output                   = output;
  this->ZoltanCallbackData.InputPointsData          = inPoints ? inPoints->GetVoidPointer(0) : NULL;
//  this->ZoltanCallbackData.OutPointCount            = 0;
  this->ZoltanCallbackData.self                     = this;

  float ver;
  int zoltan_error = Zoltan_Initialize(0, NULL, &ver);
  if (zoltan_error != ZOLTAN_OK){
    printf("Zoltan initialization failed ...\n");
    return 0;
  }
  vtkDebugMacro(<<"Zoltan Initialized");

  //
  // if a process has zero points, we need to make dummy data arrays to allow 
  // space for when data gets sent in from other processes in the zoltan unpack function 
  // This also stops hangs during collective operations by ensuring all ranks) participate
  //
  vtkSmartPointer<vtkPointData> PointDataCopy = inputCopy->GetPointData();
  this->AllocateFieldArrays(PointDataCopy);
  vtkDebugMacro(<<"FieldArrayPointers (point) Initialized");

  //
  // When particle partitioning and exchanging halos, we need to track points
  // by Id. Mesh partitioning only needs Id offsets per process (no halos yet)
  this->SetupGlobalIds(inputCopy);

  // 
  // Set all the callbacks and user config parameters that will be used during the loadbalance
  //
  this->InitializeZoltanLoadBalance();

  //
  // Check the input to see if it has a bounds translator already initialized
  // with partition info. If this is set, it tells us that the data has already been loadbalanced
  // somewhere else and we can skip main partitioning step. (example: H5Part reader can generate this)
  //
  vtkExtentTranslator *translator = inInfo ? vtkExtentTranslator::SafeDownCast(
    inInfo->Get(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR())) : NULL;
  this->InputExtentTranslator = vtkBoundsExtentTranslator::SafeDownCast(translator);
  // if the extent translator has not been initialized well - don't use it
  if (this->InputExtentTranslator && this->InputExtentTranslator->GetNumberOfPieces()==0) {
    this->InputExtentTranslator = NULL;
  }
  this->ExtentTranslator->SetNumberOfPieces(this->UpdateNumPieces);

  //
  // if the input had a BoundsExtentTranslator, we can assume that all points and cells
  // are already well partitioned. We will use bounds info for halo exchange
  //
  if (this->InputExtentTranslator) {
    //
    // If we skipped the zoltan repartitioning, then do a copy (setup everything)
    // just like would have taken place at the start of the load balancing
    //
    if (this->ZoltanCallbackData.PointType==VTK_FLOAT) {
      zprem_fn  f4 = zoltan_pre_migrate_function_points<float>; 
      f4(&this->ZoltanCallbackData, 0, 0, 0, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL); 
    }
    else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
      zprem_fn  f4 = zoltan_pre_migrate_function_points<double>;
      f4(&this->ZoltanCallbackData, 0, 0, 0, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL); 
    }

    //
    // Get bounding boxes input ExtentTranslator
    //
    this->BoxList.clear();
    for (int p=0; p<this->UpdateNumPieces; p++) {
      vtkBoundingBox box;  
      box.SetBounds(this->InputExtentTranslator->GetBoundsForPiece(p));
      this->BoxList.push_back(box);
      this->ExtentTranslator->SetBoundsForPiece(p, this->InputExtentTranslator->GetBoundsForPiece(p));
    }
  }
  else {
    //
    // Zoltan can now partition our points. 
    // After this returns, we have redistributed points and the Output holds
    // the list of correct points/fields etc for each process
    //
    int zoltan_error = Zoltan_LB_Partition(
      this->ZoltanData,                            // input (all remaining fields are output)
      &this->LoadBalanceData.changes,              // 1 if partitioning was changed, 0 otherwise 
      &this->LoadBalanceData.numGidEntries,        // Number of integers used for a global ID
      &this->LoadBalanceData.numLidEntries,        // Number of integers used for a local ID
      &this->LoadBalanceData.numImport,            // Number of vertices to be sent to me
      &this->LoadBalanceData.importGlobalGids,     // Global IDs of vertices to be sent to me
      &this->LoadBalanceData.importLocalGids,      // Local IDs of vertices to be sent to me
      &this->LoadBalanceData.importProcs,          // Process rank for source of each incoming vertex
      &this->LoadBalanceData.importToPart,         // New partition for each incoming vertex
      &this->LoadBalanceData.numExport,            // Number of vertices I must send to other processes
      &this->LoadBalanceData.exportGlobalGids,     // Global IDs of the vertices I must send
      &this->LoadBalanceData.exportLocalGids,      // Local IDs of the vertices I must send
      &this->LoadBalanceData.exportProcs,          // Process to which I send each of the vertices
      &this->LoadBalanceData.exportToPart);        // Partition to which each vertex will belong

    if (zoltan_error != ZOLTAN_OK){
      printf("Zoltan_LB_Partition NOT OK...\n");
      MPI_Finalize();
      Zoltan_Destroy(&this->ZoltanData);
      exit(0);
    }

//#ifdef EXTRA_ZOLTAN_DEBUG
    vtkDebugMacro(<<"Partitioning complete on " << this->UpdatePiece << 
      " pack_count : " << pack_count <<
      " size_count : " << size_count <<
      " unpack_count : " << unpack_count 
     );
//#endif

    //
    // Get bounding boxes from zoltan and set them in the ExtentTranslator
    //
    this->BoxList.clear();
    for (int p=0; p<this->UpdateNumPieces; p++) {
      double bounds[6];
      int ndim;
      if (ZOLTAN_OK==Zoltan_RCB_Box(this->ZoltanData, p, &ndim, &bounds[0], &bounds[2], &bounds[4], &bounds[1], &bounds[3], &bounds[5])) {
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
    this->ExtentTranslator->InitWholeBounds();
  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::RequestData(vtkInformation* info,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // This is an abstract class and subclasses must override this function
  return 0;
}
//----------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::add_Id_to_interval_map(CallbackData *data, vtkIdType GID, vtkIdType LID) {
  vtkIdType diff = GID-LID;
  std::map<vtkIdType,vtkIdType>::reverse_iterator prev = data->ReceivedGlobalToLocalIdMap.rbegin();
  if (prev!=data->ReceivedGlobalToLocalIdMap.rend()) {
    vtkIdType key = prev->first;
    vtkIdType delta = prev->second;
    if (diff==delta) {
      // no need to add this value to the map as it has the same delta as the previous entry
    }
    else {
      // add a new entry to the 'interval' map
     data->ReceivedGlobalToLocalIdMap[GID] = diff;
    }
  }
  else {
    // add a new entry to the 'interval' map
    data->ReceivedGlobalToLocalIdMap[GID] = diff;
  }
}
//----------------------------------------------------------------------------
vtkIdType vtkZoltanV1PartitionFilter::global_to_local_Id(vtkIdType GID) {
  std::map<vtkIdType,vtkIdType>::iterator ub = this->ZoltanCallbackData.ReceivedGlobalToLocalIdMap.upper_bound(GID);
    ub--;
    vtkIdType key = ub->first;
    vtkIdType delta = ub->second;
    return GID - delta;
}
//----------------------------------------------------------------------------
MPI_Comm vtkZoltanV1PartitionFilter::GetMPIComm() {
#ifdef VTK_USE_MPI
  vtkMPICommunicator *communicator = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  MPI_Comm mpiComm = MPI_COMM_NULL;
  if (communicator) {
    mpiComm = *(communicator->GetMPIComm()->GetHandle());
  }
#else
  int mpiComm = 0;
#endif
  return mpiComm;
}
//----------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::InitializeFieldDataArrayPointers(
  CallbackData *callbackdata, 
  vtkFieldData *infielddata, 
  vtkFieldData *outfielddata,
  vtkIdType Nfinal) 
{
  // make sure pointers for copying field data are set to each data array start
  // and the size of each data tuple is recorded so we can quickly do a memcpy
  // for each ID that is transferred to/from this process
  callbackdata->TotalSizePerId = 0;
  callbackdata->MemoryPerTuple.clear();
  callbackdata->InputArrayPointers.clear();
  callbackdata->OutputArrayPointers.clear();
  callbackdata->NumberOfFields = infielddata->GetNumberOfArrays();
  for (int i=0; i<callbackdata->NumberOfFields; i++) {
    vtkDataArray *iarray = infielddata->GetArray(i);
    vtkDataArray *oarray = outfielddata->GetArray(i);
    oarray->Resize(Nfinal);
    oarray->SetNumberOfTuples(Nfinal);
    callbackdata->InputArrayPointers.push_back(iarray->GetVoidPointer(0));
    callbackdata->OutputArrayPointers.push_back(oarray->GetVoidPointer(0));
    // we need to know the amount of data to copy for each array tuple
    int Nc = iarray->GetNumberOfComponents();
    int Ns = iarray->GetDataTypeSize();
    callbackdata->MemoryPerTuple.push_back(Nc*Ns);
    callbackdata->TotalSizePerId += Nc*Ns;
  }
}
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::ManualPointMigrate(PartitionInfo &partitioninfo, bool useoutput)
{
  //
  // Ask zoltan to create the inverse map of whos sends/receives from who
  //
  int           num_known        = partitioninfo.GlobalIds.size(); 
  int           num_found        = 0;
  ZOLTAN_ID_PTR found_global_ids = NULL;
  ZOLTAN_ID_PTR found_local_ids  = NULL;
  int          *found_procs      = NULL;
  int          *found_to_part    = NULL;
  //
  int zoltan_error = Zoltan_Invert_Lists(this->ZoltanData, 
    num_known,
    num_known>0 ? &partitioninfo.GlobalIds[0] : NULL,
    /*num_known>0 ? &partitioninfo.LocalIds[0]  : */NULL,
    num_known>0 ? &partitioninfo.Procs[0]     : NULL,
    /*num_known>0 ? &partitioninfo.Procs[0]     : */NULL,
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

  this->ZoltanCallbackData.CopyFromOutput = useoutput;
  //
  // Before sending, set the pre-migrate function to the add function
  // which assumes we are receiving new points, but not removing any existing ones
  //
  if (this->ZoltanCallbackData.PointType==VTK_FLOAT) {
    zprem_fn f4 = zoltan_pre_migrate_function_points_add<float>;
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &this->ZoltanCallbackData);
  }
  else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
    zprem_fn f4 = zoltan_pre_migrate_function_points_add<double>;
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &this->ZoltanCallbackData);
  }

  //
  // Now let zoltan perform the send/receive exchange of particles
  //
  zoltan_error = Zoltan_Migrate (this->ZoltanData,
    num_found,
    found_global_ids,
    found_local_ids,
    found_procs,
    found_to_part,
    num_known,
    num_known>0 ? &partitioninfo.GlobalIds[0] : NULL,
    /*num_known>0 ? &partitioninfo.LocalIds[0]  : */NULL,
    num_known>0 ? &partitioninfo.Procs[0]     : NULL,
    /*num_known>0 ? &partitioninfo.Procs[0]     : */NULL
    );

  //
  // Release the arrays allocated during Zoltan_Invert_Lists
  //
  Zoltan_LB_Free_Part(
    &found_global_ids, 
    &found_local_ids, 
    &found_procs, 
    &found_to_part);

  return num_found;
}


  // Description:
  //   Initialize the cuts with arrays of information.  This type of
  //   information would be obtained from a graph partitioning software
  //   package like Zoltan.
  //
  //   bounds - the bounds (xmin, xmax, ymin, ymax, zmin, zmax) of the
  //             space being partitioned
  //   ncuts - the number cuts, also the size of the following arrays
  //   dim   - the dimension along which the cut is made (x/y/z - 0/1/2)
  //   coord - the location of the cut along the axis
  //   lower - array index for the lower region bounded by the cut
  //   upper - array index for the upper region bounded by the cut
  //   lowerDataCoord - optional upper bound of the data in the lower region
  //   upperDataCoord - optional lower bound of the data in the upper region
  //   npoints - optional number of points in the spatial region

//----------------------------------------------------------------------------
vtkSmartPointer<vtkPKdTree> vtkZoltanV1PartitionFilter::CreatePkdTree()
{
  vtkSmartPointer<vtkBSPCuts> cuts = vtkSmartPointer<vtkBSPCuts>::New();
  //
  RCB_STRUCT *rcb = (RCB_STRUCT *) (this->ZoltanData->LB.Data_Structure);
  struct rcb_tree *treept = rcb->Tree_Ptr;
  //
  if (treept[0].dim < 0) {
    vtkErrorMacro(<<"RCB Tree invalid");
    return NULL;
  }

  // Zoltan might have remapped the partitions from the default depth first ordering
  std::vector<int> remapping(this->UpdateNumPieces, -1);
  for (int i=0; i<this->ZoltanData->LB.Num_Global_Parts; i++) {
    if (this->ZoltanData->LB.Remap) {
      remapping[this->ZoltanData->LB.Remap[i]] = i;
    }
    else {
      remapping[i] = i;
    }
  }

  // list we will pass to CreateCuts
  std::vector<int>    cut_axis;
  std::vector<double> cut_position;
  std::vector<int>    cut_lower;
  std::vector<int>    cut_upper;

  // use a stack to traverse the RCB cut tree
  typedef std::pair<struct rcb_tree*, int> cutpair;
  std::stack<cutpair> tree_stack;
  struct rcb_tree *nodept = &treept[0];
  
  // push the root onto the stack
  if (nodept->right_leaf) {
    tree_stack.push(cutpair(&treept[nodept->right_leaf],-1));
  }

  // Traverse the RCB cut tree creating nodes for the the vtk cuts
  int RegionId = 0;
  while (!tree_stack.empty()) {
    nodept = tree_stack.top().first;
    int node_parent = tree_stack.top().second;
    tree_stack.pop();

    // create a node for this cut, append it to the list
    cut_position.push_back(nodept->cut);
    cut_axis.push_back(nodept->dim);

    // we don't know where child nodes will be in the list yet
    cut_lower.push_back(-1); 
    cut_upper.push_back(-1); 

    // the new node we've just added to the list has this index
    int new_node_index = cut_position.size()-1;

    // if the parent of this node had a missing child pointer/index
    // and this is the child required, set the child index
    if (node_parent>=0) {
      // the left node will always be done first, so check it
      if (cut_lower[node_parent] == -1) {
        cut_lower[node_parent] = new_node_index;
      }
      // if it has already been filled in, we must be the right node
      else if (cut_upper[node_parent] == -1) {
        cut_upper[node_parent] = new_node_index;
      }
    }

    // if this is a new 'normal' parent node with 2 children
    if (nodept->right_leaf > 0) {
      // push two new nodes onto the stack : always push right(upper) then left(lower) 
      // to make sure the left is processed/popped first and the array offset/indexes are correct
      // (since we don't know the child locations when we create them, save the parent ID
      // so that we can fill in the missing info when the children are created)
      tree_stack.push(cutpair(&(treept[nodept->right_leaf]),new_node_index));
      tree_stack.push(cutpair(&(treept[nodept->left_leaf]), new_node_index));
    }
    // if the node is a leaf, hanging from an unbalanced node (e.g. 7 leaf tree)
    else if (nodept->left_leaf > 0) {
      // push the node which has more cuts
      tree_stack.push(cutpair(&(treept[nodept->left_leaf]),new_node_index));
      // create one leaf node, set the (region) Id for BSPCuts to use
      cut_position.push_back(0.0);
      cut_axis.push_back(0);
      cut_lower.push_back(-remapping[RegionId]); 
      cut_upper.push_back(-1); 
      // point the parent node to this leaf
      cut_upper[new_node_index] = new_node_index+1; 
      RegionId++;
    }
    else {
      // create two leaf nodes, set the (region) Id for BSPCuts to use
      cut_position.push_back(0.0);
      cut_axis.push_back(0);
      cut_lower.push_back(-remapping[RegionId]); 
      cut_upper.push_back(-1); 
      RegionId++;
      cut_position.push_back(0.0);
      cut_axis.push_back(0);
      cut_lower.push_back(-remapping[RegionId]); 
      cut_upper.push_back(-1); 
      RegionId++;
      // the children are created now, so we know what the indices are
      cut_lower[new_node_index] = new_node_index+1; 
      cut_upper[new_node_index] = new_node_index+2; 
    }
  };
 
  cuts->CreateCuts(
    this->ExtentTranslator->GetWholeBounds(),
    cut_axis.size(),
    &cut_axis[0],
    &cut_position[0],
    &cut_lower[0],
    &cut_upper[0],
    NULL,NULL,NULL);

  this->KdTree = vtkSmartPointer<vtkPKdTree2>::New();
  this->KdTree->SetController(this->Controller);
  this->KdTree->SetCuts(cuts);
  vtkPKdTree2::SafeDownCast(this->KdTree)->BuildLocator(
    this->ExtentTranslator->GetWholeBounds(), this->ZoltanData->LB.Remap /* &remapping[0]*/, this->UpdateNumPieces);

  // in order for the rest of paraview to work properly, we need to 'build' the locator internal
  // using some fake cells to ge things going. Each process will use its bounding box centre
  // as a single cell and partition it.
/*
  vtkNew<vtkUnstructuredGrid> grid; 
  vtkBoundingBox *box = this->GetPartitionBoundingBox(this->UpdatePiece);
  box->Inflate( -box->GetDiagonalLength()/1000.0 );

  //
  vtkNew<vtkCellArray> cells;
#define HEXAHEDRON_CELLS
#ifdef HEXAHEDRON_CELLS
  const double *minpt = box->GetMinPoint();
  const double *maxpt = box->GetMaxPoint();
  //
  double P0[3] = {minpt[0], minpt[1], minpt[2]};
  double P1[3] = {maxpt[0], minpt[1], minpt[2]};
  double P2[3] = {maxpt[0], maxpt[1], minpt[2]};
  double P3[3] = {minpt[0], maxpt[1], minpt[2]};
  double P4[3] = {minpt[0], minpt[1], maxpt[2]};
  double P5[3] = {maxpt[0], minpt[1], maxpt[2]};
  double P6[3] = {maxpt[0], maxpt[1], maxpt[2]};
  double P7[3] = {minpt[0], maxpt[1], maxpt[2]};

  vtkNew<vtkPoints> points;
  points->InsertNextPoint(P0);
  points->InsertNextPoint(P1);
  points->InsertNextPoint(P2);
  points->InsertNextPoint(P3);
  points->InsertNextPoint(P4);
  points->InsertNextPoint(P5);
  points->InsertNextPoint(P6);
  points->InsertNextPoint(P7);

  vtkIdType pts[8] = { 0,1,2,3,4,5,6,7};
  cells->InsertNextCell(8, pts);
  grid->SetCells(VTK_HEXAHEDRON, cells.GetPointer());
#else
  double P0[3];
  box->GetCenter(P0);
  vtkNew<vtkPoints> points;
  points->InsertNextPoint(P0);
  vtkIdType pts[1] = { 0 };
  cells->InsertNextCell(1, pts);
  grid->SetCells(VTK_VERTEX, cells.GetPointer());
#endif

  grid->SetPoints(points.GetPointer());
*/
//  this->KdTree->AssignRegionsRoundRobin();
//  this->KdTree->SetDataSet(grid.GetPointer());
//  this->KdTree->BuildLocator();


//  const int *regionmap = this->KdTree->GetRegionAssignmentMap();
//


  return KdTree;
}
