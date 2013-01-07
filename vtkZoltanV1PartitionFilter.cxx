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
vtkStandardNewMacro(vtkZoltanV1PartitionFilter);
vtkCxxSetObjectMacro(vtkZoltanV1PartitionFilter, Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
// Zoltan callback which returns number of objects participating in exchange
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::get_number_of_objects_points(void *data, int *ierr)
{
  int res = static_cast<CallbackData*>(data)->InputNumberOfLocalPoints; 
  *ierr = (res < 0) ? ZOLTAN_FATAL : ZOLTAN_OK;
  return res;
}
//----------------------------------------------------------------------------
// Zoltan callback which fills the Ids for each object in the exchange
//----------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::get_object_list_points(void *data, int sizeGID, int sizeLID,
  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr)
{
  CallbackData *mesh = (CallbackData*)data;
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
  CallbackData *mesh = (CallbackData*)data;
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
int vtkZoltanV1PartitionFilter::zoltan_obj_size_func_points(void *data, 
  int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, 
  ZOLTAN_ID_PTR local_id, int *ierr)
{
  INC_SIZE_COUNT
  CallbackData *mesh = (CallbackData*)data;
  *ierr = ZOLTAN_OK;
  return mesh->TotalSizePerId + sizeof(T)*3;
}
//----------------------------------------------------------------------------
// Zoltan callback to pack all the data for one point into a buffer
//----------------------------------------------------------------------------
template<typename T>
void vtkZoltanV1PartitionFilter::zoltan_pack_obj_func_points(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size, char *buf, int *ierr)
{
  INC_PACK_COUNT
  CallbackData *mesh = (CallbackData*)data;
  vtkIdType GID = *global_id;
  vtkIdType LID = *local_id;
  //
  for (int i=0; i<mesh->NumberOfFields; i++) {
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
void vtkZoltanV1PartitionFilter::zoltan_unpack_obj_func_points(void *data, int num_gid_entries,
  ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr)
{
  INC_UNPACK_COUNT
  if (num_gid_entries != 1) {
    *ierr = ZOLTAN_FATAL;
    return;
  }
  CallbackData *mesh = (CallbackData*)data;
  vtkIdType GID = *global_id;
  //
  vtkPointData *inPD  = mesh->Input->GetPointData();
  vtkPointData *outPD = mesh->Output->GetPointData();
  //
  for (int i=0; i<mesh->NumberOfFields; i++) {
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
void vtkZoltanV1PartitionFilter::zoltan_pre_migrate_func_points(void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  CallbackData *mesh = (CallbackData*)data;
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
  for (int i=0; i<mesh->NumberOfFields; i++) {
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
// vtkZoltanV1PartitionFilter :: implementation 
//----------------------------------------------------------------------------
vtkZoltanV1PartitionFilter::vtkZoltanV1PartitionFilter()
{
  this->UpdatePiece               = 0;
  this->UpdateNumPieces           = 1;
  this->IdChannelArray            = NULL;
  this->MaxAspectRatio            = 5.0;
  this->ExtentTranslator          = vtkBoundsExtentTranslator::New();
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
vtkSmartPointer<vtkIdTypeArray> vtkZoltanV1PartitionFilter::GenerateGlobalIds(vtkIdType Npoints, vtkIdType Ncells, const char *ptidname)
{
  vtkSmartPointer<vtkIdTypeArray> ptIds = vtkSmartPointer<vtkIdTypeArray>::New();

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

  //
  // Global point IDs generated here
  //
  vtkIdType offset = this->ZoltanCallbackData.ProcessOffsetsPointId[this->UpdatePiece];
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
//    databox.Inflate(this->GhostCellOverlap);
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
void vtkZoltanV1PartitionFilter::SetupFieldArrayPointers(vtkDataSetAttributes *fields) 
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
  Zoltan_Set_Param(this->ZoltanData, "NUM_LID_ENTRIES", "1");

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

  //
  // Query functions, to provide geometry to Zoltan 
  //
  Zoltan_Set_Num_Obj_Fn(this->ZoltanData,    get_number_of_objects_points, &this->ZoltanCallbackData);
  Zoltan_Set_Obj_List_Fn(this->ZoltanData,   get_object_list_points,       &this->ZoltanCallbackData);
  Zoltan_Set_Num_Geom_Fn(this->ZoltanData,   get_num_geometry,             &this->ZoltanCallbackData);
  if (this->ZoltanCallbackData.PointType==VTK_FLOAT) {
    vtkDebugMacro(<<"Using float data pointers ");
    Zoltan_Set_Geom_Multi_Fn(this->ZoltanData, get_geometry_list<float>, &this->ZoltanCallbackData);
  }
  else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
    vtkDebugMacro(<<"Using double data pointers ");
    Zoltan_Set_Geom_Multi_Fn(this->ZoltanData, get_geometry_list<double>, &this->ZoltanCallbackData);
  }

  //
  // Register functions for packing and unpacking data
  // by migration tools.  
  if (this->ZoltanCallbackData.PointType==VTK_FLOAT) {
    zsize_fn  f1 = zoltan_obj_size_func_points<float>;
    zpack_fn  f2 = zoltan_pack_obj_func_points<float>;
    zupack_fn f3 = zoltan_unpack_obj_func_points<float>;
    zprem_fn  f4 = zoltan_pre_migrate_func_points<float>; 
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &this->ZoltanCallbackData); 
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &this->ZoltanCallbackData); 
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &this->ZoltanCallbackData); 
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &this->ZoltanCallbackData); 
  }
  else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
    zsize_fn  f1 = zoltan_obj_size_func_points<double>;
    zpack_fn  f2 = zoltan_pack_obj_func_points<double>;
    zupack_fn f3 = zoltan_unpack_obj_func_points<double>;
    zprem_fn  f4 = zoltan_pre_migrate_func_points<double>;
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
  int ghostLevel        = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
  vtkDebugMacro(<<"Partition filter " << this->UpdatePiece << " Ghost level " << ghostLevel);
  //
  Timer = vtkSmartPointer<vtkTimerLog>::New();
  Timer->StartTimer();

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

  // if input had 0 points, make sure output is still setup correctly (float/double?)
  // collective exchanges will break if this is wrong as we may still receive data
  this->ZoltanCallbackData.PointType = this->GatherDataTypeInfo(input->GetPoints());
  outPoints->SetDataType(this->ZoltanCallbackData.PointType);
  output->SetPoints(outPoints);

  vtkDebugMacro(<<"Initializing Zoltan on " << this->UpdatePiece);
  float ver;
  int zoltan_error = Zoltan_Initialize(0, NULL, &ver);
  if (zoltan_error != ZOLTAN_OK){
    printf("Zoltan initialization failed ...\n");
    return 0;
  }
  vtkDebugMacro(<<"Zoltan Initialized on " << this->UpdatePiece);

  //
  // Setup mesh structure as a user parameter to zoltan 
  // with partition info (H5Part reader can generate this)
  //
  this->ZoltanCallbackData.ProcessRank              = this->UpdatePiece;
  this->ZoltanCallbackData.Input                    = inputCopy;
  this->ZoltanCallbackData.Output                   = output;
  this->ZoltanCallbackData.InputNumberOfLocalPoints = numPoints;
  this->ZoltanCallbackData.InputNumberOfLocalCells  = numCells;
  this->ZoltanCallbackData.InputPointsData          = inPoints ? inPoints->GetVoidPointer(0) : NULL;
  this->ZoltanCallbackData.OutputPoints             = outPoints;
  this->ZoltanCallbackData.OutPointCount            = 0;

  //
  // if a process has zero points, we need to make dummy point data arrays to allow 
  // space for when data gets sent in from other processes in the zoltan unpack function 
  // This also stops hangs during collective operations by ensuring all ranks participate
  //
  vtkSmartPointer<vtkPointData> PointDataCopy = inputCopy->GetPointData();
  this->SetupFieldArrayPointers(PointDataCopy);

  //
  // Global Ids : always do them after other point arrays 
  //
  if (this->IdChannelArray) {
    this->IdsName = this->IdChannelArray;
  }
  if (this->IdsName.empty() || this->IdsName==std::string("Not available")) {
    this->IdsName = "PPF_PointIds";
  } 

  vtkSmartPointer<vtkDataArray> Ids = NULL;
  Ids = PointDataCopy->GetArray(this->IdsName.c_str());
  if (!Ids) {
    // Try loading the user supplied global ids.
    Ids = PointDataCopy->GetGlobalIds();
  }
  if (!Ids) {
    // Generate our own since none exist
    Ids = this->GenerateGlobalIds(numPoints, numCells, this->IdsName.c_str());
    inputCopy->GetPointData()->AddArray(Ids);
    // and increment the mesh field count
    this->ZoltanCallbackData.NumberOfFields++;
  }

  //
  this->InitializeZoltanLoadBalance();

  //
  // Check the input to see if it has a bounds translator already initialized
  // with partition info (example: H5Part reader can generate this)
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
      zprem_fn  f4 = zoltan_pre_migrate_func_points<float>; 
      f4(&this->ZoltanCallbackData, 0, 0, 0, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL); 
    }
    else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
      zprem_fn  f4 = zoltan_pre_migrate_func_points<double>;
      f4(&this->ZoltanCallbackData, 0, 0, 0, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL); 
    }

    //
    // Get bounding boxes input ExtentTranslator
    //
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

#ifdef EXTRA_ZOLTAN_DEBUG
    vtkDebugMacro(<<"Partitioning complete on " << this->UpdatePiece << 
      " pack_count : " << pack_count <<
      " size_count : " << size_count <<
      " unpack_count : " << unpack_count 
     );
#endif

    //
    // Get bounding boxes from zoltan and set them in the ExtentTranslator
    //
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
  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkZoltanV1PartitionFilter::RequestData(vtkInformation* info,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
/*
    //
    //
    //
    vtkDebugMacro(<<"Process " << this->UpdatePiece << " Points Output : " << this->ZoltanCallbackData.OutPointCount);
    //  for (int i=0; i<this->ZoltanCallbackData.NumberOfFields; i++) {
    //    vtkDataArray *darray = output->GetPointData()->GetArray(i);
    //    vtkDebugMacro(<<"Process " << this->UpdatePiece << " Array Output : " << darray->GetNumberOfTuples());
    //  }

    //
    // If polydata create Vertices for each point
    //
    if (vtkPolyData::SafeDownCast(output)) {
      vtkIdType *arraydata = cells->WritePointer(this->ZoltanCallbackData.OutPointCount, 2*this->ZoltanCallbackData.OutPointCount);
      for (int i=0; i<this->ZoltanCallbackData.OutPointCount; i++) {
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
    Zoltan_Destroy(&this->ZoltanData);

  }

  this->Controller->Barrier();
  Timer->StopTimer();
  vtkDebugMacro(<<"Particle partitioning : " << Timer->GetElapsedTime() << " seconds");
*/
  return 0;
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

