/*=========================================================================

  Module : vtkZoltanBasePartitionFilter

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
#include "vtkInformationDoubleKey.h"
#include "vtkInformationDoubleVectorKey.h"
#include "vtkInformationIntegerKey.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
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
#include "vtkZoltanBasePartitionFilter.h"
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>
#include <algorithm>
#include <map>
#include <stack>
#include <iostream>
#include <ostream>
#include <sstream>
#include <iterator>
//
#include "zz_const.h"
#include "rcb.h"
//
//----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkZoltanBasePartitionFilter, Controller, vtkMultiProcessController);
vtkInformationKeyMacro(vtkZoltanBasePartitionFilter, ZOLTAN_SAMPLE_RESOLUTION, DoubleVector);
vtkInformationKeyMacro(vtkZoltanBasePartitionFilter, ZOLTAN_SAMPLE_ORIGIN,     DoubleVector);
//----------------------------------------------------------------------------

#ifdef ZOLTAN_DEBUG_OUTPUT
int vtkZoltanBasePartitionFilter::pack_count = 0;
int vtkZoltanBasePartitionFilter::unpack_count = 0;
int vtkZoltanBasePartitionFilter::size_count = 0;
#endif

//----------------------------------------------------------------------------
#if defined ZOLTAN_DEBUG_OUTPUT && !defined VTK_WRAPPING_CXX

# undef vtkDebugMacro
# define vtkDebugMacro(msg)  \
   DebugSynchronized(this->UpdatePiece, this->UpdateNumPieces, this->Controller, msg);

# undef  vtkErrorMacro
# define vtkErrorMacro(a) vtkDebugMacro(a)
#endif
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#include "vtkZoltanBasePartitionFilter.txx"
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
// Zoltan callback which returns number of objects participating in exchange
//----------------------------------------------------------------------------
int vtkZoltanBasePartitionFilter::get_number_of_objects_points(void *data, int *ierr)
{
  int res = static_cast<CallbackData*>(data)->Input->GetNumberOfPoints();
  *ierr = (res < 0) ? ZOLTAN_FATAL : ZOLTAN_OK;
  return res;
}

//----------------------------------------------------------------------------
// Zoltan callback which fills the Ids/weights for each object in the exchange
//----------------------------------------------------------------------------
void vtkZoltanBasePartitionFilter::get_object_list_points(void *data, int sizeGID, int sizeLID,
    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
    int wgt_dim, float *obj_wgts, int *ierr)
{
    CallbackData *callbackdata = static_cast<CallbackData*>(data);

    *ierr = ZOLTAN_OK;
    vtkIdType N = callbackdata->Input->GetNumberOfPoints();
    for (vtkIdType i=0; i<N; ++i) {
        globalID[i] = i + callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];
        if (wgt_dim && callbackdata->self->weights_data_ptr) {
            obj_wgts[i] = ((float*)callbackdata->self->weights_data_ptr)[i];
        }
    }
}

//----------------------------------------------------------------------------
// Zoltan callback which returns the dimension of geometry (3D for us)
//----------------------------------------------------------------------------
int vtkZoltanBasePartitionFilter::get_num_geometry(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 3;
}

//----------------------------------------------------------------------------
// Zoltan callback which does nothing, we register this during load balance
// when we do not want any pre migration operations (we manually handle it)
//----------------------------------------------------------------------------
void vtkZoltanBasePartitionFilter::zoltan_pre_migrate_function_null(
  void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  *ierr = ZOLTAN_OK;
}

//----------------------------------------------------------------------------
// vtkZoltanBasePartitionFilter :: implementation
//----------------------------------------------------------------------------
vtkZoltanBasePartitionFilter::vtkZoltanBasePartitionFilter()
{
  this->UpdatePiece                    = 0;
  this->UpdateNumPieces                = 1;
  this->MaxAspectRatio                 = 5.0;
  this->GhostHaloSize                  = 0.0;
  this->ExtentTranslator               = vtkSmartPointer<vtkBoundsExtentTranslator>::New();
  this->InputExtentTranslator          = NULL;
  this->ZoltanData                     = NULL;
  this->InputDisposable                = 0;
  this->KeepInversePointLists          = 0;
  this->PointWeightsArrayName          = NULL;
  this->weights_data_ptr               = NULL;
  this->ImbalanceValue                 =-1.0; // invalid
  this->Controller                     = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
}
//----------------------------------------------------------------------------
vtkZoltanBasePartitionFilter::~vtkZoltanBasePartitionFilter()
{
  // make sure we do not leave any memory unfreed
  if (this->MigrateLists.num_found!=-1) {
    //vtkDebugMacro("Deleting InverseMigrationLists");
    Zoltan_LB_Free_Part(
      &this->MigrateLists.found_global_ids,
      &this->MigrateLists.found_local_ids,
      &this->MigrateLists.found_procs,
      &this->MigrateLists.found_to_part);
    // set to zero so we know data has been deleted
    this->MigrateLists.num_found = -1;
  }
  //
  this->SetPointWeightsArrayName(NULL);
  //
  this->SetController(NULL);
}

//----------------------------------------------------------------------------
void vtkZoltanBasePartitionFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkZoltanBasePartitionFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter supports any vtkPointSet type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkZoltanBasePartitionFilter::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // output should be the same as the input (a subclass of vtkPointSet)
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
vtkBoundingBox *vtkZoltanBasePartitionFilter::GetPartitionBoundingBox(int partition)
{
  if (partition<this->BoxList.size()) {
    return &this->BoxList[partition];
  }
  vtkErrorMacro("Partition not found in Bounding Box list");
  return NULL;
}

//----------------------------------------------------------------------------
vtkBoundingBox *vtkZoltanBasePartitionFilter::GetPartitionBoundingBoxHalo(int partition)
{
  if (partition<this->BoxListWithHalo.size()) {
    return &this->BoxListWithHalo[partition];
  }
  vtkErrorMacro("Partition not found in Bounding Box list");
  return NULL;
}

//----------------------------------------------------------------------------
void vtkZoltanBasePartitionFilter::ComputeIdOffsets(vtkIdType Npoints, vtkIdType Ncells)
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
}

//----------------------------------------------------------------------------
struct vtkZPF_datainfo {
  int  datatype;
  int  numC;
  int  attributes;
  char name[256];
  vtkZPF_datainfo() : datatype(-1), numC(-1), attributes(-1) {};
};

struct vtkZPF_datasetinfo {
  int  point_data_type;
  int  polydata_types;
  vtkZPF_datasetinfo() : point_data_type(-1), polydata_types(-1) {};
};

//----------------------------------------------------------------------------
bool vtkZoltanBasePartitionFilter::GatherDataArrayInfo(
    vtkDataArray *data, vtkDataSetAttributes *attribs,
  int &datatype, std::string &dataname, int &numComponents, int &attrib)
{
#ifdef VTK_USE_MPI
  std::vector< vtkZPF_datainfo > datatypes(this->UpdateNumPieces);
  if (data) {
    vtkZPF_datainfo &info = datatypes[this->UpdatePiece];
    info.datatype   = data->GetDataType();
    info.numC       = data->GetNumberOfComponents();
    info.attributes = -1;
    strncpy(((vtkZPF_datainfo*)&datatypes[this->UpdatePiece])->name, data->GetName(), 64);
    // we need the index of the array to get attribute type if present
    int index = -1;
    vtkAbstractArray *dummy = attribs->GetAbstractArray(data->GetName(), index);
    if (dummy==data) {
        info.attributes = attribs->IsArrayAnAttribute(index);
        if (info.attributes>=0) {
            debug_no_sync(<< data->GetName() << " is attribute " << info.attributes);
        }
    }

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
      attrib = newdata.attributes;
    }
  }
  return (result == 1) ;
#else
  return 1;
#endif
}

//----------------------------------------------------------------------------
int vtkZoltanBasePartitionFilter::GatherDataTypeInfo(vtkDataSet *input, vtkPoints *points)
{
#ifdef VTK_USE_MPI
  if (this->UpdateNumPieces==1) {
      return points->GetDataType();
  }
  std::vector< vtkZPF_datasetinfo > datatypes(this->UpdateNumPieces, vtkZPF_datasetinfo());
  datatypes[this->UpdatePiece].point_data_type = points ?
      (points->GetNumberOfPoints() ? points->GetDataType() : -1) : -1;
  int datatype = datatypes[this->UpdatePiece].point_data_type;
  //
  vtkPolyData *pdata = vtkPolyData::SafeDownCast(input);
  if (pdata) {
      int v = pdata->GetNumberOfVerts()>0  ? 1 : 0;
      int l = pdata->GetNumberOfLines()>0  ? 2 : 0;
      int p = pdata->GetNumberOfPolys()>0  ? 4 : 0;
      int s = pdata->GetNumberOfStrips()>0 ? 8 : 0;
      datatypes[this->UpdatePiece].polydata_types = v + l + p + s;
  }
  vtkMPICommunicator* com = vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  int result = com->AllGather((const char *)MPI_IN_PLACE, (char *)&datatypes[0], sizeof(vtkZPF_datasetinfo));
  this->polydata_types = 0;
  for (int i=0; i<this->UpdateNumPieces; i++) {
    int &newdatatype = datatypes[i].point_data_type;
    if (datatype==-1 && newdatatype!=-1) {
      datatype = newdatatype;
    }
    else if (datatype!=-1 && newdatatype!=-1 && newdatatype!=datatype) {
      vtkErrorMacro("Fatal datatype error in Point DataType Gather");
    }
    this->polydata_types |= datatypes[i].polydata_types;
  }
  return datatype;
#else
  return pInput->GetPoints()->GetDataType();
#endif
}

//-------------------------------------------------------------------------
vtkBoundingBox vtkZoltanBasePartitionFilter::GetGlobalBounds(vtkDataSet *input)
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
void vtkZoltanBasePartitionFilter::AllocateFieldArrays(vtkDataSetAttributes *fields, vtkDataSetAttributes *fieldcopy)
{
  int NumberOfFieldArrays = fields->GetNumberOfArrays();
  this->Controller->AllReduce(&NumberOfFieldArrays, &this->ZoltanCallbackData.NumberOfFields, 1, vtkCommunicator::MAX_OP);
  vtkDebugMacro("NULL data found, NumberOfFieldArrays " << this->ZoltanCallbackData.NumberOfFields);
  for (int i=0; i<this->ZoltanCallbackData.NumberOfFields; i++) {
    vtkSmartPointer<vtkDataArray> darray = fields->GetArray(i);
    //
    int correctType = -1, numComponents = -1, attrib = -1;
    std::string correctName;
    this->GatherDataArrayInfo(darray, fields, correctType, correctName, numComponents, attrib);
    if (!darray) {
      vtkDebugMacro("NULL data found, used MPI_Gather to find :"
        << " DataType " << correctType
        << " Name " << correctName.c_str()
        << " NumComponents " << numComponents);
      darray.TakeReference(vtkDataArray::CreateDataArray(correctType));
      darray->SetNumberOfComponents(numComponents);
      darray->SetName(correctName.c_str());
      fieldcopy->AddArray(darray);
      if (attrib!=-1) {
          vtkDebugMacro("setting " << correctName.c_str() << " attrib " << attrib);
          fieldcopy->SetActiveAttribute(correctName.c_str(), attrib);
      }
      else {
          vtkDebugMacro(<< correctName.c_str() << " no attrib ");
      }
    }
    else {
        // twice because they are synchronized with the debug messages above
        vtkDebugMacro("Data Ok " << darray->GetName());
        vtkDebugMacro("Data Ok " << darray->GetName());
    }
  }
  vtkDebugMacro("AllocateFieldArrays completed");
}

//-------------------------------------------------------------------------
void vtkZoltanBasePartitionFilter::SetupPointWeights(vtkDataSetAttributes *fields) {
    // get the array that is used for weights, check it the same as coords because
    // @TODO, zoltan only allows one template param for coords and weights
    // so we can't use different types yet
    vtkDataArray *weightsArray = this->PointWeightsArrayName ?
        fields->GetArray(this->PointWeightsArrayName) : NULL;
    this->weights_data_ptr = NULL;
    //
    if (weightsArray && (VTK_FLOAT != weightsArray->GetDataType())) {
        vtkWarningMacro(<<"Weights datatype must be the same as coordinate type");
        weightsArray = NULL;
    }
    if (weightsArray) {
        // get the pointer to the actual data, if the array is empty, set the pointer
        // to some non NULL value so that later checks for weights present don't fail
        this->weights_data_ptr = weightsArray->GetVoidPointer(0) ? weightsArray->GetVoidPointer(0) : (void*)0xFFFFFFFF;
    }
}
//-------------------------------------------------------------------------
int vtkZoltanBasePartitionFilter::RequestUpdateExtent(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // We require preceding filters to refrain from creating ghost cells.
  int piece, numPieces, ghostLevels = 0;

  piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(), piece);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(), numPieces);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(), ghostLevels);
  inInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);

  return 1;
}

//----------------------------------------------------------------------------
int vtkZoltanBasePartitionFilter::RequestInformation(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
#ifdef VTK_USE_MPI
  vtkMPICommunicator *communicator =
    vtkMPICommunicator::SafeDownCast(this->Controller->GetCommunicator());
  int maxpieces = communicator ? communicator->GetNumberOfProcesses() : 1;
  //
  this->UpdatePiece     = this->Controller->GetLocalProcessId();
  this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
#else
  int maxpieces = 1;
#endif

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  outInfo->Set(vtkBoundsExtentTranslator::META_DATA(), this->ExtentTranslator);
  //
  //outInfo->Set(vtkBoundsExtentTranslator::META_DATA(),
  //             inInfo->Get(vtkBoundsExtentTranslator::META_DATA()));
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
               inInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()),
               6);

  if (outInfo->Has(vtkZoltanBasePartitionFilter::ZOLTAN_SAMPLE_ORIGIN())) {
    std::cout << "Zoltan received an origin request " << std::endl;
  };


  return 1;
}

//----------------------------------------------------------------------------
void vtkZoltanBasePartitionFilter::InitializeZoltanLoadBalance()
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
  Zoltan_Set_Param(this->ZoltanData, "LB_METHOD", "RCB");
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
//  Zoltan_Set_Param(this->ZoltanData, "NUM_GLOBAL_PARTS", global.str().c_str());
//  Zoltan_Set_Param(this->ZoltanData, "NUM_LOCAL_PARTS",  local.str().c_str());

  // if we have weights, turn on weight
  if(this->weights_data_ptr) {
      vtkDebugMacro("Setting zoltan weights dimension 1");
      Zoltan_Set_Param(this->ZoltanData, "OBJ_WEIGHT_DIM", "1");
  }
  else {
      vtkDebugMacro("Setting zoltan weights dimension 0");
      Zoltan_Set_Param(this->ZoltanData, "OBJ_WEIGHT_DIM", "0");
  }

  // we need the import and export lists
  Zoltan_Set_Param(this->ZoltanData, "RETURN_LISTS", "IMPORT AND EXPORT");

  // RCB parameters
  // Zoltan_Set_Param(this->ZoltanData, "PARMETIS_METHOD", "PARTKWAY");
  Zoltan_Set_Param(this->ZoltanData, "RCB_RECOMPUTE_BOX", "1");
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

  // Let Zoltan compute the load balance step, but we will manually
  // do the migrations of points/cells afterwards
  Zoltan_Set_Param(this->ZoltanData, "AUTO_MIGRATE", "0");
  Zoltan_Set_Param(this->ZoltanData, "MIGRATE_ONLY_PROC_CHANGES", "1");

  //
  // Query functions, to provide geometry to Zoltan
  //
  Zoltan_Set_Num_Obj_Fn(this->ZoltanData, get_number_of_objects_points,
      &this->ZoltanCallbackData);
  Zoltan_Set_Obj_List_Fn(this->ZoltanData, get_object_list_points,
      &this->ZoltanCallbackData);

  Zoltan_Set_Num_Geom_Fn(this->ZoltanData, get_num_geometry, &this->ZoltanCallbackData);
  if (this->ZoltanCallbackData.PointType == VTK_FLOAT) {
//    vtkDebugMacro("Using float data pointers ");
    Zoltan_Set_Geom_Multi_Fn(this->ZoltanData, get_geometry_list < float > ,
        &this->ZoltanCallbackData);
  } else if (this->ZoltanCallbackData.PointType == VTK_DOUBLE) {
//    vtkDebugMacro("Using double data pointers ");
    Zoltan_Set_Geom_Multi_Fn(this->ZoltanData, get_geometry_list < double > ,
        &this->ZoltanCallbackData);
  }
}

//----------------------------------------------------------------------------
int vtkZoltanBasePartitionFilter::PartitionPoints(vtkInformation*,
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
  this->UpdatePiece     = this->Controller->GetLocalProcessId();
  this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
  //
//  this->UpdatePiece     = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
//  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
//  int ghostLevel        = inInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
//  vtkDebugMacro("Partition filter " << this->UpdatePiece << " Ghost level " << ghostLevel);
  //
  Timer = vtkSmartPointer<vtkTimerLog>::New();
  Timer->StartTimer();

  // Get input
  vtkIdType       numPoints = input->GetNumberOfPoints();
  vtkIdType        numCells = input->GetNumberOfCells();
  vtkDataArray    *inPoints = numPoints>0 ? input->GetPoints()->GetData() : NULL;
  vtkDebugMacro("Partitioning { Points : " << numPoints << " } | { Cells " << numCells << " }");

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
  vtkSmartPointer<vtkPoints>   outPoints;
  if (this->UpdateNumPieces>1) {
    // create a new output points array to be filled
    outPoints = vtkSmartPointer<vtkPoints>::New();
    // if input had 0 points, make sure output is still setup correctly (float/double?)
    // collective exchanges will break if this is wrong as we may still receive data from another process
    // even though we are not sending any
    this->ZoltanCallbackData.PointType = this->GatherDataTypeInfo(input, input->GetPoints());
    outPoints->SetDataType(this->ZoltanCallbackData.PointType);
    output->SetPoints(outPoints);
  }
  else {
    // if only one process, we can just pass data through
    output->ShallowCopy(input);
    // vertex generation will fail if we don't set certain values
    this->ZoltanCallbackData.Output = output;
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
  this->ZoltanCallbackData.Input                    = input;
  this->ZoltanCallbackData.Output                   = output;
  this->ZoltanCallbackData.InputPointsData          = inPoints ? inPoints->GetVoidPointer(0) : NULL;
  this->ZoltanCallbackData.self                     = this;

  float ver;
  int zoltan_error = Zoltan_Initialize(0, NULL, &ver);
  if (zoltan_error != ZOLTAN_OK){
    printf("Zoltan initialization failed ...\n");
    return 0;
  }
  vtkDebugMacro("Zoltan Initialized");

  //
  // if a process has zero points, we need to make dummy data arrays to allow
  // space for when data gets sent in from other processes in the zoltan unpack function
  // This also stops hangs during collective operations by ensuring all ranks participate
  // we must not modfiy the input, so make a copy
  //
  vtkDebugMacro("Creating copy of input PointData")
  vtkPointData *inputPointData = input->GetPointData();
  this->ZoltanCallbackData.InputPointData = inputPointData->NewInstance();
  this->ZoltanCallbackData.InputPointData->ShallowCopy(inputPointData);
  this->AllocateFieldArrays(inputPointData, this->ZoltanCallbackData.InputPointData);
  vtkDebugMacro("FieldArrayPointers (point) Initialized");

  //
  // To convert Global to Local IDs we need to know the offsets of points/cells
  // relative to a zero based index on rank 0
  //
  this->ComputeIdOffsets(input->GetNumberOfPoints(), input->GetNumberOfCells());

  //
  // if weights are supplied, configure them
  //
  vtkDebugMacro("Setting up weights array");
  this->SetupPointWeights(this->ZoltanCallbackData.InputPointData);

  //
  // Set all the callbacks and user config parameters that will be used during the loadbalance
  //
  vtkDebugMacro("InitializeZoltanLoadBalance");
  this->InitializeZoltanLoadBalance();

  //
  // Check the input to see if it has a bounds translator already initialized
  // with partition info. If this is set, it tells us that the data has already been loadbalanced
  // somewhere else and we can skip main partitioning step. (example: H5Part reader can generate this)
  //
  vtkExtentTranslator *translator = inInfo ? vtkExtentTranslator::SafeDownCast(
      inInfo->Get(vtkBoundsExtentTranslator::META_DATA())) : NULL;
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
    this->ExecuteZoltanPartition(output, input);

    vtkDebugMacro("Partitioning "  <<
        " numImport : " << this->LoadBalanceData.numImport <<
        " numExport : " << this->LoadBalanceData.numExport
    );

    this->GetZoltanBoundingBoxes(globalBounds);

  }
  return 1;
}

//----------------------------------------------------------------------------
int vtkZoltanBasePartitionFilter::RequestData(vtkInformation* info,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // This is an abstract class and subclasses must override this function
  return 0;
}

//----------------------------------------------------------------------------
void vtkZoltanBasePartitionFilter::add_Id_to_interval_map(CallbackData *data, vtkIdType GID, vtkIdType LID) {
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
vtkIdType vtkZoltanBasePartitionFilter::global_to_local_Id(vtkIdType GID) {
  std::map<vtkIdType,vtkIdType>::iterator ub = this->ZoltanCallbackData.ReceivedGlobalToLocalIdMap.upper_bound(GID);
    ub--;
    vtkIdType key = ub->first;
    vtkIdType delta = ub->second;
    return GID - delta;
}

//----------------------------------------------------------------------------
MPI_Comm vtkZoltanBasePartitionFilter::GetMPIComm() {
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
void vtkZoltanBasePartitionFilter::InitializeFieldDataArrayPointers(
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
  callbackdata->NumberOfFields = outfielddata->GetNumberOfArrays();
  vtkDebugMacro("InitializeFieldDataArrayPointers "
                << ", in " <<infielddata->GetNumberOfArrays()
                << " out " << outfielddata->GetNumberOfArrays());
  for (int i=0; i<outfielddata->GetNumberOfArrays(); i++) {
    vtkDataArray *iarray = infielddata->GetArray(i);
    vtkDataArray *oarray = outfielddata->GetArray(i);
    oarray->Resize(Nfinal);
    oarray->SetNumberOfTuples(Nfinal);
    callbackdata->InputArrayPointers.push_back(iarray ? iarray->GetVoidPointer(0) : NULL);
    callbackdata->OutputArrayPointers.push_back(oarray->GetVoidPointer(0));
    // we need to know the amount of data to copy for each array tuple
    int Nc = oarray->GetNumberOfComponents();
    int Ns = oarray->GetDataTypeSize();
    callbackdata->MemoryPerTuple.push_back(Nc*Ns);
    callbackdata->TotalSizePerId += Nc*Ns;
    vtkDebugMacro("Got out array " << " " << oarray->GetName() << " sized " << Nfinal
        << " TotalSizePerId " << callbackdata->TotalSizePerId);
  }
  vtkDebugMacro("InitializeFieldDataArrayPointers completed");
}

//----------------------------------------------------------------------------
void vtkZoltanBasePartitionFilter::ComputeInvertLists(MigrationLists &migrationLists)
{
  //
  // Ask zoltan to create the inverse map of who sends/receives from who
  //
  int num_known                = static_cast<int>(migrationLists.known.nIDs ? migrationLists.known.nIDs : migrationLists.known.GlobalIds.size());
  ZOLTAN_ID_TYPE *GlobalIdsPtr = migrationLists.known.GlobalIdsPtr ? migrationLists.known.GlobalIdsPtr : (migrationLists.known.GlobalIds.size()>0 ? &migrationLists.known.GlobalIds[0] : NULL);
  int            *ProcsPtr     = migrationLists.known.ProcsPtr     ? migrationLists.known.ProcsPtr     : (migrationLists.known.Procs.size()>0     ? &migrationLists.known.Procs[0]     : NULL);
  //
  migrationLists.num_found        = 0;
  migrationLists.found_global_ids = NULL;
  migrationLists.found_local_ids  = NULL;
  migrationLists.found_procs      = NULL;
  migrationLists.found_to_part    = NULL;
  //
  int zoltan_error = Zoltan_Invert_Lists(this->ZoltanData,
    num_known,
    num_known>0 ? GlobalIdsPtr : NULL,
    NULL,
    num_known>0 ? ProcsPtr : NULL,
    NULL,
    &migrationLists.num_found,
    &migrationLists.found_global_ids,
    &migrationLists.found_local_ids,
    &migrationLists.found_procs,
    &migrationLists.found_to_part);

  if ( (zoltan_error != ZOLTAN_OK) )
  {
    printf("Zoltan_LB_Partition NOT OK\n");
    MPI_Finalize();
    Zoltan_Destroy(&this->ZoltanData);
    exit(0);
  }
  else if ( migrationLists.found_global_ids==NULL )
  {
    vtkDebugMacro("migrationLists.found_global_ids==NULL");
  }
  else if ( migrationLists.found_procs==NULL )
  {
    vtkDebugMacro("migrationLists.found_procs==NULL");
    MPI_Finalize();
    Zoltan_Destroy(&this->ZoltanData);
    exit(0);
  }
  else {
    vtkDebugMacro("migrationLists.found_procs or found_global_ids");
  }

  vtkDebugMacro("ComputeInvertLists "  <<
    " numImport : " << migrationLists.num_found <<
    " numExport : " << num_known
    );
}

//----------------------------------------------------------------------------
int vtkZoltanBasePartitionFilter::ManualPointMigrate(MigrationLists &migrationLists, bool keepinformation)
{

  int num_known                = static_cast<int>(migrationLists.known.nIDs ? migrationLists.known.nIDs : migrationLists.known.GlobalIds.size());
  ZOLTAN_ID_TYPE *GlobalIdsPtr = migrationLists.known.GlobalIdsPtr ? migrationLists.known.GlobalIdsPtr : (migrationLists.known.GlobalIds.size()>0 ? &migrationLists.known.GlobalIds[0] : NULL);
  int            *ProcsPtr     = migrationLists.known.ProcsPtr     ? migrationLists.known.ProcsPtr     : (migrationLists.known.Procs.size()>0     ? &migrationLists.known.Procs[0]     : NULL);

  //
  // After invert lists is called, we now know
  // 1) how many points we are sending
  // 2) how many points we are receiving
  // 3) how many points we are sending - but also keeping locally (computed locally earlier)
  // so now we can allocate just once the whole set of final arrays
  // instead of using a pre_migrate callback, we'll manually do it here
  // because we can use some local info and then delete some lists prior to the main exchange
  //
  int err;
  if (this->ZoltanCallbackData.PointType==VTK_FLOAT) {
    this->CopyPointsToSelf<float>(
      migrationLists.known.LocalIdsToKeep, migrationLists.num_reserved,
      &this->ZoltanCallbackData, 0, 0,
      migrationLists.num_found, migrationLists.found_global_ids, migrationLists.found_local_ids,
      migrationLists.found_procs, migrationLists.found_to_part,
      num_known,
      (num_known>0 ? GlobalIdsPtr : NULL),
      NULL,
      (num_known>0 ? ProcsPtr : NULL),
      NULL,
      &err
    );
  }
  else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
    this->CopyPointsToSelf<double>(
      migrationLists.known.LocalIdsToKeep, migrationLists.num_reserved,
      &this->ZoltanCallbackData, 0, 0,
      migrationLists.num_found, migrationLists.found_global_ids, migrationLists.found_local_ids,
      migrationLists.found_procs, migrationLists.found_to_part,
      num_known,
      (num_known>0 ? GlobalIdsPtr : NULL),
      NULL,
      (num_known>0 ? ProcsPtr : NULL),
      NULL,
      &err
    );
  }

  vtkDebugMacro("ManualPointMigrate (CopyPointsToSelf) ");
  return this->ZoltanPointMigrate(migrationLists, keepinformation);
}

//----------------------------------------------------------------------------
int vtkZoltanBasePartitionFilter::ZoltanPointMigrate(MigrationLists &migrationLists, bool keepinformation)
{
  int num_known                = static_cast<int>(migrationLists.known.nIDs ? migrationLists.known.nIDs : migrationLists.known.GlobalIds.size());
  ZOLTAN_ID_TYPE *GlobalIdsPtr = migrationLists.known.GlobalIdsPtr ? migrationLists.known.GlobalIdsPtr : (migrationLists.known.GlobalIds.size()>0 ? &migrationLists.known.GlobalIds[0] : NULL);
  int            *ProcsPtr     = migrationLists.known.ProcsPtr     ? migrationLists.known.ProcsPtr     : (migrationLists.known.Procs.size()>0     ? &migrationLists.known.Procs[0]     : NULL);

  //
  // Register functions for packing and unpacking data by migration tools.
  //
  if (this->ZoltanCallbackData.PointType==VTK_FLOAT) {
    zsize_fn  f1 = zoltan_obj_size_function_pointdata<float>;
    zpack_fn  f2 = zoltan_pack_obj_function_pointdata<float>;
    zupack_fn f3 = zoltan_unpack_obj_function_pointdata<float>;
    zprem_fn  f4 = zoltan_pre_migrate_function_null;
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &this->ZoltanCallbackData);
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &this->ZoltanCallbackData);
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &this->ZoltanCallbackData);
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &this->ZoltanCallbackData);
  }
  else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
    zsize_fn  f1 = zoltan_obj_size_function_pointdata<double>;
    zpack_fn  f2 = zoltan_pack_obj_function_pointdata<double>;
    zupack_fn f3 = zoltan_unpack_obj_function_pointdata<double>;
    zprem_fn  f4 = zoltan_pre_migrate_function_null;
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &this->ZoltanCallbackData);
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &this->ZoltanCallbackData);
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &this->ZoltanCallbackData);
    Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &this->ZoltanCallbackData);
  }

  CLEAR_ZOLTAN_DEBUG

  //
  // Now let zoltan perform the send/receive exchange of data
  //

  int zoltan_error = Zoltan_Migrate (this->ZoltanData,
    migrationLists.num_found,
    migrationLists.found_global_ids,
    migrationLists.found_local_ids,
    migrationLists.found_procs,
    migrationLists.found_to_part,
    num_known,
    (num_known>0 ? GlobalIdsPtr : NULL),
    NULL,
    (num_known>0 ? ProcsPtr : NULL),
    NULL
    );


#ifdef ZOLTAN_DEBUG_OUTPUT
    vtkDebugMacro("Partitioning complete " <<
      " pack_count : " << pack_count <<
      " size_count : " << size_count <<
      " unpack_count : " << unpack_count
     );
#endif

  int number_found = migrationLists.num_found;

  //
  // Release the arrays allocated during Zoltan_Invert_Lists
  //
  if (!keepinformation) {
    Zoltan_LB_Free_Part(
      &migrationLists.found_global_ids,
      &migrationLists.found_local_ids,
      &migrationLists.found_procs,
      &migrationLists.found_to_part);
    // set to zero so we know data has been deleted
    migrationLists.num_found = -1;
  }
  return number_found;
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
vtkSmartPointer<vtkPKdTree> vtkZoltanBasePartitionFilter::CreatePkdTree()
{
  vtkSmartPointer<vtkBSPCuts> cuts = vtkSmartPointer<vtkBSPCuts>::New();
  //
  RCB_STRUCT *rcb = (RCB_STRUCT *) (this->ZoltanData->LB.Data_Structure);
  struct rcb_tree *treept = rcb->Tree_Ptr;
  //
  if (treept[0].dim < 0) {
    vtkErrorMacro("RCB Tree invalid");
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

  return KdTree;
}


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int vtkZoltanBasePartitionFilter::zoltan_obj_size_function_pointdata(void *data,
  int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id,
  ZOLTAN_ID_PTR local_id, int *ierr)
{
  INC_SIZE_COUNT
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  *ierr = ZOLTAN_OK;
  return callbackdata->TotalSizePerId;
}

//----------------------------------------------------------------------------
// Zoltan callback to pack all the data for one point into a buffer
//----------------------------------------------------------------------------
void vtkZoltanBasePartitionFilter::zoltan_pack_obj_function_pointdata(void *data, int num_gid_entries, int num_lid_entries,
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
  *ierr = ZOLTAN_OK;
  return;
}

//----------------------------------------------------------------------------
// Zoltan callback to unpack all the data for one point from a buffer
//----------------------------------------------------------------------------
void vtkZoltanBasePartitionFilter::zoltan_unpack_obj_function_pointdata(void *data, int num_gid_entries,
  ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr)
{
  INC_UNPACK_COUNT
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  vtkIdType GID = *global_id;
  vtkIdType LID = GID - callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];
  //
  for (int i=0; i<callbackdata->NumberOfFields; i++) {
    int asize = callbackdata->MemoryPerTuple[i];
    char *dataptr = (char*)(callbackdata->OutputArrayPointers[i]) + asize*(callbackdata->MigrationPointCount);
    memcpy(dataptr, buf, asize);
    buf += asize;
  }
  callbackdata->MigrationPointCount++;
  *ierr = ZOLTAN_OK;
  return;
}

//----------------------------------------------------------------------------
void vtkZoltanBasePartitionFilter::zoltan_pre_migrate_function_pointdata(
  void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  //
  *ierr = ZOLTAN_OK;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
bool vtkZoltanBasePartitionFilter::MigratePointData(vtkDataSetAttributes *inPointData, vtkDataSetAttributes *outPointData)
{
  // we cannot operate if we don't have the right lists
  if (!this->KeepInversePointLists || this->MigrateLists.num_found == -1) {
    return false;
  }
  //
  // register new callback functions to pack/unpack data
  //
  zsize_fn  f1 = zoltan_obj_size_function_pointdata;
  zpack_fn  f2 = zoltan_pack_obj_function_pointdata;
  zupack_fn f3 = zoltan_unpack_obj_function_pointdata;
  zprem_fn  f4 = zoltan_pre_migrate_function_pointdata;
  Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_OBJ_SIZE_FN_TYPE,       (void (*)()) f1, &this->ZoltanCallbackData);
  Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PACK_OBJ_FN_TYPE,       (void (*)()) f2, &this->ZoltanCallbackData);
  Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_UNPACK_OBJ_FN_TYPE,     (void (*)()) f3, &this->ZoltanCallbackData);
  Zoltan_Set_Fn(this->ZoltanData, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, (void (*)()) f4, &this->ZoltanCallbackData);

  int N1 = this->ZoltanCallbackData.OutPointCount;
  // this sets internal flags used by CopyData to ensure arrays are marked for copying
  vtkDebugMacro("Setting up point data with "
      << inPointData->GetNumberOfArrays() << " "
      << outPointData->GetNumberOfArrays());
  outPointData->CopyAllocate(outPointData, N1);
  this->InitializeFieldDataArrayPointers(&this->ZoltanCallbackData, inPointData, outPointData, N1);

  vtkIdType N2 = inPointData->GetNumberOfTuples();

  this->ZoltanCallbackData.MigrationPointCount = 0;

  //
  // Pointdata for points which are not moving or being duplicated need to be preserved
  // before the migration step
  //
  vtkIdType maxID = 0, points_ok = 0;
  for (vtkIdType i=0; i<N2; i++) {
    // for each point that is staying on this process
    if (this->ZoltanCallbackData.LocalToLocalIdMap[i]>=0) {
      vtkIdType newID = this->ZoltanCallbackData.LocalToLocalIdMap[i];
      outPointData->CopyData(inPointData, i, newID);
      this->ZoltanCallbackData.MigrationPointCount++;
      maxID = std::max(maxID, newID);
    }
  }
  if (maxID>this->ZoltanCallbackData.MigrationPointCount) {
    vtkErrorMacro("Local ID mapped to new ID outside of permitted range");
  }

  //
  // perform a migration
  //
  CLEAR_ZOLTAN_DEBUG

  //
  // Now let zoltan perform the send/receive exchange of data between processes
  //
  int num_known = this->MigrateLists.known.GlobalIds.size();
  int zoltan_error = Zoltan_Migrate (this->ZoltanData,
    this->MigrateLists.num_found,
    this->MigrateLists.found_global_ids,
    this->MigrateLists.found_local_ids,
    this->MigrateLists.found_procs,
    this->MigrateLists.found_to_part,
    num_known,
    num_known>0 ? &this->MigrateLists.known.GlobalIds[0] : NULL,
    NULL,
    num_known>0 ? &this->MigrateLists.known.Procs[0]     : NULL,
    NULL
    );

#ifdef ZOLTAN_DEBUG_OUTPUT
    vtkDebugMacro("MigratePointData complete on " << this->UpdatePiece <<
      " pack_count : " << pack_count <<
      " size_count : " << size_count <<
      " unpack_count : " << unpack_count
     );
#endif

  vtkDebugMacro( "Expected " << N1 << " Points , found " << this->ZoltanCallbackData.MigrationPointCount);

  return true;
}

//----------------------------------------------------------------------------
void vtkZoltanBasePartitionFilter::AddHaloToBoundingBoxes(double GhostCellOverlap)
{
  //
  // Set the halo/ghost regions we need around each process bounding box
  //
  this->BoxListWithHalo.clear();
  if (this->InputExtentTranslator && this->InputExtentTranslator->GetBoundsHalosEnabled()) {
    this->ExtentTranslator->SetBoundsHalosEnabled(1);
    for (int p=0; p<this->UpdateNumPieces; p++) {
      vtkBoundingBox box;
      box.SetBounds(this->InputExtentTranslator->GetBoundsHaloForPiece(p));
      this->BoxListWithHalo.push_back(box);
      this->ExtentTranslator->SetBoundsHaloForPiece(p,this->InputExtentTranslator->GetBoundsHaloForPiece(p));
    }
  }
  else {
    this->ExtentTranslator->SetBoundsHalosEnabled(1);
    // @todo : extend this to handle AMR ghost regions etc.
    for (int p=0; p<this->UpdateNumPieces; p++) {
      vtkBoundingBox box = this->BoxList[p];
      box.Inflate(GhostCellOverlap);
      this->BoxListWithHalo.push_back(box);
      this->ExtentTranslator->SetBoundsHaloForPiece(p, box);
    }
  }
}
