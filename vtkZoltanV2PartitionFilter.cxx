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
#include "vtkZoltanV2PartitionFilter.h"
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
//#include "vtkZoltanBasePartitionFilter.txx"

#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>


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
vtkStandardNewMacro(vtkZoltanV2PartitionFilter);
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// vtkZoltanV2PartitionFilter :: implementation 
//----------------------------------------------------------------------------
vtkZoltanV2PartitionFilter::vtkZoltanV2PartitionFilter()
{
  this->ZoltanParams = new Teuchos::ParameterList;
}
//----------------------------------------------------------------------------
vtkZoltanV2PartitionFilter::~vtkZoltanV2PartitionFilter()
{
  if (this->ZoltanParams) delete this->ZoltanParams;
}

//----------------------------------------------------------------------------
void vtkZoltanV2PartitionFilter::InitializeZoltanLoadBalance()
{
    // TODO : To be removed later
#ifdef HAVE_ZOLTAN2_MPI
    int rank, nprocs;
    MPI_Comm_size(this->GetMPIComm(), &nprocs);
    MPI_Comm_rank(this->GetMPIComm(), &rank);
#else
    int rank=0, nprocs=1;
#endif

    // Zoltan 2 parameters
    double tolerance = 1.1;
    if (this->ZoltanParams) delete this->ZoltanParams;
    this->ZoltanParams = new Teuchos::ParameterList("test params");
    this->ZoltanParams->set("debug_level", "basic_status");
    this->ZoltanParams->set("debug_procs", "0");
    this->ZoltanParams->set("error_check_level", "debug_mode_assertions");
    this->ZoltanParams->set("compute_metrics", "true");
    this->ZoltanParams->set("algorithm", "multijagged");
    this->ZoltanParams->set("imbalance_tolerance", tolerance);
    this->ZoltanParams->set("num_global_parts", nprocs);
    this->ZoltanParams->set("bisection_num_test_cuts", 1);
    this->ZoltanParams->set("mj_keep_part_boxes", 1);

    /*
    Zoltan2 all parameters
    {
      "error_check_level" : string = basic_assertions
      "debug_level" : string = basic_status
      "timer_type" : string = no_timers
      "debug_output_stream" : string = cout
      "timer_output_stream" : string = cout
      "memory_output_stream" : string = cout
      "debug_output_file" : string = /dev/null
      "timer_output_file" : string = /dev/null
      "memory_output_file" : string = /dev/null
      "debug_procs" : string = 0
      "mj_parts" : string = 0
      "memory_procs" : string = 0
      "speed_versus_quality" : string = balance
      "memory_versus_speed" : string = balance
      "random_seed" : string = 0.5
      "order_method" : string = rcm
      "order_package" : string = amd
      "compute_metrics" : string = no
      "topology" : string =
      "randomize_input" : string = no
      "partitioning_objective" : string = balance_object_weight
      "imbalance_tolerance" : string = 1.1
      "num_global_parts" : string = 0
      "num_local_parts" : string = 0
      "partitioning_approach" : string = partition
      "objects_to_partition" : string = graph_vertices
      "model" : string = graph
      "algorithm" : string = random
      "rectilinear" : string = no
      "average_cuts" : string = no
      "bisection_num_test_cuts" : int = 1
      "symmetrize_input" : string = no
      "subset_graph" : string = no
      "mj_concurrent_part_count" : int = 1
      "mj_minimum_migration_imbalance" : string = 1.1
      "mj_migration_option" : int = 1
      "remap_parts" : string = no
      "mapping_type" : int = -1
      "mj_keep_part_boxes" : int = -1
      "mj_enable_rcb" : int = 0
      "mj_recursion_depth" : int = -1
      "color_method" : string = rcm
      "color_choice" : string = amd
    }

    */
    vtkZoltanBasePartitionFilter::InitializeZoltanLoadBalance();
}

//----------------------------------------------------------------------------
template <typename scalar_t>
struct vtkZoltan2Helper
{
    typedef Zoltan2::BasicUserTypes<scalar_t, globalId_t, localId_t> myTypes;
    typedef Zoltan2::BasicVectorAdapter<myTypes> inputAdapter_t;
    typedef typename inputAdapter_t::part_t part_t;
    typedef Zoltan2::PartitioningProblem<inputAdapter_t> result_type;

    static result_type *SolveZoltan2Partition(
        vtkDataArray *datarray, vtkIdType localCount,
        globalId_t *globalIds, const scalar_t *weightarray,
        vtkZoltanV2PartitionFilter *self, vtkPointSet *input)
    {
//#define ZERO_COPY_DATA
#ifdef ZERO_COPY_DATA
        // points are {x,y,z} in a single array
        const scalar_t *x = static_cast<const scalar_t *>(datarray->GetVoidPointer(0));
        const scalar_t *y = static_cast<const scalar_t *>(datarray->GetVoidPointer(1));
        const scalar_t *z = static_cast<const scalar_t *>(datarray->GetVoidPointer(2));
        const int stride = 3;
        debug_2("USING ZERO COPY ");
#else
        scalar_t *coords = new scalar_t[3*localCount];
        scalar_t *x = coords;
        scalar_t *y = x + localCount;
        scalar_t *z = y + localCount;
        const int stride = 1;
        const scalar_t *p = static_cast<const scalar_t *>(datarray->GetVoidPointer(0));
        for (int i = 0; i < localCount; ++i)
        {
            x[i] = p[i*3 + 0];
            y[i] = p[i*3 + 1];
            z[i] = p[i*3 + 2];
        }
#endif
        // create an adapter that will point to the correct coordinates
        inputAdapter_t *InputAdapter = nullptr;
        result_type *problem1 = nullptr;

        if (weightarray==NULL) {
            InputAdapter = new inputAdapter_t(localCount, globalIds, x, y, z, stride, stride, stride);
            problem1     = new result_type(InputAdapter, self->ZoltanParams);
        }
        else {
            // coordinates
            std::vector<const scalar_t *> coordVec = {x, y, z};
            std::vector<int> coordStrides = {stride, stride, stride};
            // weights
            std::vector<const scalar_t*> weightVec = {weightarray};
            std::vector<int> weightStrides = {1};

            InputAdapter = new inputAdapter_t(
                localCount, globalIds,
                coordVec, coordStrides,
                weightVec, weightStrides);
            problem1    = new Zoltan2::PartitioningProblem<inputAdapter_t>(InputAdapter, self->ZoltanParams);
        }
        // Solve the problem
        problem1->solve();

        // we may query this from outside the filter
        self->ImbalanceValue = problem1->getWeightImbalance();

        // get the solution object
        const Zoltan2::PartitioningSolution<inputAdapter_t> &solution1 = problem1->getSolution();

        int numExport = 0;
        const part_t *partd = solution1.getPartListView();
        for (int i = 0; i < localCount; ++i)
        {
            if (partd[i]!=self->UpdatePiece)
            {
                numExport++;
            }
        }

        unsigned int *exportLocalGids = new unsigned int[numExport];
        unsigned int *exportGlobalGids = new unsigned int[numExport];
        int *exportProcs = new int[numExport];
        int k = 0;
        int offset = self->ZoltanCallbackData.ProcessOffsetsPointId[self->ZoltanCallbackData.ProcessRank];
        for (int i = 0; i < localCount; ++i)
        {
            if (partd[i]!=self->UpdatePiece){
                exportLocalGids[k] = i;
                exportGlobalGids[k] = i + offset;
                exportProcs[k] = partd[i];
                k++;
            }
        }
        self->LoadBalanceData.numExport = numExport;
        self->LoadBalanceData.exportGlobalGids = exportGlobalGids;
        self->LoadBalanceData.exportLocalGids  = exportLocalGids;
        self->LoadBalanceData.exportProcs = exportProcs;
        self->LoadBalanceData.exportToPart = exportProcs;

        // Zoltan 2 bounding box code
        self->BoxList.clear();
        std::vector<Zoltan2::coordinateModelPartBox<scalar_t, part_t> > &boxView = solution1.getPartBoxesView();
        for (int i=0; i<boxView.size(); i++) {
            double bounds[6];
            scalar_t *minss = boxView[i].getlmins();
            scalar_t *maxss = boxView[i].getlmaxs();
            bounds[0] = minss[0];
            bounds[1] = maxss[0];
            bounds[2] = minss[1];
            bounds[3] = maxss[1];
            bounds[4] = minss[2];
            bounds[5] = maxss[2];
            vtkBoundingBox box(bounds);
            self->BoxList.push_back(box);
            self->ExtentTranslator->SetBoundsForPiece(i, bounds);
        }
#ifndef ZERO_COPY_DATA
        delete []coords;
#endif
        return problem1;
    }
};

//----------------------------------------------------------------------------
#undef vtkTemplateMacroCase
#define vtkTemplateMacroCase(typeN, type, call)     \
  case typeN: { typedef type VTK_TT; call; }; break
#define vtkZoltanTemplateMacro(call)                \
  vtkTemplateMacroCase(VTK_DOUBLE, double, call);   \
  vtkTemplateMacroCase(VTK_FLOAT, float, call);

//----------------------------------------------------------------------------
void vtkZoltanV2PartitionFilter::ExecuteZoltanPartition(
    vtkPointSet *output,
    vtkPointSet *input)
{
    //
    // Get the input arrays, points coordinates and weights
    //

    // get the array that is used for coordinates, it might be float or double
    vtkPoints    *myInPoints = input->GetPoints();
    vtkDataArray *coordArray = (myInPoints!=nullptr)?myInPoints->GetData():nullptr;

    //////////////////////////////////////////////////////////////////////
    // Zoltan 2 partitioning
    //////////////////////////////////////////////////////////////////////
    int rank = this->UpdatePiece;
    int nprocs = this->UpdateNumPieces;

    int localCount = (coordArray)?coordArray->GetNumberOfTuples():0;
    // global Ids should be optional, fix this when they are
    globalId_t *globalIds = new globalId_t [localCount];
    globalId_t offset = this->ZoltanCallbackData.ProcessOffsetsPointId[this->ZoltanCallbackData.ProcessRank];
    for (size_t i=0; i < localCount; i++)
        globalIds[i] = offset++;

    switch(/*coordArray->GetDataType()*/ VTK_FLOAT)
    {
        vtkZoltanTemplateMacro(
            vtkZoltan2Helper<VTK_TT>::SolveZoltan2Partition(
                coordArray, localCount, globalIds, static_cast<VTK_TT*>(this->weights_data_ptr), this, input));
    }


    //////////////////////////////////////////////////////////////////////
    // Zoltan 2 partitioning ends here
    //////////////////////////////////////////////////////////////////////

}

//----------------------------------------------------------------------------
void vtkZoltanV2PartitionFilter::GetZoltanBoundingBoxes(vtkBoundingBox &globalBounds)
{
/*
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
        // cout<<"###"<<" ("<<bounds[0]<<", "<<bounds[1]<<")\t"<<" ("<<bounds[2]<<", "<<bounds[3]<<")\t"<<" ("<<bounds[4]<<", "<<bounds[5]<<")\t";
      }
      // cout<<endl;
    }
*/
    this->ExtentTranslator->InitWholeBounds();
}
