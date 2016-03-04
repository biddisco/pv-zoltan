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
#include "vtkUnsignedCharArray.h"
#include "vtkBoundingBox.h"
#include "vtkMath.h"
#include "vtkPointLocator.h"
#include "vtkPKdTree.h"
#include "vtkCellTreeLocator.h"
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
#if defined ZOLTAN_DEBUG_OUTPUT && !defined VTK_WRAPPING_CXX

# undef vtkDebugMacro
# define vtkDebugMacro(msg)  \
   DebugSynchronized(this->UpdatePiece, this->UpdateNumPieces, this->Controller, msg);

# undef  vtkErrorMacro
# define vtkErrorMacro(a) vtkDebugMacro(a)
#endif
//----------------------------------------------------------------------------
namespace debug {
    template<typename T>
    void output(const std::string &name, const std::vector<T> &v)
    {
#ifdef ZOLTAN_DEBUG_OUTPUT
        std::cout << name.c_str() << "\t : {" << v.size() << "} : ";
        std::copy(std::begin(v), std::end(v), std::ostream_iterator<T>(std::cout, ", "));
        std::cout << "\n" << std::flush;
#endif
    }

    template <class T>
    struct PrintPair : public std::unary_function<T, void>
    {
        std::ostream &os;
        PrintPair(std::ostream &strm) : os(strm) {}

        void operator()(const T& elem) const {
            os << "{" << elem.first << ", " << elem.second << "}, ";
        }
    };

    template<>
    void output(const std::string &name, const std::vector<std::pair<vtkIdType, int>> &v)
    {
#ifdef ZOLTAN_DEBUG_OUTPUT
        std::cout << name.c_str() << "\t : {" << v.size() << "} : ";
        std::for_each(v.begin(), v.end(), PrintPair<std::pair<vtkIdType, int>>(std::cout));
        std::cout << "\n" << std::flush;
#endif
    }

    template<typename Iter>
    void output(const std::string &name, Iter begin, Iter end)
    {
#ifdef ZOLTAN_DEBUG_OUTPUT
        std::cout << name.c_str() << "\t : {" << std::distance(begin, end) << "} : ";
        std::copy(begin, end,
            std::ostream_iterator<typename std::iterator_traits<Iter>::value_type>(
                std::cout, ", "));
        std::cout << "\n" << std::flush;
#endif
    }

    template <typename T>
    void output_sync(const std::string &name, T data, int N, int n, vtkMultiProcessController *ctrlr) {
#ifdef ZOLTAN_DEBUG_OUTPUT
        ctrlr->Barrier();
        for (int i=0; i<N; i++) {
            if (i==n) {
                std::cout << "Process " << n << " : ";
                debug::output(name, data);
            }
            ctrlr->Barrier();
        }
#endif
    }

};


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
void vtkMeshPartitionFilter::zoltan_pre_migrate_function_cell(
  void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);

  // How many cells will we have at the end
  vtkIdType OutputNumberOfLocalCells = callbackdata->Input->GetNumberOfCells();

  // our size estimates of the final number of cells can be messed up because
  // the list of cells being sent away contains some which are sent to multiple
  // remote processes. We can't therefore use the size of this list to subtract away
  // from the original cells to get final list size.
  // We therefore first need to count the number of unique cells being sent away
  //
  // We will do that whilst we ...
  // Mark cells to be moved with -1 so we know which local cells will still be local after the exchange
  debug_2("num_export cells " << num_export);
  callbackdata->LocalToLocalCellMap.assign(OutputNumberOfLocalCells, 0);
  vtkIdType uniqueSends = 0;
  for (vtkIdType i=0; i<num_export; i++) {
    vtkIdType GID = export_global_ids[i];
    vtkIdType LID = GID - callbackdata->ProcessOffsetsCellId[callbackdata->ProcessRank];
    if (callbackdata->LocalToLocalCellMap[LID]==0) {
      callbackdata->LocalToLocalCellMap[LID] = -1;
      uniqueSends++;
    }
  }

  // now mark those we are keeping as local again
  for (vtkIdType i=0; i<callbackdata->LocalIdsToKeep.size(); i++) {
    vtkIdType LID = callbackdata->LocalIdsToKeep[i];
    if (callbackdata->LocalToLocalCellMap[LID]==-1) {
      callbackdata->LocalToLocalCellMap[LID] = 0;
      uniqueSends--;
    }
  }

  // final cell count will be
  vtkIdType OutputNumberOfFinalCells = OutputNumberOfLocalCells + num_import - uniqueSends;

  debug_2("Copying cells to self"
            <<  " FinalCells:"<<OutputNumberOfFinalCells
            << "\tLocalCells:"<<OutputNumberOfLocalCells
            << "\tnum_import:"<<num_import
            << "\tuniqueSends:"<<uniqueSends
            << "\tLocalIdsToKeep:"<<callbackdata->LocalIdsToKeep.size());

  callbackdata->OutputPointsData = callbackdata->Output->GetPoints()->GetData()->GetVoidPointer(0);

  // get the output dataset pointer
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(callbackdata->Input);
  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(callbackdata->Input);
  vtkPolyData         *pdata2 = vtkPolyData::SafeDownCast(callbackdata->Output);
  vtkUnstructuredGrid *udata2 = vtkUnstructuredGrid::SafeDownCast(callbackdata->Output);

  // we'll need enough space to handle the largest cells,
  // for PolyData we must create a cell array for each of points/lines/strips/polys
  // for UnstructuredGrid, we will manage CellArray and CellTypeArray ourselves
  callbackdata->MaxCellSize     = callbackdata->Input->GetMaxCellSize();
  callbackdata->OutputUnstructuredCellTypes = NULL;
  if (pdata) {

      if ((callbackdata->self->polydata_types & 1) == 1) {
          vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
          verts->Allocate(OutputNumberOfFinalCells * 2);
          pdata2->SetVerts(verts);
      }
      if ((callbackdata->self->polydata_types & 2) == 2) {
          vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
          lines->Allocate(OutputNumberOfFinalCells * (callbackdata->MaxCellSize + 1));
          pdata2->SetLines(lines);
      }
      if ((callbackdata->self->polydata_types & 4) == 4) {
          vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
          polys->Allocate(OutputNumberOfFinalCells * (callbackdata->MaxCellSize + 1));
          pdata2->SetPolys(polys);
      }
      if ((callbackdata->self->polydata_types & 8) == 8) {
          vtkSmartPointer<vtkCellArray> strips = vtkSmartPointer<vtkCellArray>::New();
          strips->Allocate(OutputNumberOfFinalCells * (callbackdata->MaxCellSize + 1));
          pdata2->SetStrips(strips);
      }
      debug_2("Poly data types is " << callbackdata->self->polydata_types);
  }
  else if (udata) {
      callbackdata->OutputUnstructuredCellArray = vtkSmartPointer<vtkCellArray>::New();
      callbackdata->OutputUnstructuredCellArray->Allocate(OutputNumberOfFinalCells*(callbackdata->MaxCellSize+1));
      callbackdata->OutputUnstructuredCellArray->InitTraversal();
      callbackdata->OutputUnstructuredCellTypes = new int[OutputNumberOfFinalCells];
  }

  //
  vtkCellData *outCD = callbackdata->Output->GetCellData();
  //
  debug_2("Setting up cell data with "
      << callbackdata->InputCellData->GetNumberOfArrays() << " "
      << outCD->GetNumberOfArrays());
  callbackdata->self->InitializeFieldDataArrayPointers(callbackdata, callbackdata->InputCellData, outCD, OutputNumberOfFinalCells);

  vtkMeshPartitionFilter *mpf = dynamic_cast<vtkMeshPartitionFilter*>(callbackdata->self);
  if (mpf->ghost_cell_out_rank) {
      mpf->ghost_cell_out_rank->SetNumberOfTuples(OutputNumberOfFinalCells);
      for (int i=0; i<OutputNumberOfFinalCells; ++i) {
          mpf->ghost_cell_out_rank->SetValue(i,mpf->UpdatePiece+1);
      }
  }

  //
  vtkIdType npts, *pts, newPts[32];
  //
  callbackdata->OutCellCount = 0;
  for (vtkIdType cellId=0; cellId<OutputNumberOfLocalCells; cellId++) {
      if (callbackdata->LocalToLocalCellMap[cellId]!=-1) {
          // copy cell data from old to new datasets
          outCD->CopyData(callbackdata->InputCellData, cellId, callbackdata->OutCellCount);
          // copy cell point Ids to new dataset,
          int ctype;
          if (pdata) {
              ctype = pdata->GetCellType(cellId);
              pdata->GetCellPoints(cellId, npts, pts);
          }
          else if (udata) {
              ctype = udata->GetCellType(cellId);
              udata->GetCellPoints(cellId, npts, pts);
          }
          for (int i=0; i<npts; i++) {
              if (callbackdata->LocalToLocalIdMap[pts[i]]!=-1) {
                  newPts[i] = callbackdata->LocalToLocalIdMap[pts[i]];
              }
              else {
                  error_2("cell " << cellId << " point " << i << " assignment " << pts[i]);
              }
          }
          if (pdata) {
              //debug_2("Inserting cell of type " << ctype << " npts " << npts);
              pdata2->InsertNextCell(ctype, npts, newPts);
          }
          else if (udata) {
              callbackdata->OutputUnstructuredCellArray->InsertNextCell(npts, newPts);
              callbackdata->OutputUnstructuredCellTypes[callbackdata->OutCellCount] = ctype;
          }
          callbackdata->LocalToLocalCellMap[cellId] = callbackdata->OutCellCount;
          callbackdata->OutCellCount++;
      }
  }
  debug_2("completed zoltan_pre_migrate_function_cell");
}
//----------------------------------------------------------------------------
// Zoltan callback function : returns size of each cell and all its data
//----------------------------------------------------------------------------
int vtkMeshPartitionFilter::zoltan_obj_size_function_cell(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
{
  INC_SIZE_COUNT
  vtkZoltanBasePartitionFilter::CallbackData *callbackdata = (vtkZoltanBasePartitionFilter::CallbackData*)data;
  *ierr = ZOLTAN_OK;
  // return the size of the cell data + number of points in cell + point Ids
  vtkIdType GID = *global_id;
  vtkIdType LID = GID - callbackdata->ProcessOffsetsCellId[callbackdata->ProcessRank];
  //
  vtkIdType npts, *pts;
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(callbackdata->Input);
  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(callbackdata->Input);
  //
  if (pdata) {
      pdata->GetCellPoints(LID, npts, pts);
  }
  else if (udata) {
      udata->GetCellPoints(LID, npts, pts);
  }
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
  vtkZoltanBasePartitionFilter::CallbackData *callbackdata = (vtkZoltanBasePartitionFilter::CallbackData*)data;
  vtkIdType GID = *global_id;
  vtkIdType LID = GID - callbackdata->ProcessOffsetsCellId[callbackdata->ProcessRank];
  //
  vtkIdType npts, *pts, newPts[32];
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(callbackdata->Input);
  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(callbackdata->Input);

  int  ghost_temp;
  int *ghost;
  vtkIntArray *ghost_cell_rank =
      dynamic_cast<vtkMeshPartitionFilter*>(callbackdata->self)->ghost_cell_rank;
  if (ghost_cell_rank) {
      ghost = ghost_cell_rank->GetPointer(0);
      ghost_temp = ghost[LID];
      if (ghost[LID] > 0) {
          if (ghost[LID] == (dest + 1)) {
              ghost[LID] = 0;
          }
      }
  }

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
  int ctype;
  if (pdata) {
      ctype = pdata->GetCellType(LID);
      pdata->GetCellPoints(LID, npts, pts);
  }
  else if (udata) {
      ctype = udata->GetCellType(LID);
      udata->GetCellPoints(LID, npts, pts);
  }

  // copy the number of points and cell type so we know what's been sent when we unpack
  newPts[0] = npts;
  newPts[1] = ctype;
  // and the points Ids converted to global Ids
  for (int i=0; i<npts; i++) {
    newPts[i+2] = pts[i] + callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];
  }
  memcpy(buf, newPts, sizeof(vtkIdType)*(npts+2));
  //debug_2("Sending cell of type " << ctype << " npts " << npts << " to " << dest);

  if (ghost_cell_rank) {
      ghost[LID] = ghost_temp;
  }
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
  vtkZoltanBasePartitionFilter::CallbackData *callbackdata = (vtkZoltanBasePartitionFilter::CallbackData*)data;
  //
  vtkIdType GID = *global_id;
//  vtkIdType LID = GID - callbackdata->ProcessOffsetsCellId[callbackdata->ProcessRank];
  //
  vtkIdType npts, *pts, ctype, newPts[32];
  vtkPolyData         *pdata = vtkPolyData::SafeDownCast(callbackdata->Input);
  vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(callbackdata->Input);
  vtkPolyData         *pdata2 = vtkPolyData::SafeDownCast(callbackdata->Output);
  vtkUnstructuredGrid *udata2 = vtkUnstructuredGrid::SafeDownCast(callbackdata->Output);
  //
  for (int i=0; i<callbackdata->NumberOfFields; i++) {
    int asize = callbackdata->MemoryPerTuple[i];
    //  debug_2("receiving with asize at " << asize);
    char *dataptr = (char*)(callbackdata->OutputArrayPointers[i]) + asize*(callbackdata->OutCellCount);
    memcpy(dataptr, buf, asize);
    buf += asize;
  }
  // we have received a cell with a list of Global point Ids,
  // We need to convert each of the global Ids to the local Ids they became
  // when they were copied into the local point list
  npts  =  reinterpret_cast<vtkIdType*>(buf)[0];
  ctype =  reinterpret_cast<vtkIdType*>(buf)[1];
  pts   = &reinterpret_cast<vtkIdType*>(buf)[2];
  //debug_2("Received cell of type " << ctype << " npts " << npts);

  for (int i=0; i<npts; i++) {
    newPts[i] = callbackdata->self->global_to_local_Id(pts[i]);
  }
  if (pdata) {
      //debug_2("Inserting cell of type " << ctype << " npts " << npts);
      pdata2->InsertNextCell(ctype, npts, newPts);
  }
  else if (udata) {
      callbackdata->OutputUnstructuredCellArray->InsertNextCell(npts, newPts);
      callbackdata->OutputUnstructuredCellTypes[callbackdata->OutCellCount] = ctype;
  }

  callbackdata->OutCellCount++;
  *ierr = ZOLTAN_OK;
  return;
}
//----------------------------------------------------------------------------
// vtkMeshPartitionFilter :: implementation
//----------------------------------------------------------------------------
vtkMeshPartitionFilter::vtkMeshPartitionFilter()
{
  this->GhostMode           = vtkMeshPartitionFilter::None;
  this->BoundaryMode        = vtkMeshPartitionFilter::First;
  this->NumberOfGhostLevels = 0;
  this->ghost_cell_rank     = NULL;
  this->ghost_cell_flags    = NULL;
  this->ghost_cell_out_rank = NULL;
  this->KeepGhostRankArray  = 0;
  //this->DebugOn();
}
//----------------------------------------------------------------------------
vtkMeshPartitionFilter::~vtkMeshPartitionFilter()
{
  if (this->ZoltanData) {
    Zoltan_Destroy(&this->ZoltanData);
    this->ZoltanData = NULL;
  }
  // free up SmartPointers
  this->ghost_cell_rank     = NULL;
  this->ghost_cell_flags    = NULL;
  this->ghost_cell_out_rank = NULL;
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
  // Calculate even distribution of points across processes
  // This step only performs the load balance analysis,
  // no actual sending of data takes place yet.
  //
  this->PartitionPoints(info, inputVector, outputVector);

  if (this->UpdateNumPieces==1) {
    // input has been copied to output during PartitionPoints
    return 1;
  }

  // Make sure output cell arrays are setup before ghost ranks are added
  // otherwise array ordering is incorrect when ranks with empty data are present
  vtkDebugMacro(<<"Setting up in/out cell data arrays");
  vtkCellData *inCD  = this->ZoltanCallbackData.Input->GetCellData();
  vtkCellData *outCD = this->ZoltanCallbackData.Output->GetCellData();

  // we must not modify the input, so copy the cell data
  this->ZoltanCallbackData.InputCellData = inCD->NewInstance();
  this->ZoltanCallbackData.InputCellData->ShallowCopy(inCD);

  // and make sure that all ranks have the same/correct arrays on the input
  this->AllocateFieldArrays(inCD, this->ZoltanCallbackData.InputCellData);

  // initialize the output with correct arrays
  outCD->CopyAllOn();
  outCD->CopyAllocate(this->ZoltanCallbackData.InputCellData);

  vtkDebugMacro(<<"Output cell data " << outCD->GetNumberOfArrays());

  // if we are generating ghost cells for the mesh, then we must allocate a new array
  // on the output to store the ghost cell information (level 0,1,2...N ) etc
  vtkIdType numCells = this->ZoltanCallbackData.Input->GetNumberOfCells();
  if (this->GhostMode==vtkMeshPartitionFilter::Boundary ||
      this->GhostMode==vtkMeshPartitionFilter::BoundingBox) {
      vtkDebugMacro("Adding vtkGhostRanks");
      // ghost_cell_rank stores our internal rank info about ghosts
      this->ghost_cell_rank = vtkSmartPointer<vtkIntArray>::New();
      this->ghost_cell_rank->SetName("vtkGhostRanks");
      this->ghost_cell_rank->SetNumberOfTuples(numCells);
      this->ZoltanCallbackData.InputCellData->AddArray(this->ghost_cell_rank);
      this->ghost_cell_out_rank = vtkSmartPointer<vtkIntArray>::New();
      this->ghost_cell_out_rank->SetName("vtkGhostRanks");
      outCD->AddArray(this->ghost_cell_out_rank);
      // make sure the internals know that there is an extra array present
      this->ZoltanCallbackData.NumberOfFields++;
  }
  else {
    vtkDebugMacro("No need for vtkGhostRanks :"
                  << " inCD " << inCD->GetNumberOfArrays()
                  << " outCD " << outCD->GetNumberOfArrays());
  }

  //
  // based on the point partition, decide which cells need to be sent away
  // sending some cells may imply sending a few extra points too
  //
  PartitionInfo cell_partitioninfo;

  vtkDebugMacro("Calling BuildCellToProcessList");
  if (this->ZoltanCallbackData.PointType==VTK_FLOAT) {
    this->BuildCellToProcessList<float>(
      this->ZoltanCallbackData.Input,
      cell_partitioninfo,       // lists of which cells to send to which process
      this->MigrateLists.known, // list of which points to send to which process
      this->LoadBalanceData     // the partition information generated during PartitionPoints
    );
  }
  else if (this->ZoltanCallbackData.PointType==VTK_DOUBLE) {
    this->BuildCellToProcessList<double>(
      this->ZoltanCallbackData.Input,
      cell_partitioninfo,       // lists of which cells to send to which process
      this->MigrateLists.known, // list of which points to send to which process
      this->LoadBalanceData     // the partition information generated during PartitionPoints
    );
  }

  //*****************************************************************
  // Free the arrays allocated by Zoltan_LB_Partition
  // before we do a manual migration.
  //*****************************************************************
  vtkDebugMacro("Freeing Zoltan LB arrays");
  // Zoltan_LB_Free_Part(&this->LoadBalanceData.importGlobalGids, &this->LoadBalanceData.importLocalGids, &this->LoadBalanceData.importProcs, &this->LoadBalanceData.importToPart);
  // Zoltan_LB_Free_Part(&this->LoadBalanceData.exportGlobalGids, &this->LoadBalanceData.exportLocalGids, &this->LoadBalanceData.exportProcs, &this->LoadBalanceData.exportToPart);

  //
  // Based on the original partition and our extra cell point allocations
  // perform the main point exchange between all processes
  //
  this->ComputeInvertLists(this->MigrateLists);

  this->ManualPointMigrate(this->MigrateLists, this->KeepInversePointLists==1);


  if (!this->KeepInversePointLists) {
    vtkDebugMacro("Release point exchange data");
    this->MigrateLists.known.GlobalIds.clear();
    this->MigrateLists.known.Procs.clear();
    this->MigrateLists.known.LocalIdsToKeep.clear();
  }

  if (this->InputDisposable) {
    vtkDebugMacro("Disposing of input points and point data");
    this->ZoltanCallbackData.Input->SetPoints(NULL);
    this->ZoltanCallbackData.Input->GetPointData()->Initialize();
  }

  // after deleting memory, add a barrier to let ranks free as much as possible before the next big allocation
  this->Controller->Barrier();

  //
  // Distribute cells based on the usage of the points already distributed
  //
  this->PartitionCells(cell_partitioninfo);

  //
  // Distribute cells based on the usage of the points already distributed
  //
  this->UnmarkInvalidGhostCells(this->ZoltanCallbackData.Output);

  //
  // build a tree of bounding boxes to use for rendering info/hints or other spatial tests
  //
#ifdef VTK_ZOLTAN1_PARTITION_FILTER
  vtkDebugMacro("Create KdTree");
  this->CreatePkdTree();
  this->ExtentTranslator->SetKdTree(this->GetKdtree());
#endif

  //*****************************************************************
  // Free the storage allocated for the Zoltan structure.
  //*****************************************************************
  if (!this->KeepInversePointLists) {
    vtkDebugMacro("Zoltan_Destroy");
    Zoltan_Destroy(&this->ZoltanData);
    this->ZoltanData = NULL;
  }

  this->Timer->StopTimer();
  vtkDebugMacro("Mesh partitioning : " << this->Timer->GetElapsedTime() << " seconds");
  return 1;
}
//----------------------------------------------------------------------------
int vtkMeshPartitionFilter::PartitionCells(PartitionInfo &cell_partitioninfo)
{
  vtkDebugMacro("Entering PartitionCells");

  //
  // now we have a map of cells to processId, so do a collective 'invert lists'
  // operation to compute the global exchange map of who sends cells to who
  //
  size_t        num_known = cell_partitioninfo.GlobalIds.size();
  int           num_found = 0;
  ZOLTAN_ID_PTR found_global_ids = NULL;
  ZOLTAN_ID_PTR found_local_ids  = NULL;
  int          *found_procs      = NULL;
  int          *found_to_part    = NULL;
  //
  vtkDebugMacro("About to invert lists (cell migration)");
  int zoltan_error = Zoltan_Invert_Lists(this->ZoltanData,
    (int)num_known,
    num_known>0 ? &cell_partitioninfo.GlobalIds[0] : NULL,
    NULL,
    num_known>0 ? &cell_partitioninfo.Procs[0]     : NULL,
    NULL,
    &num_found,
    &found_global_ids,
    &found_local_ids,
    &found_procs,
    &found_to_part);
  //
  vtkDebugMacro("After invert lists : num_export " << num_known
      << " num_found " << num_found);
  //
  if (zoltan_error != ZOLTAN_OK){
    printf("Zoltan_LB_Partition NOT OK...\n");
    MPI_Finalize();
    Zoltan_Destroy(&this->ZoltanData);
    exit(0);
  }

  this->ZoltanCallbackData.LocalIdsToKeep = std::move(cell_partitioninfo.LocalIdsToKeep);

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
  vtkDebugMacro("About to Zoltan_Migrate (cells)");
  zoltan_error = Zoltan_Migrate (this->ZoltanData,
    (int)num_found,
    found_global_ids,
    found_local_ids,
    found_procs,
    found_to_part,
    (int)num_known,
    num_known>0 ? &cell_partitioninfo.GlobalIds[0] : NULL,
    NULL,
    num_known>0 ? &cell_partitioninfo.Procs[0]     : NULL,
    NULL
    );

  //
  // Release the arrays allocated during Zoltan_Invert_Lists
  //
  vtkDebugMacro("About to Free Zoltan_Migrate (cells)");
  Zoltan_LB_Free_Part(
    &found_global_ids,
    &found_local_ids,
    &found_procs,
    &found_to_part);
  vtkDebugMacro("Done Migration (cells)");

  vtkDebugMacro("Release pre-invert arrays (cells)");
  cell_partitioninfo.GlobalIds.clear();
  cell_partitioninfo.Procs.clear();
  this->ZoltanCallbackData.LocalIdsToKeep.clear();

  // For UnstructuredGrids, we must put the cells into the actual output dataset
  // It happens at the end with UnstructuredGrid because we manaully set the cell
  // array/types, whereas with polydata we use InsertNextCell at the dataset layer.
  vtkUnstructuredGrid *udata2 = vtkUnstructuredGrid::SafeDownCast(ZoltanCallbackData.Output);
  if (udata2) {
    udata2->SetCells(ZoltanCallbackData.OutputUnstructuredCellTypes, ZoltanCallbackData.OutputUnstructuredCellArray);
  }

  return 1;
}

//----------------------------------------------------------------------------
//
// Build a list of cell to process Ids based on the already performed
// point redistribution. Points which are being sent away will be taking
// their cells with them. Caution, some points are shared between cells and may
// go to different processes. For these cells, we must send another copy
// of the points to the relevant process.
//
template <typename T>
void vtkMeshPartitionFilter::BuildCellToProcessList(
  vtkPointSet *data,
  PartitionInfo &cell_partitioninfo,
  PartitionInfo &point_partitioninfo,
  ZoltanLoadBalanceData &loadBalanceData)
{
    vtkIdType numPts = data->GetNumberOfPoints();
    vtkIdType numCells = data->GetNumberOfCells();

    // we will put the results of the cell tests in these arrays
    cell_partitioninfo.Procs.reserve(numCells/this->UpdateNumPieces);
    cell_partitioninfo.GlobalIds.reserve(numCells/this->UpdateNumPieces);

    // we know that some points on this process will be exported to remote processes
    // so build a point to process map to quickly lookup the process Id from the point Id
    // 1) initially set all points as local to this process
    std::vector<int> localId_to_process_map(numPts, this->UpdatePiece);
    // 2) loop over all to be exported and note the destination
    for (vtkIdType i=0; i<loadBalanceData.numExport; i++) {
        vtkIdType id = loadBalanceData.exportGlobalGids[i] -
            this->ZoltanCallbackData.ProcessOffsetsPointId[this->ZoltanCallbackData.ProcessRank];
        localId_to_process_map[id] = loadBalanceData.exportProcs[i];
    }

    // the list now has the mapping from local point ID to remote process,
    // now we must examine each cell and see if they can be sent to a remote process
    // or if they are split between processes, when we find cells split across (multiple) remote processes,
    // we may need to add some points to the send lists so that the whole cell is correctly represented.
    // This vector will store {Id,process} pairs whilst the list is being scanned
    typedef std::pair<vtkIdType, int> process_tuple;
    std::vector<process_tuple> process_vector;

    // some points are marked to be sent away, but we will also need to keep a copy
    // reserve a little space to get things going (1% of total exports to start)
    point_partitioninfo.LocalIdsToKeep.reserve(loadBalanceData.numExport/100);

    // we must handle PolyData and UnstructuredGrid slightly differently
    vtkIdType npts, *pts;
    vtkPolyData         *pdata = vtkPolyData::SafeDownCast(data);
    vtkUnstructuredGrid *udata = vtkUnstructuredGrid::SafeDownCast(data);
    if (!pdata && !udata) {
        vtkErrorMacro("Only PolyData and UnstructuredGrid supported so far");
    }
    // polydata requires a cell map (verts/lines/polys/strips) to be present before we traverse cells
    if (pdata) pdata->BuildCells();

    // a simple bitmask with one entry per process, used for each cell to count processes for the cell
    std::vector<unsigned int> process_flag(this->UpdateNumPieces,0);
    std::vector<unsigned int> process_ghost(this->UpdateNumPieces,0);

    // sending a cellid to a process
    std::vector< process_tuple > cellDestProcesses;
    cellDestProcesses.reserve(this->UpdateNumPieces);


    // before iterating over cells, compute BBoxes if we need them
    double bounds[6];
    data->GetBounds(bounds);
    vtkBoundingBox localBoundingBox(bounds);
    vtkSmartPointer<vtkCellTreeLocator> cell_tree;
    vtkSmartPointer<vtkUnsignedCharArray> ghost_possible;
    if (this->GhostMode==vtkMeshPartitionFilter::BoundingBox) {
        this->AddHaloToBoundingBoxes(this->GhostCellOverlap);
        //
        cell_tree = vtkSmartPointer<vtkCellTreeLocator>::New();
        cell_tree->SetCacheCellBounds(1);
        cell_tree->SetNumberOfCellsPerNode(32);
        cell_tree->SetMaxLevel(20);
        cell_tree->SetLazyEvaluation(0);
        cell_tree->SetAutomatic(0);
        cell_tree->SetDataSet(data);
        cell_tree->BuildLocator();
        ghost_possible = vtkSmartPointer<vtkUnsignedCharArray>::New();
        ghost_possible->SetNumberOfTuples(numCells);
        for (vtkIdType cellId=0; cellId<numCells; ++cellId) {
            ghost_possible->SetValue(cellId,0);
        }
        vtkSmartPointer<vtkIdList> Ids = vtkSmartPointer<vtkIdList>::New();
        for (vtkIdType p=0; p<this->UpdateNumPieces; ++p) {
            cell_tree->FindCellsWithinBounds(this->BoxListWithHalo[p], Ids);
            for (vtkIdType i=0; i<Ids->GetNumberOfIds(); ++i) {
                ghost_possible->SetValue(Ids->GetId(i),1);
            }
            Ids->Reset();
        }
    }

    //
    // for each cell, find if any/all points are designated as remote and cell needs to be sent away
    //
    for (vtkIdType cellId=0; cellId<numCells; ++cellId) {
        // mark non-ghost status initially
        if (this->GhostMode!=vtkMeshPartitionFilter::None) {
            this->ghost_cell_rank->SetValue(cellId, 0);
        }

        // get a pointer to the cell points
        if (pdata) { pdata->GetCellPoints(cellId, npts, pts); }
        else if (udata) { udata->GetCellPoints(cellId, npts, pts); }

        //
        // we need to examine all points of the cell and classify it : there are several possible actions
        //
        // 1) all points are local                     : keep cell
        // 2) all points on same remote process        : send cell to remote process
        // 3) some points are local, some remote       : send or keep cell depending on count, make sure any local points sent are kept locally too
        // 4) all points on different remote processes : send cell to one remote process, also add other points to send list for that process
        //
        // step 1: increment a counter for each rank to mark/mask processes receiving points
        process_flag.assign(this->UpdateNumPieces,0);
        int points_remote = 0;
        int process_count = 0;
        for (int j=0; j<npts; ++j) {
            int process = localId_to_process_map[pts[j]];
            // if point will be sent away
            if (process!=this->UpdatePiece)   points_remote++;
            // if dest process has not been counted, count it
            if (process_flag[process]++ == 0) process_count++;
        }

        // step 2: classify the cell based on where all the points are going
        int cellDestProcess = this->UpdatePiece;
        bool cell_being_sent = false;
        //
        CellStatus cellstatus = UNDEFINED;
        if (points_remote==0) {
            cellstatus = LOCAL;
        }
        else {
            // all points on the same remote process
            if (process_count==1) {
                cellstatus = SAME;
            }
            // some local, others all on remote process(es)
            else if (process_count>1 && points_remote<npts) {
                cellstatus = SPLIT;
            }
            // all on remote processes, but not all the same process
            else if (process_count>1 && points_remote==npts) {
                cellstatus = SCATTERED;
            }
            else {
                throw std::string("This should not be possible");
            }

            // step 3: decide which rank will become the destination process for this cell
            if (this->BoundaryMode==vtkMeshPartitionFilter::First) {
                cellDestProcess = localId_to_process_map[pts[0]];
            }
            else if (this->BoundaryMode==vtkMeshPartitionFilter::Most) {
                cellDestProcess = std::distance(process_flag.begin(),
                    std::max_element(process_flag.begin(), process_flag.end()));
            }
            else if (this->BoundaryMode==vtkMeshPartitionFilter::Centroid) {
                T centroid[3];
                this->FindCentroid<T>(npts, pts, &this->ZoltanCallbackData, centroid);
                cellDestProcess = this->FindProcessFromPoint<T>(centroid);
            }

            // step 4: mark cell for sending if necessary
            if (cellDestProcess!=this->UpdatePiece) {
                cell_being_sent = true;
                // cell is going to be sent away, add it to temporary send list
                cellDestProcesses.push_back( process_tuple(cellId, cellDestProcess) );
            }
        }

        // step 5: any *_points_* which belong to a cell marked for sending, but are not already
        // marked for sending themselves, must be marked for both local and remote handling
        if (cellstatus==SAME && this->GhostMode==vtkMeshPartitionFilter::Boundary) {
            // all points are marked as remote, to the same destination
            // so nothing else needs to be done with the points
            continue;
        }
        //
        // cells of type SPLIT and SCATTERED need special treatment to handle
        // the cell split over multiple processes
        //
        else {
            // Ghost cell duplication : any SPLIT or SCATTERED cell is a boundary cell
            bool ghost_cell = false;
            process_ghost.assign(this->UpdateNumPieces,0);
            // In boundary mode, all SPLIT and SCATTERED cells are ghost cells
            if (this->GhostMode==vtkMeshPartitionFilter::Boundary &&
                (cellstatus==SPLIT || cellstatus==SCATTERED))
            {
                // cells of type SPLIT/SCATTERED must be duplicated on all processes receiving points
                ghost_cell = true;
                for (int p=0; p<this->UpdateNumPieces; p++) {
                    if (cellDestProcess != p && process_flag[p]) {
                        process_ghost[p] = 1;
                    }
                }
                this->ghost_cell_rank->SetValue(cellId, cellDestProcess+1);
            }
            else if (this->GhostMode==vtkMeshPartitionFilter::BoundingBox) {
                if (ghost_possible->GetValue(cellId)!=0) {
                    for (int p=0; p<this->UpdateNumPieces; p++) {
                        if (cellDestProcess!=p) {
                            for (int j=0; process_ghost[p]==0 && j<npts; ++j) {
                                double *pt = data->GetPoint(pts[j]);
                                if (BoxListWithHalo[p].ContainsPoint(pt)) {
                                    ghost_cell = true;
                                    process_ghost[p] = 1;
                                    this->ghost_cell_rank->SetValue(cellId, cellDestProcess+1);
                                }
                            }
                        }
                    }
                }
            }

            // should we keep a copy of this cell
            bool cell_keep_copy = ghost_cell && (process_flag[this->UpdatePiece]!=0 || process_ghost[this->UpdatePiece]!=0);
            if (cell_being_sent && cell_keep_copy) {
                // keep a ghost copy of the cell locally as well as sending it
                cell_partitioninfo.LocalIdsToKeep.push_back(cellId);
            }

            bool cell_sent = false;

            // for each point in the cell
            for (int j=0; j<npts; ++j) {
                const vtkIdType &ptId = pts[j];

                // which process is this point assigned to
                const int &pointDestProcess = localId_to_process_map[ptId];
                // is it going away
                bool point_migrating = (pointDestProcess!=this->UpdatePiece);
                bool point_send = false;
                bool point_keep = false;

                // Note 1) if the cell is being sent away, but some points are going
                // to a different location, then send copies to the cell destination
                // to ensure it receives all points for the cell
                if (cell_being_sent && pointDestProcess!=cellDestProcess) {
                    process_vector.push_back( process_tuple(ptId, cellDestProcess) );
                    point_send = true;
                }

                // if the cell is partially on this process and partially elsewhere
                if (cellstatus==SPLIT || cellstatus==SCATTERED) {
                    // if the cell is staying, but this point exported, keep a local copy
                    if (!cell_being_sent && point_migrating) {
                        point_partitioninfo.LocalIdsToKeep.push_back(ptId);
                        point_keep = true;
                    }
                    // if the cell is going, but this point is staying, keep a copy
                    // (it will have been marked for export above in Note 1)
                    else if (cell_being_sent && !point_migrating) {
                        point_partitioninfo.LocalIdsToKeep.push_back(ptId);
                        point_keep = true;
                    }
                }

                if (ghost_cell) {
                    for (int p=0; p<this->UpdateNumPieces; ++p) {
                        if (process_ghost[p]!=0) {
                            if (p!=pointDestProcess && p!=this->UpdatePiece) {
                                process_vector.push_back( process_tuple(ptId, p) );
                                point_send = true;
                            }
                            if (!point_keep) {
                                if (point_migrating) {
                                    if (p==this->UpdatePiece) {
                                        point_partitioninfo.LocalIdsToKeep.push_back(ptId);
                                    }
                                }
                                if (!point_migrating) {
                                    if (point_send && !cell_being_sent && pointDestProcess==this->UpdatePiece) {
                                        point_partitioninfo.LocalIdsToKeep.push_back(ptId);
                                    }
                                }
                            }
                            // send the cell if it wasn't already sent
                            if (p!=this->UpdatePiece) cellDestProcesses.push_back( process_tuple(cellId, p) );
                            if (!cell_being_sent) {
                                cell_partitioninfo.LocalIdsToKeep.push_back(cellId);
                            }
                        }
                    }
                }
            }
        }
    }

    vtkDebugMacro("Cleaning up the cell_partitioninfo struct ");
    // add the flagged cells to the final migration send list
    if (cellDestProcesses.size()>0) {
        std::sort(cellDestProcesses.begin(), cellDestProcesses.end());
        cellDestProcesses.resize(
            std::unique(cellDestProcesses.begin(), cellDestProcesses.end()) -
            cellDestProcesses.begin());
        for (auto &p : cellDestProcesses) {
            cell_partitioninfo.Procs.push_back(p.second);
            cell_partitioninfo.GlobalIds.push_back(
                p.first + this->ZoltanCallbackData.ProcessOffsetsCellId[this->UpdatePiece]);
        }
        cellDestProcesses.clear();
    }

    //
    // remove any duplicated export ids (note that these are pairs, so equality is both of <id,process>
    // some points might be sent to multiple locations, that's allowed.
    //
    std::sort(process_vector.begin(), process_vector.end());
    process_vector.resize(std::unique(process_vector.begin(), process_vector.end()) - process_vector.begin());
    //debug::output_sync("process vector", process_vector, this->UpdateNumPieces, this->UpdatePiece, this->Controller);
    //
    // remove any duplicated ids of points we are keeping as well as sending, not pairs here
    //
    std::sort(point_partitioninfo.LocalIdsToKeep.begin(), point_partitioninfo.LocalIdsToKeep.end());
    point_partitioninfo.LocalIdsToKeep.resize(
        std::unique(point_partitioninfo.LocalIdsToKeep.begin(), point_partitioninfo.LocalIdsToKeep.end()) -
        point_partitioninfo.LocalIdsToKeep.begin());
    //debug::output_sync("point_partitioninfo LocalIdsToKeep", point_partitioninfo.LocalIdsToKeep, this->UpdateNumPieces, this->UpdatePiece, this->Controller);

    std::sort(cell_partitioninfo.LocalIdsToKeep.begin(), cell_partitioninfo.LocalIdsToKeep.end());
    cell_partitioninfo.LocalIdsToKeep.resize(
        std::unique(cell_partitioninfo.LocalIdsToKeep.begin(), cell_partitioninfo.LocalIdsToKeep.end()) -
        cell_partitioninfo.LocalIdsToKeep.begin());
    //debug::output_sync("cell_partitioninfo LocalIdsToKeep", cell_partitioninfo.LocalIdsToKeep, this->UpdateNumPieces, this->UpdatePiece, this->Controller);

    //
    // After examining cells, we found some points we need to send away,
    // we must add these to the points the original load balance already declared in export lists
    //
    point_partitioninfo.GlobalIds.reserve(loadBalanceData.numExport + process_vector.size());
    point_partitioninfo.Procs.reserve(loadBalanceData.numExport + process_vector.size());

    // 1) add the points from original zoltan load balance to send list
    for (int i=0; i<loadBalanceData.numExport; ++i) {
        point_partitioninfo.GlobalIds.push_back(loadBalanceData.exportGlobalGids[i]);
        point_partitioninfo.Procs.push_back(loadBalanceData.exportProcs[i]);
    }
    // @todo : delete the arrays from zoltan above

    // 2) add the points from cell tests just performed to send list
    for (std::vector<process_tuple>::iterator x=process_vector.begin(); x!=process_vector.end(); ++x)
    {
        point_partitioninfo.GlobalIds.push_back(x->first + this->ZoltanCallbackData.ProcessOffsetsPointId[this->UpdatePiece]);
        point_partitioninfo.Procs.push_back(x->second);
    }
    //debug::output_sync("point_partitioninfo GlobalIds", point_partitioninfo.GlobalIds, this->UpdateNumPieces, this->UpdatePiece, this->Controller);
    //debug::output_sync("point_partitioninfo Procs", point_partitioninfo.Procs, this->UpdateNumPieces, this->UpdatePiece, this->Controller);

    vtkDebugMacro("BuildCellToProcessList "  <<
        " numImport : " << this->LoadBalanceData.numImport <<
        " numExport : " << point_partitioninfo.GlobalIds .size()
    );
    //debug::output_sync("cell_partitioninfo GlobalIds", cell_partitioninfo.GlobalIds, this->UpdateNumPieces, this->UpdatePiece, this->Controller);
    //debug::output_sync("cell_partitioninfo Procs", cell_partitioninfo.Procs, this->UpdateNumPieces, this->UpdatePiece, this->Controller);
}

//----------------------------------------------------------------------------
void vtkMeshPartitionFilter::UnmarkInvalidGhostCells(vtkPointSet *outdata)
{
    if (!this->ghost_cell_rank) {
        return;
    }
    vtkIdType numCells = outdata->GetNumberOfCells();
    // ghost_cell_flags is passed downstream and holds the real ghost type flags
    this->ghost_cell_flags = vtkUnsignedCharArray::New();
    this->ghost_cell_flags->SetName("vtkGhostType");
    this->ghost_cell_flags->SetNumberOfTuples(numCells);
    //

    int *ghostrankout = this->ghost_cell_out_rank->GetPointer(0);

    vtkDebugMacro("LocalToLocalCellMap " << this->ZoltanCallbackData.LocalToLocalCellMap.size());
    for (vtkIdType cellId=0; cellId<this->ZoltanCallbackData.LocalToLocalCellMap.size(); ++cellId) {
        //
        vtkIdType LID = this->ZoltanCallbackData.LocalToLocalCellMap[cellId];
        if (LID==-1) continue;
        //
        ghostrankout[LID] = 0;
        //
        int ghost_rank = this->ghost_cell_rank->GetValue(cellId);
        if (ghost_rank!=this->UpdatePiece+1) {
            ghostrankout[LID] = ghost_rank;
            this->ghost_cell_flags->SetValue(LID, 4); // vtkDataSetAttributes::DUPLICATECELL);
        }
    }

    //
    vtkDebugMacro("Unmark ghost flags " << numCells
        << " of ghost_cell_flags " << this->ghost_cell_flags->GetNumberOfTuples()
        << " from ghost_cell_rank_out " << this->ghost_cell_out_rank->GetNumberOfTuples());
    for (int i=0; i<numCells; ++i) {
        int ghost_rank = this->ghost_cell_out_rank->GetValue(i);
        this->ghost_cell_flags->SetValue(i, rand()%2);
        if (ghost_rank>0) {
            if (ghost_rank!=this->UpdatePiece+1) {
                this->ghost_cell_flags->SetValue(i, 4); //vtkDataSetAttributes::DUPLICATECELL);
            }
            else {
                this->ghost_cell_flags->SetValue(i, 0);
            }
        }
        else {
            this->ghost_cell_flags->SetValue(i, 0);
        }
    }
    //
    // for each cell, unmark any cells marked as ghosts that stayed local
    //

    //
    outdata->GetCellData()->AddArray(this->ghost_cell_flags);
    if (!this->KeepGhostRankArray) {
        this->ZoltanCallbackData.Output->GetCellData()->RemoveArray("vtkGhostRanks");
    }
    vtkDebugMacro("Completed unmark invalid ghost cells ");
}

//----------------------------------------------------------------------------

