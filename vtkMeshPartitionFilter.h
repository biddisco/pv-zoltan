/*=========================================================================

  Module                  : vtkMeshPartitionFilter.h

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
// .NAME vtkMeshPartitionFilter Efficiently distribute PolyData or UnstructuredGrid datasets in parallel
// .SECTION Description
// vtkMeshPartitionFilter is a parallel load balancing/partitioning 
// filter for PolyData or UnstructuredGrid datasets. 
// It uses the Zoltan library from the Trilinos package to perform the redistribution.

#ifndef __vtkMeshPartitionFilter_h
#define __vtkMeshPartitionFilter_h

#include "vtkZoltanVersion.h"

#include "vtkBoundingBox.h"
#include "vtkSmartPointer.h"
//
#include <vector>
//
#include "zoltan.h"

class vtkMultiProcessController;
class vtkPoints;
class vtkIdTypeArray;
class vtkUnsignedCharArray;
class vtkBoundsExtentTranslator;
class vtkPointSet;

class VTK_EXPORT vtkMeshPartitionFilter : public VTK_ZOLTAN_PARTITION_FILTER
{
  public:
    static vtkMeshPartitionFilter *New();
    vtkTypeMacro(vtkMeshPartitionFilter, VTK_ZOLTAN_PARTITION_FILTER);
    void PrintSelf(ostream& os, vtkIndent indent);

    template <typename T>
    static void zoltan_pre_migrate_function_cell(void *data, int num_gid_entries, int num_lid_entries,
      int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
      int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
      ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr);

    static int zoltan_obj_size_function_cell(void *data, int num_gid_entries, int num_lid_entries,
      ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr);

    static void zoltan_pack_obj_function_cell(void *data, int num_gid_entries, int num_lid_entries,
      ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size, char *buf, int *ierr);

    static void zoltan_unpack_obj_function_cell(void *data, int num_gid_entries,
      ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr);

    int PartitionCells(PartitionInfo &cell_partitioninfo);

    template <typename T>
    void BuildCellToProcessList(
      vtkDataSet *data, 
      PartitionInfo &cell_partitioninfo, 
      PartitionInfo &point_partitioninfo, 
      ZoltanLoadBalanceData &loadBalanceData);
  
    enum GhostAlgorithm {
        None        = 0,
        Boundary    = 1,
        BoundingBox = 2,
        Neighbour   = 3
    };

    // GhostMode is an option to control how ghost cells are generated, modes are
    // Boundary: flags only cells which straddle the boundary of a partition
    // the cell will be duplicated on all partitions that it overlaps, but on one
    // partition it will be a normal cell, on others it will be flagged as s ghost cell
    // BoundingBox: duplicates all cells which overlap the GhostCellOverlap region
    // of another process.
    // Neighbour: Any cell touching a boundary cell up to the GhostLevel depth
    // becomes a ghost cell on an adjacent process.
    vtkSetMacro(GhostMode, int);
    vtkGetMacro(GhostMode, int);
    // convenience setter/getters for GhostMode
    void SetGhostModeToNone() { this->SetGhostMode(None); }
    void SetGhostModeToBoundary() { this->SetGhostMode(Boundary); }
    void SetGhostModeToBoundingBox() { this->SetGhostMode(BoundingBox); }
    void SetGhostModeToNeighbourCells() { this->SetGhostMode(Neighbour); }

    // Description:
    // Specify the ghost level that will be used to generate ghost cells
    // a level of 1 produces one layer of touching cells, 2 produces 2 layers etc
    // Note that this variable will be overridden if the information key
    // for GHOST_LEVELS is present in the information passed upstream in
    // the pipeline.
    // When GhostMode is BoundingBox, then this value is ignored and the
    // GhostCellOverlap is used instead as a measure for ghost regions
    vtkSetMacro(NumberOfGhostLevels, int);
    vtkGetMacro(NumberOfGhostLevels, int);

    // The distance beyond a process region for which we require ghost cells
    vtkSetMacro(GhostCellOverlap, double);
    vtkGetMacro(GhostCellOverlap, double);

  protected:
     vtkMeshPartitionFilter();
    ~vtkMeshPartitionFilter();
  
    // Description:
    virtual int RequestData(vtkInformation*,
                            vtkInformationVector**,
                            vtkInformationVector*);

    int                     GhostMode;
    double                  GhostCellOverlap;
    int                     NumberOfGhostLevels;
    vtkUnsignedCharArray   *ghost_array;

  private:
    vtkMeshPartitionFilter(const vtkMeshPartitionFilter&);  // Not implemented.
    void operator=(const vtkMeshPartitionFilter&);  // Not implemented.
};

#endif
