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
// .NAME vtkParticlePartitionFilter distribute particle datasets in parallel
// .SECTION Description
// vtkParticlePartitionFilter is a parallel load balancing/partitioning 
// filter for particle datasets. It uses the Zoltan library from the Trilinos 
// package to perform the redistribution.

#ifndef __vtkParticlePartitionFilter_h
#define __vtkParticlePartitionFilter_h

#include "vtkZoltanV1PartitionFilter.h" // superclass
#include "vtkBoundingBox.h"
#include <vector>
#include "zoltan.h"
#include "vtkSmartPointer.h"

class vtkMultiProcessController;
class vtkPoints;
class vtkIdTypeArray;
class vtkIntArray;
class vtkBoundsExtentTranslator;
class vtkPointSet;

class VTK_EXPORT vtkParticlePartitionFilter : public vtkZoltanV1PartitionFilter
{
  public:
    static vtkParticlePartitionFilter *New();
    vtkTypeMacro(vtkParticlePartitionFilter,vtkDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // The thickness of the region between each partition that is used for 
    // ghost cell exchanges. Any particles within this overlap region of another
    // processor will be duplicated on neighbouring processors (possibly multiple times
    // at corner region overlaps)
    vtkSetMacro(GhostCellOverlap, double);
    vtkGetMacro(GhostCellOverlap, double);
        
    // Description:
    // Return the Bounding Box for a partition plus the extended region
    // all around the box where ghost cells might be required/present
    vtkBoundingBox *GetPartitionBoundingBoxWithHalo(int partition);

    void FillPartitionBoundingBoxWithHalo();

    template <typename T>
    static void zoltan_pre_migrate_func_halo(void *data, int num_gid_entries, int num_lid_entries,
      int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
      int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
      ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr);

    void FindPointsInHaloRegions(vtkPoints *pts, vtkIdTypeArray *IdArray, PartitionInfo &ghostinfo);

    int ExchangeHaloPoints(vtkInformation* info,
          vtkInformationVector** inputVector,
          vtkInformationVector* outputVector);

    double GhostCellOverlap;
    vtkBoundingBox             *LocalBoxHalo;
    std::vector<vtkBoundingBox> BoxListWithHalo;

  protected:
     vtkParticlePartitionFilter();
    ~vtkParticlePartitionFilter();

    // Description:
    virtual int RequestData(vtkInformation*,
                            vtkInformationVector**,
                            vtkInformationVector*);

  private:
    vtkParticlePartitionFilter(const vtkParticlePartitionFilter&);  // Not implemented.
    void operator=(const vtkParticlePartitionFilter&);  // Not implemented.
};

#endif
