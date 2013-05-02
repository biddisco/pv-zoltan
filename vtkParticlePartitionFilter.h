/*=========================================================================

  Module                  : vtkParticlePartitionFilter.h

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
// .NAME vtkParticlePartitionFilter Efficiently distribute particle datasets in parallel
// .SECTION Description
// vtkParticlePartitionFilter is a parallel load balancing/partitioning 
// filter for particle datasets. Halo regions may also be requested and particles
// will be duplicated on bordering processes of these regions.
// It uses the Zoltan library from the Trilinos package to perform the redistribution.

#ifndef __vtkParticlePartitionFilter_h
#define __vtkParticlePartitionFilter_h

#include "vtkZoltanV1PartitionFilter.h" // superclass
#include "vtkBoundingBox.h"
#include <vector>
#include "zoltan.h"

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
    vtkTypeMacro(vtkParticlePartitionFilter,vtkZoltanV1PartitionFilter);

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

  protected:
     vtkParticlePartitionFilter();
    ~vtkParticlePartitionFilter();

    // Description:
    virtual int RequestData(vtkInformation*,
                            vtkInformationVector**,
                            vtkInformationVector*);

    void AddHaloToBoundingBoxes();

    void FindPointsInHaloRegions(vtkPoints *pts, PartitionInfo &point_partitioninfo, ZoltanLoadBalanceData &loadBalanceData, PartitionInfo &ghost_info);

    double                      GhostCellOverlap;
    std::vector<vtkBoundingBox> BoxListWithHalo;

  private:
    vtkParticlePartitionFilter(const vtkParticlePartitionFilter&);  // Not implemented.
    void operator=(const vtkParticlePartitionFilter&);  // Not implemented.
};

#endif
