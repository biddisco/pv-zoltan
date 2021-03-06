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

#include "vtkZoltanVersion.h"

#include "vtkBoundingBox.h"
#include <vector>

class vtkMultiProcessController;
class vtkPoints;
class vtkIdTypeArray;
class vtkIntArray;
class vtkBoundsExtentTranslator;
class vtkPointSet;

class VTK_EXPORT vtkParticlePartitionFilter : public @VTK_ZOLTAN_PARTITION_FILTER@
{
  public:
    static vtkParticlePartitionFilter *New();
    vtkTypeMacro(vtkParticlePartitionFilter, @VTK_ZOLTAN_PARTITION_FILTER@);

    // Description:
    // Specify the point spacing on the X/Y/Z axis
    vtkSetMacro(GridSpacing, double);
    vtkGetMacro(GridSpacing, double);

    // Description
    // If Resolution[X/Y/Z]are all Non-zero, then
    // the spacing is ignored and the box defined by the points is
    // sampled using the specified resolutions
    vtkSetVector3Macro(GridOrigin, double);
    vtkGetVector3Macro(GridOrigin, double);

  protected:
     vtkParticlePartitionFilter();
    ~vtkParticlePartitionFilter();

    // Description:
    virtual int RequestData(vtkInformation*,
                            vtkInformationVector**,
                            vtkInformationVector*);

    void FindPointsInHaloRegions(vtkPoints *pts, PartitionInfo &point_partitioninfo, ZoltanLoadBalanceData &loadBalanceData, PartitionInfo &ghost_info);

    double                      GridSpacing;
    double                      GridOrigin[3];

  private:
    vtkParticlePartitionFilter(const vtkParticlePartitionFilter&);  // Not implemented.
    void operator=(const vtkParticlePartitionFilter&);  // Not implemented.
};

#endif
