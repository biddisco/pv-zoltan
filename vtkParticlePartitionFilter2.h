/*=========================================================================

  Module                  : vtkParticlePartitionFilter2.h

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
// .NAME vtkParticlePartitionFilter2 Efficiently distribute particle datasets in parallel
// .SECTION Description
// vtkParticlePartitionFilter2 is a parallel load balancing/partitioning 
// filter for particle datasets. Halo regions may also be requested and particles
// will be duplicated on bordering processes of these regions.
// It uses the Zoltan library from the Trilinos package to perform the redistribution.

#ifndef __vtkParticlePartitionFilter2_h
#define __vtkParticlePartitionFilter2_h

#include "vtkZoltanV2PartitionFilter.h" // superclass
#include "vtkBoundingBox.h"
#include <vector>

class vtkMultiProcessController;
class vtkPoints;
class vtkIdTypeArray;
class vtkIntArray;
class vtkBoundsExtentTranslator;
class vtkPointSet;

class VTK_EXPORT vtkParticlePartitionFilter2 : public vtkZoltanV2PartitionFilter
{
  public:
    static vtkParticlePartitionFilter2 *New();
    vtkTypeMacro(vtkParticlePartitionFilter2,vtkZoltanV2PartitionFilter);

    // Description:
    // The thickness of the region between each partition that is used for 
    // ghost cell exchanges. Any particles within this overlap region of another
    // processor will be duplicated on neighbouring processors (possibly multiple times
    // at corner region overlaps)
    vtkSetMacro(GhostCellOverlap, double);
    vtkGetMacro(GhostCellOverlap, double);
        
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

        
    // Description:
    // Return the Bounding Box for a partition plus the extended region
    // all around the box where ghost cells might be required/present
    vtkBoundingBox *GetPartitionBoundingBoxWithHalo(int partition);

  protected:
     vtkParticlePartitionFilter2();
    ~vtkParticlePartitionFilter2();

    // Description:
    virtual int RequestData(vtkInformation*,
                            vtkInformationVector**,
                            vtkInformationVector*);

    void AddHaloToBoundingBoxes();

    void FindPointsInHaloRegions(vtkPoints *pts, PartitionInfo &point_partitioninfo, ZoltanLoadBalanceData &loadBalanceData, PartitionInfo &ghost_info);

    double                      GhostCellOverlap;
    std::vector<vtkBoundingBox> BoxListWithHalo;
    double                      GridSpacing;
    double                      GridOrigin[3];

  private:
    vtkParticlePartitionFilter2(const vtkParticlePartitionFilter2&);  // Not implemented.
    void operator=(const vtkParticlePartitionFilter2&);  // Not implemented.
};

#endif
