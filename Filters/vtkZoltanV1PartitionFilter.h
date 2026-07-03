/*=========================================================================

  Module : vtkZoltanV1PartitionFilter.h

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
// .NAME vtkZoltanV1PartitionFilter Efficiently distribute datasets in parallel
// .SECTION Description
// vtkZoltanV1PartitionFilter is the abstract base class used for parallel 
// load balancing/partitioning using the Zoltan library from Trilinos.
//
// .SECTION See Also
// vtkParticlePartitionFilter, vtkMeshPartitionFilter
//
#ifndef __vtkZoltanV1PartitionFilter_h
#define __vtkZoltanV1PartitionFilter_h
//
#include "vtkZoltanBasePartitionFilter.h" 

//----------------------------------------------------------------------------
class VTK_EXPORT vtkZoltanV1PartitionFilter : public vtkZoltanBasePartitionFilter
{
  public:
    static vtkZoltanV1PartitionFilter *New();
    vtkTypeMacro(vtkZoltanV1PartitionFilter, vtkZoltanBasePartitionFilter);

  protected:
     vtkZoltanV1PartitionFilter();
    ~vtkZoltanV1PartitionFilter();

    virtual void ExecuteZoltanPartition(vtkPointSet *output, vtkPointSet *input);
    virtual void GetZoltanBoundingBoxes(vtkBoundingBox &globalBounds);

  private:
    vtkZoltanV1PartitionFilter(const vtkZoltanV1PartitionFilter&);  // Not implemented.
    void operator=(const vtkZoltanV1PartitionFilter&);  // Not implemented.
};

#endif
