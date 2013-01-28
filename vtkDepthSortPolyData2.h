/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDepthSortPolyData2.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkDepthSortPolyData2 - sort poly data along camera view direction
// .SECTION Description
// vtkDepthSortPolyData2 rearranges the order of cells so that certain
// rendering operations (e.g., transparency or Painter's algorithms)
// generate correct results. To use this filter you must specify the
// direction vector along which to sort the cells. You can do this by
// specifying a camera and/or prop to define a view direction; or
// explicitly set a view direction.

// .SECTION Caveats
// The sort operation will not work well for long, thin primitives, or cells
// that intersect, overlap, or interpenetrate each other.

#ifndef __vtkDepthSortPolyData2_h
#define __vtkDepthSortPolyData2_h

#include "vtkDepthSortPolyData.h"
#include "vtkSmartPointer.h"
class vtkIdTypeArray;

class VTK_EXPORT vtkDepthSortPolyData2 : public vtkDepthSortPolyData
{
public:
  // Description:
  // Instantiate object.
  static vtkDepthSortPolyData2 *New();

  vtkTypeMacro(vtkDepthSortPolyData2,vtkDepthSortPolyData);

    // Description:
  // Set/Get the sort origin. This ivar only has effect if the sort
  // direction is set to SetDirectionToSpecifiedVector(). The sort occurs
  // in the direction of the vector, with this point specifying the
  // origin.
//  vtkSetVector3Macro(Origin,double);
//  vtkGetVectorMacro(Origin,double,3);

  // Description:
  // When set we ignore all cells other than the polygon array and cache
  // arrays so that repeated calls as the camera moves will reuse the
  // previously sorted array as the start point
  vtkSetMacro(FastPolygonMode, int);
  vtkGetMacro(FastPolygonMode, int);
  vtkBooleanMacro(FastPolygonMode, int);
  
  // Description:
  // When UseCachedSortOrder is set, the sorted depth order is reused
  // to traverse cells and the sort algorithm changes to insertion sort
  // instead of quicksort (or std::sort default sort algorithm).
  // The reason is that when the points are already nearly sorted, the insertion
  // sort algorithm is ver fast and for small camera movments this is beneficial.
  // Ideal usage is that when the camera movement is less than delta, use insertion
  // otherwise, use quicksort.
  vtkSetMacro(UseCachedSortOrder, int);
  vtkGetMacro(UseCachedSortOrder, int);
  vtkBooleanMacro(UseCachedSortOrder, int);

protected:
  vtkDepthSortPolyData2();
  ~vtkDepthSortPolyData2();
  //
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  //
  int                             FastPolygonMode;
  int                             UseCachedSortOrder;
  vtkTimeStamp                    LastSortTime;
  void                           *SortingList;
  vtkSmartPointer<vtkIdTypeArray> DepthOrder;

private:
  vtkDepthSortPolyData2(const vtkDepthSortPolyData2&);  // Not implemented.
  void operator=(const vtkDepthSortPolyData2&);  // Not implemented.
};

#endif
