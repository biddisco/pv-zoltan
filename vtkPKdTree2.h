/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPKdTree2.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*----------------------------------------------------------------------------
 Copyright (c) Sandia Corporation
 See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.
----------------------------------------------------------------------------*/

// .NAME vtkPKdTree2 - Build a k-d tree decomposition of a list of points.
//
// .SECTION Description
//      Build, in parallel, a k-d tree decomposition of one or more
//      vtkDataSets distributed across processors.  We assume each
//      process has read in one portion of a large distributed data set.
//      When done, each process has access to the k-d tree structure,
//      can obtain information about which process contains
//      data for each spatial region, and can depth sort the spatial
//      regions.
//
//      This class can also assign spatial regions to processors, based
//      on one of several region assignment schemes.  By default
//      a contiguous, convex region is assigned to each process.  Several
//      queries return information about how many and what cells I have
//      that lie in a region assigned to another process.
//
// .SECTION See Also
//      vtkKdTree

#ifndef __vtkPKdTree2_h
#define __vtkPKdTree2_h

#include "vtkFiltersParallelModule.h" // For export macro
#include "vtkPKdTree.h"

class vtkMultiProcessController;
class vtkCommunicator;
class vtkSubGroup;
class vtkIntArray;
class vtkKdNode;

class VTK_EXPORT vtkPKdTree2 : public vtkPKdTree
{
public:
  vtkTypeMacro(vtkPKdTree2, vtkPKdTree);
  static vtkPKdTree2 *New();

  void BuildLocator(double *bounds, int *remapping, int numregions);

  // Description:
  // Create a polydata representation of the boundaries of
  // the k-d tree regions.  If level equals GetLevel(), the
  // leaf nodes are represented.
  virtual void GenerateRepresentation(int level, vtkPolyData *pd);
  void GenerateBoxes(int level, vtkPolyData *pd);

   vtkSetMacro(InflateFactor, double);
   vtkGetMacro(InflateFactor, double);
//  void SetRank(int r) { this->Rank=r; }
//  void SetNumberOfRanks(int nr) { this->NumRanks=nr; }

protected:

   vtkPKdTree2();
  ~vtkPKdTree2();

  double InflateFactor;
private:

  vtkPKdTree2(const vtkPKdTree2&); // Not implemented
  void operator=(const vtkPKdTree2&); // Not implemented
};

#endif
