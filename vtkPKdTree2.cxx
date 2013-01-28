/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPKdTree2.cxx

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

#include "vtkPKdTree2.h"
#include "vtkKdNode.h"
#include "vtkDataSet.h"
#include "vtkCellCenters.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGrid.h"
#include "vtkObjectFactory.h"
#include "vtkMultiProcessController.h"
#include "vtkSocketController.h"
#include "vtkTimerLog.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkIntArray.h"
#include "vtkIdList.h"
#include "vtkSubGroup.h"
#include "vtkCommand.h"
//
#include "vtkSmartPointer.h"
#include "vtkAppendPolyData.h"
#include "vtkOutlineSource.h"
#include "vtkBoundingBox.h"

#include <stack>
#include <algorithm>


//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPKdTree2);
//----------------------------------------------------------------------------
vtkPKdTree2::vtkPKdTree2()
{
}
//----------------------------------------------------------------------------
vtkPKdTree2::~vtkPKdTree2()
{
  this->Rank = 0;
  this->NumRanks = 1;
}
//----------------------------------------------------------------------------
void vtkPKdTree2::GenerateRepresentation(int level, vtkPolyData *pd)
{
  this->GenerateBoxes(level, pd);
}
//----------------------------------------------------------------------------
void vtkPKdTree2::GenerateBoxes(int level, vtkPolyData *pd)
{
  if ((level < 0) || (level > this->Level)) {
    level = this->Level;
  }
  //
  vtkSmartPointer<vtkAppendPolyData> polys = vtkSmartPointer<vtkAppendPolyData>::New();

  vtkSmartPointer<vtkIntArray> processIds = vtkSmartPointer<vtkIntArray>::New();
  processIds->SetName("ProcessId");

  vtkKdNode    *kd    = this->Top;
  double *min = kd->GetMinBounds();
  double *max = kd->GetMaxBounds();

  std::stack<vtkKdNode*> tree_stack;
  tree_stack.push(this->Top);
  //
  while (!tree_stack.empty()) {
    vtkKdNode *node = tree_stack.top();
    tree_stack.pop();
    //
    if (node->GetLeft()) {
      tree_stack.push(node->GetLeft());
      tree_stack.push(node->GetRight());
    }
    else {
      double bounds[6];
      node->GetBounds(bounds);
      //
      vtkSmartPointer<vtkOutlineSource> cube = vtkSmartPointer<vtkOutlineSource>::New();
      cube->SetBounds(bounds);
      cube->Update();
      polys->AddInputData(cube->GetOutput());
      //    for (int p=0; p<8; p++) processIds->InsertNextValue(i);
    }
  }

  polys->Update();
  pd->SetPoints(polys->GetOutput()->GetPoints());
  pd->SetLines(polys->GetOutput()->GetLines());
//  pd->GetPointData()->AddArray(processIds);
//  pd->GetPointData()->AddArray(quantity);
}
//----------------------------------------------------------------------------

