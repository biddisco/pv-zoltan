/*=========================================================================

Program:   Visualization Toolkit
Module:    vtkDepthSortPainter.cxx

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// .NAME vtkDepthSortPainter
// .SECTION Thanks
// <verbatim>
//
//  This file is part of the PointSprites plugin developed and contributed by
//
//  Copyright (c) CSCS - Swiss National Supercomputing Centre
//                EDF - Electricite de France
//
//  John Biddiscombe, Ugo Varetto (CSCS)
//  Stephane Ploix (EDF)
//
// </verbatim>

#include "vtkDepthSortPainter.h"
//
#include "vtkObjectFactory.h"
#include "vtkGarbageCollector.h"
#include "vtkDataSet.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkPolyData.h"
#include "vtkCompositeDataSet.h"
#include "vtkCompositeDataIterator.h"
#include "vtkProperty.h"
#include "vtkDepthSortPolyData2.h"

//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkDepthSortPainter)
//-----------------------------------------------------------------------------
vtkCxxSetObjectMacro(vtkDepthSortPainter,DepthSortPolyData,vtkDepthSortPolyData2);
vtkCxxSetObjectMacro(vtkDepthSortPainter,OutputData,vtkDataObject);
//-----------------------------------------------------------------------------
vtkDepthSortPainter::vtkDepthSortPainter()
{
  this->DepthSortEnableMode             = ENABLE_SORT_IF_NO_DEPTH_PEELING;
  this->DepthSortMode                   = VTK_SORT_FIRST_POINT;
  this->UseCachedSortOrder              = true;
  this->Direction                       = VTK_DIRECTION_BACK_TO_FRONT;
  this->DepthSortPolyData               = vtkDepthSortPolyData2::New();
  this->OutputData                      = NULL;
  this->DepthSortOverrideFlag           = 0;
}
//-----------------------------------------------------------------------------
vtkDepthSortPainter::~vtkDepthSortPainter()
{
  this->SetDepthSortPolyData(NULL);
  this->SetOutputData(NULL);
}
//----------------------------------------------------------------------------
void vtkDepthSortPainter::ReportReferences(vtkGarbageCollector *collector)
{
  this->Superclass::ReportReferences(collector);

  //  vtkGarbageCollectorReport(collector, this->DepthSortPolyData, "DepthSortPolyData");
  //  vtkGarbageCollectorReport(collector, this->OutputData,        "OutputData");
}
//-----------------------------------------------------------------------------
void vtkDepthSortPainter::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
//-----------------------------------------------------------------------------
void vtkDepthSortPainter::SetDepthSortRequired(int required)
{
  this->DepthSortOverrideFlag = required; 
}
//-----------------------------------------------------------------------------
void vtkDepthSortPainter::PrepareForRendering(vtkRenderer* renderer, vtkActor* actor)
{
  // first set the DepthSortPolyData ivars
  if (this->DepthSortPolyData != NULL) {
    this->DepthSortPolyData->SetCamera(renderer->GetActiveCamera());
    this->DepthSortPolyData->SetProp3D(actor);
  }

  // check if we need to update
  if (this->GetMTime() < this->SortTime && this->DepthSortPolyData->GetMTime()
    < this->SortTime && this->GetInput()->GetMTime() < this->SortTime)
  {
    return;
  }

  // update the OutputData, initialize it with a shallow copy of the input
  this->SetOutputData(NULL);
  vtkDataObject * input = this->GetInput();
  if (!input) {
    return;
  }
  vtkDataObject* output = input->NewInstance();
  output->ShallowCopy(input);
  this->SetOutputData(output);
  output->FastDelete();

  if (this->DepthSortPolyData != NULL && this->NeedSorting(renderer, actor))
  {
    if (input->IsA("vtkCompositeDataSet")) {
      vtkCompositeDataSet* cdInput = vtkCompositeDataSet::SafeDownCast(input);
      vtkCompositeDataSet* cdOutput = vtkCompositeDataSet::SafeDownCast(this->OutputData);
      vtkCompositeDataIterator* iter = cdInput->NewIterator();
      for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem()) {
        vtkDataSet* pdInput = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
        vtkDataSet* pdOutput = vtkDataSet::SafeDownCast(cdOutput->GetDataSet(iter));
        if (pdInput && pdOutput) {
          this->Sort(pdOutput, pdInput, renderer, actor);
        }
      }
      iter->Delete();
    }
    else {
      this->Sort(vtkDataSet::SafeDownCast(this->OutputData),
        vtkDataSet::SafeDownCast(input), renderer, actor);
    }
    this->SortTime.Modified();
  }
}
//-----------------------------------------------------------------------------
void vtkDepthSortPainter::Sort(vtkDataSet* output,
  vtkDataSet* input,
  vtkRenderer* vtkNotUsed(renderer),
  vtkActor* vtkNotUsed(actor))
{
  this->DepthSortPolyData->SetInputData(input);
  this->DepthSortPolyData->SetDirectionToBackToFront();
  this->DepthSortPolyData->SetDepthSortMode(this->DepthSortMode);
  this->DepthSortPolyData->SetUseCachedSortOrder(this->UseCachedSortOrder);
  this->DepthSortPolyData->SetDirection(this->Direction);
  this->DepthSortPolyData->Update();
  output->ShallowCopy(this->DepthSortPolyData->GetOutput());
}
//-----------------------------------------------------------------------------
int vtkDepthSortPainter::NeedSorting(vtkRenderer* renderer, vtkActor* actor)
{
  // exit immediately if invalid or disabled 
  if (!actor || !renderer || (this->GetDepthSortEnableMode() == ENABLE_SORT_NEVER))
    return false;

  // hopefully not used very often
  if (this->GetDepthSortEnableMode() == ENABLE_SORT_ALWAYS)
    return true;

  // default mode
  if (this->GetDepthSortEnableMode() == ENABLE_SORT_IF_NO_DEPTH_PEELING
    && renderer->GetUseDepthPeeling())
    return false;

  // this flag might be set by the Representation if a vtkTwoScalarsToColors painter 
  // is being used -  it saves us checking color arrays by hand
  if (this->DepthSortOverrideFlag) 
    return true;

  if (actor->GetProperty()->GetOpacity() < 1)
    return true;

  // if the color table is translucent, then whoever is calling
  // this painter, should set DepthSortOverride flag to true

  return actor->HasTranslucentPolygonalGeometry();
}
//-----------------------------------------------------------------------------
vtkDataObject* vtkDepthSortPainter::GetOutput()
{
  return vtkDataObject::SafeDownCast(this->OutputData);
}
//-----------------------------------------------------------------------------
