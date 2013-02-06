/*=========================================================================

Program:   ParaView
Module:    $RCSfile$

Copyright (c) Kitware, Inc.
All rights reserved.
See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkDepthSortRepresentation.h"

#include "vtkObjectFactory.h"
#include "vtkGarbageCollector.h"
#include "vtkCompositePolyDataMapper2.h"
#include "vtkDepthSortDefaultPainter.h"
#include "vtkPVLODActor.h"
#include "vtkTexture.h"
#include "vtkPVCacheKeeper.h"
#include "vtkPVView.h"
#include "vtkPVRenderView.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkAlgorithmOutput.h"
#include "vtkCommunicator.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkMultiProcessController.h"
#include "vtkDummyController.h"
#include "vtkOutputWindow.h"
#include "vtkPVGeometryFilter.h"
//
#include "vtkDepthSortPainter.h"
#include "vtkDepthSortPolygonsPainter.h"
#include "vtkTwoScalarsToColorsPainter.h"
#include "vtkBoundsExtentTranslator.h"
#include "vtkPistonPolygonsPainter.h"
//
#include <vtksys/SystemTools.hxx>
#include <vtksys/SystemInformation.hxx>
#include <vtksys/RegularExpression.hxx>
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkDepthSortRepresentation);
vtkCxxSetObjectMacro(vtkDepthSortRepresentation, Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
vtkDepthSortRepresentation::vtkDepthSortRepresentation()
{
  vtkOutputWindow *win = vtkOutputWindow::GetInstance();
  win->SetPromptUser(0);
  //
  vtkPVGeometryFilter::SafeDownCast(this->GeometryFilter)->SetUseOutline(0);
  //
  this->UseDataPartitions         = 1;
  this->DepthSortDefaultPainter   = vtkDepthSortDefaultPainter::New();
  this->DepthSortPainter          = this->DepthSortDefaultPainter->GetDepthSortPainter();
  this->TwoScalarsToColorsPainter = this->DepthSortDefaultPainter->GetTwoScalarsToColorsPainter();
  this->DepthSortPolygonsPainter  = this->DepthSortDefaultPainter->GetDepthSortPolygonsPainter();
  this->PistonPolygonsPainter     = this->DepthSortDefaultPainter->GetPistonPolygonsPainter();

  //
  vtkMath::UninitializeBounds(this->GlobalDataBounds);
  //
  vtkCompositePolyDataMapper2* compositeMapper = vtkCompositePolyDataMapper2::SafeDownCast(this->Mapper);
  this->DepthSortDefaultPainter->SetDelegatePainter(compositeMapper->GetPainter()->GetDelegatePainter());
  compositeMapper->SetPainter(this->DepthSortDefaultPainter);
  //
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
  if (this->Controller == NULL) {
    this->SetController(vtkSmartPointer<vtkDummyController>::New());
  }
}
//----------------------------------------------------------------------------
vtkDepthSortRepresentation::~vtkDepthSortRepresentation()
{
  this->DepthSortDefaultPainter->FastDelete();
  this->SetController(NULL);
}
//----------------------------------------------------------------------------
void vtkDepthSortRepresentation::ReportReferences(vtkGarbageCollector *collector)
{
  this->Superclass::ReportReferences(collector);
}
//----------------------------------------------------------------------------
bool vtkDepthSortRepresentation::AddToView(vtkView* view)
{
  vtkPVRenderView* rview = vtkPVRenderView::SafeDownCast(view);
  if (rview && !vtkPistonPolygonsPainter::IsEnabledCudaGL()) {
    
    // if we have no hints as to the display VAR, we'll use the rank (modulus GPU count)
    vtkMultiProcessController *controller = vtkMultiProcessController::GetGlobalController();
    int rank = controller->GetLocalProcessId();

    // see if a DISPLAY env var was set for us to get the GPU from
    // match DISPLAY vars of the form :0.0, :0.1, :0.2, :0.N - and get N
    std::string display;
    vtksys::SystemTools::GetEnv("DISPLAY", display);
    vtksys::RegularExpression regex(".*:.*\\.([0-9]+)");
    regex.find(display.c_str());
    int displaynum = -1;
    if (regex.match(1).size()>0) {
      displaynum = atoi(regex.match(1).c_str());
    }
    // gather memory info about this host
    vtksys::SystemInformation sysInfo;
    std::string hostname = sysInfo.GetHostname();
    std::cout <<"Rank " << std::setw(3) <<  rank << " Hostname " << hostname.c_str() << " Initializing CudaGL with GPU " << displaynum << std::endl;
    //
    vtkPistonPolygonsPainter::InitCudaGL(rview->GetRenderWindow(), rank, displaynum);
  }
  return this->Superclass::AddToView(view);
}
//----------------------------------------------------------------------------
void vtkDepthSortRepresentation::SetUseDataPartitions(bool val)
{
  if (this->UseDataPartitions != val)
  {
    this->UseDataPartitions = val;
    this->MarkModified();
  }
}
//----------------------------------------------------------------------------
void vtkDepthSortRepresentation::SetOpacityArrayName(const char* opacity)
{
  this->TwoScalarsToColorsPainter->SetOpacityArrayName(opacity);
  this->PistonPolygonsPainter->SetOpacityArrayName(opacity);
  this->MarkModified();
}
//----------------------------------------------------------------------------
void vtkDepthSortRepresentation::SetEnableOpacity(int enable)
{
  this->TwoScalarsToColorsPainter->SetEnableOpacity(enable);
  this->DepthSortPainter->SetDepthSortRequired(enable);
  this->PistonPolygonsPainter->SetEnableOpacity(enable);
  this->MarkModified();
}
//----------------------------------------------------------------------------
int vtkDepthSortRepresentation::GetEnableOpacity()
{
  return this->TwoScalarsToColorsPainter->GetEnableOpacity();
}
//----------------------------------------------------------------------------
void vtkDepthSortRepresentation::SetEnablePiston(int mode)
{
  this->DepthSortDefaultPainter->SetEnablePiston(mode);
  if (mode) {
    this->PistonPolygonsPainter->SetScalarsToColors(this->TwoScalarsToColorsPainter);
    // swap the painter chain
  }
  this->MarkModified();
}
//----------------------------------------------------------------------------
void vtkDepthSortRepresentation::SetDepthSortEnableMode(int mode)
{
  this->DepthSortPainter->SetDepthSortEnableMode(mode);
  this->MarkModified();
}
//----------------------------------------------------------------------------
void vtkDepthSortRepresentation::SetDepthSortMode(int mode)
{
  this->DepthSortPainter->SetDepthSortMode(mode);
  this->MarkModified();
}
//----------------------------------------------------------------------------
void vtkDepthSortRepresentation::SetUseCachedSortOrder(int mode)
{
  this->DepthSortPainter->SetUseCachedSortOrder(mode);
  this->MarkModified();
}
//----------------------------------------------------------------------------
void vtkDepthSortRepresentation::SetDirection(int dir)
{
  if (dir>=0) {
    this->DepthSortPainter->SetDepthSortEnableMode(1); // if no depth peeling
    this->DepthSortPainter->SetDirection(dir);
  }
  else {
    this->DepthSortPainter->SetDepthSortEnableMode(2); // never
  }
  this->PistonPolygonsPainter->SetDirection(dir);
  this->MarkModified();
}
//----------------------------------------------------------------------------
int vtkDepthSortRepresentation::RequestData(vtkInformation* request,
                                            vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
  if (!this->Superclass::RequestData(request, inputVector, outputVector))
  {
    return 0;
  }

  if (inputVector[0]->GetNumberOfInformationObjects()==1 && this->UseDataPartitions) {
    // reduce bounds across processes in parallel.
    if (this->Controller->GetNumberOfProcesses() > 1)
    {
      double bounds_max[3] = {VTK_DOUBLE_MIN, VTK_DOUBLE_MIN};
      double bounds_min[3] = {VTK_DOUBLE_MAX, VTK_DOUBLE_MAX};
      if (vtkMath::AreBoundsInitialized(this->DataBounds))
      {
        bounds_min[0] = this->DataBounds[0];
        bounds_min[1] = this->DataBounds[2];
        bounds_min[2] = this->DataBounds[4];
        bounds_max[0] = this->DataBounds[1];
        bounds_max[1] = this->DataBounds[3];
        bounds_max[2] = this->DataBounds[5];
      }

      double reduced_bounds_max[3], reduced_bounds_min[3];
      this->Controller->AllReduce(bounds_max, reduced_bounds_max, 3, vtkCommunicator::MAX_OP);
      this->Controller->AllReduce(bounds_min, reduced_bounds_min, 3, vtkCommunicator::MIN_OP);

      this->GlobalDataBounds[0] = reduced_bounds_min[0];
      this->GlobalDataBounds[2] = reduced_bounds_min[1];
      this->GlobalDataBounds[4] = reduced_bounds_min[2];
      this->GlobalDataBounds[1] = reduced_bounds_max[0];
      this->GlobalDataBounds[3] = reduced_bounds_max[1];
      this->GlobalDataBounds[5] = reduced_bounds_max[2];
    }
  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkDepthSortRepresentation::ProcessViewRequest(
  vtkInformationRequestKey* request_type,
  vtkInformation* inInfo, vtkInformation* outInfo)
{
  if (!this->Superclass::ProcessViewRequest(request_type, inInfo, outInfo))
  {
    return 0;
  }

  if (request_type == vtkPVView::REQUEST_UPDATE())
    {

    if (this->GetEnableOpacity()) {
      outInfo->Set(vtkPVRenderView::NEED_ORDERED_COMPOSITING(), 1);
    }

    // override the setting made in the geometry representation, we do not allow any redistribution
    // because we assume that the data has been partitioned
    vtkPVRenderView::MarkAsRedistributable(inInfo, this, false);

    if (this->GetNumberOfInputConnections(0) == 1 && this->UseDataPartitions) {
      vtkAlgorithmOutput* connection = this->GetInputConnection(0, 0);
      vtkAlgorithm* inputAlgo = connection->GetProducer();
      vtkStreamingDemandDrivenPipeline* sddp = vtkStreamingDemandDrivenPipeline::SafeDownCast(inputAlgo->GetExecutive());
      vtkExtentTranslator* translator = sddp->GetExtentTranslator(connection->GetIndex());
      //sddp->GetWholeExtent(sddp->GetOutputInformation(connection->GetIndex()), whole_extent);

      //
      // Check the input to see if it has a bounds translator already initialized
      // with partition info 
      //
      this->BoundsTranslator = vtkBoundsExtentTranslator::SafeDownCast(translator);
      // if the extent translator has not been initialized well - don't use it
      if (this->BoundsTranslator && this->BoundsTranslator->GetNumberOfPieces()==0) {
        this->BoundsTranslator = NULL;
      }
      //
      // we should have enough information to generate a bounds translator if one was not supplied
      //
      if (!this->BoundsTranslator) {
/*
        this->BoundsTranslator = vtkSmartPointer<vtkBoundsExtentTranslator>::New();
        vtkBoundingBox box(this->DataBounds);
        // shrink the box by 25% of the shortest side.
        double lengths[3];
        box.GetLengths(lengths);
        double min_side = std::min(lengths[0], std::min(lengths[1], lengths[2]));
        box.Inflate(min_side*0.25*0.5);
        double shrunkenbounds[6];
        box.GetBounds(shrunkenbounds);
        this->BoundsTranslator->ExchangeBoundsForAllProcesses(this->Controller, shrunkenbounds);
        this->BoundsTranslator->InitWholeBounds();

        int whole_extent[6] = {0, 8191, 0, 8191, 0, 8191};
        this->BoundsTranslator->SetWholeExtent(whole_extent);
        //
        double *whole_bounds = this->BoundsTranslator->GetWholeBounds();
        // pick an arbitrary resolution of 4096^3 
        double spacing[3] = { 
          (whole_bounds[1]-whole_bounds[0])/(whole_extent[1]+1.0), 
          (whole_bounds[3]-whole_bounds[2])/(whole_extent[3]+1.0),  
          (whole_bounds[5]-whole_bounds[4])/(whole_extent[5]+1.0)
        };
        this->BoundsTranslator->SetSpacing(spacing);
        double origin[3] = { this->GlobalDataBounds[0], this->GlobalDataBounds[2], this->GlobalDataBounds[4] };
        vtkPVRenderView::SetOrderedCompositingInformation(inInfo, this, this->BoundsTranslator, whole_extent, origin, spacing);
*/
      }
      else {
        // 
        // The user supplied a translator of the right kind, lets use it
        //
        int whole_extent[6] = {1, -1, 1, -1, 1, -1};
        this->BoundsTranslator->GetWholeExtent(whole_extent);
        //
        double origin[3] = { this->GlobalDataBounds[0], this->GlobalDataBounds[2], this->GlobalDataBounds[4] };
        double spacing[3] = {
          (this->GlobalDataBounds[1] - this->GlobalDataBounds[0]) / (whole_extent[1] - whole_extent[0] + 1),
          (this->GlobalDataBounds[3] - this->GlobalDataBounds[2]) / (whole_extent[3] - whole_extent[2] + 1),
          (this->GlobalDataBounds[5] - this->GlobalDataBounds[4]) / (whole_extent[5] - whole_extent[4] + 1)
        };

        vtkPVRenderView::SetOrderedCompositingInformation(inInfo, this,
          translator, whole_extent, origin, spacing, this->BoundsTranslator->GetKdTree());
      }
    }
    else if (!this->UseDataPartitions) {
      double origin[3] = {0, 0, 0};
      double spacing[3] = {1, 1, 1};
      int whole_extent[6] = {1, -1, 1, -1, 1, -1};

      // Unset the ordered compositing info, so that vtkPVRenderView will
      // redistribute the unstructured grid as needed to volume render it.
      vtkPVRenderView::SetOrderedCompositingInformation(inInfo, this,
        NULL, whole_extent, origin, spacing);
    }
  }
  return 1;
}
//----------------------------------------------------------------------------
void vtkDepthSortRepresentation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

