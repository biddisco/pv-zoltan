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
// .NAME vtkDepthSortRepresentation
// .SECTION Description
// vtkDepthSortRepresentation is an extension for vtkGeometryRepresentation
// that renders point-sprites at all point locations.

#ifndef __vtkDepthSortRepresentation_h
#define __vtkDepthSortRepresentation_h

#include "vtkGeometryRepresentation.h"
#include "vtkSmartPointer.h"

class vtkDepthSortDefaultPainter;
class vtkDepthSortPainter;
class vtkDepthSortPolygonsPainter;
class vtkTwoScalarsToColorsPainter;
class vtkMultiProcessController;
class vtkBoundsExtentTranslator;
class vtkPistonPolygonsPainter;

class VTK_EXPORT vtkDepthSortRepresentation : public vtkGeometryRepresentation
{
public:
  static vtkDepthSortRepresentation* New();
  vtkTypeMacro(vtkDepthSortRepresentation, vtkGeometryRepresentation);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // When on (default), the representation tells the view to use the
  // partitioning information from the input structured grid for ordered
  // compositing. When off we let the view build its own ordering and
  // redistribute data as needed.
  void SetUseDataPartitions(bool);
  vtkGetMacro(UseDataPartitions, bool);

  // Description:
  // By default this filter uses the global controller,
  // but this method can be used to set another instead.
  virtual void SetController(vtkMultiProcessController*);

  // Description:
  // Set/Get the name of the second array to blend with.
  void SetOpacityArrayName(const char* opacity);

  void SetEnableOpacity(int enable);
  int  GetEnableOpacity();

  void SetDepthSortEnableMode(int mode);
  void SetDepthSortMode(int mode);
  void SetUseCachedSortOrder(int mode);
  void SetDirection(int dir);
  void SetEnablePiston(int mode);

  // Description:
  // Adds the representation to the view.  This is called from
  // vtkView::AddRepresentation().  Subclasses should override this method.
  // Returns true if the addition succeeds.
  virtual bool AddToView(vtkView* view);

protected:
  vtkDepthSortRepresentation();
  ~vtkDepthSortRepresentation();

  virtual void ReportReferences(vtkGarbageCollector *collector);

  // Description:
  virtual int RequestData(vtkInformation*,
    vtkInformationVector**, vtkInformationVector*);

  virtual int ProcessViewRequest(
    vtkInformationRequestKey* request_type,
    vtkInformation* inInfo, vtkInformation* outInfo);

  //
  vtkDepthSortDefaultPainter                   *DepthSortDefaultPainter;
  vtkSmartPointer<vtkDepthSortPainter>          DepthSortPainter;
  vtkSmartPointer<vtkDepthSortPolygonsPainter>  DepthSortPolygonsPainter;
  vtkSmartPointer<vtkTwoScalarsToColorsPainter> TwoScalarsToColorsPainter;
  vtkSmartPointer<vtkBoundsExtentTranslator>    BoundsTranslator;
#ifdef PV_ZOLTAN_USE_PISTON
  vtkSmartPointer<vtkPistonPolygonsPainter>     PistonPolygonsPainter;
#endif
 
  int UseDataPartitions;
  //
  vtkMultiProcessController *Controller;
  //
  double GlobalDataBounds[6];

private:
  vtkDepthSortRepresentation(const vtkDepthSortRepresentation&); // Not implemented
  void operator=(const vtkDepthSortRepresentation&); // Not implemented

};

#endif
