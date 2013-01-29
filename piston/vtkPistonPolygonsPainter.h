/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPistonPolygonsPainter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPistonPolygonsPainter - this painter paints polygons.
// .SECTION Description
// This painter renders Polys in vtkPolyData. It can render the polys
// in any representation (VTK_POINTS, VTK_WIREFRAME, VTK_SURFACE).

#ifndef __vtkPistonPolygonsPainter_h
#define __vtkPistonPolygonsPainter_h

#include "vtkPolygonsPainter.h"
#include "vtkSmartPointer.h"
#include "vtkDataSetToPiston.h"
#include "vtkTwoScalarsToColorsPainter.h"

class vtkRenderWindow;
class vtkCamera;
class vtkActor;

#define VTK_DIRECTION_NO_SORT      -1
#define VTK_DIRECTION_BACK_TO_FRONT 0
#define VTK_DIRECTION_FRONT_TO_BACK 1

class VTK_EXPORT vtkPistonPolygonsPainter : public vtkPolygonsPainter
{
public:
  static vtkPistonPolygonsPainter* New();
  vtkTypeMacro(vtkPistonPolygonsPainter, vtkPolygonsPainter);

  void SetScalarsToColors(vtkTwoScalarsToColorsPainter *painter);

  // Description:
  // Manually call this before any cuda filters are created
  // to use direct GPU rendering.
  static void InitCudaGL(vtkRenderWindow *rw, int rank, int displayId);

  // Description:
  // Return true if using cuda interop feature otherwise false.
  inline static bool IsEnabledCudaGL()
    {
    return CudaGLInitted;
    }

  // Description:
  // Release any graphics resources that are being consumed by this mapper.
  // The parameter window could be used to determine which graphic
  // resources to release.
  virtual void ReleaseGraphicsResources(vtkWindow *) {};

  void RenderOnGPU(vtkCamera *cam, vtkActor *act);

  void ComputeProjectionVector(
    vtkCamera *cam, vtkActor *act, 
    double vector[3], double origin[3]);

  // Description:
  // Set/Get the name of the second array to blend with.
  vtkSetStringMacro(OpacityArrayName);
  vtkGetStringMacro(OpacityArrayName);

  // Description:
  // Enable/disble this painter
  vtkSetMacro(EnableOpacity, int);
  vtkGetMacro(EnableOpacity, int);

  // Description:
  // Swap sort order 
  // VTK_DIRECTION_NO_SORT      -1
  // VTK_DIRECTION_BACK_TO_FRONT 0
  // VTK_DIRECTION_FRONT_TO_BACK 1

  vtkSetMacro(Direction, int);
  vtkGetMacro(Direction, int);
  vtkBooleanMacro(Direction, int);

protected:
  vtkPistonPolygonsPainter();
  ~vtkPistonPolygonsPainter();

  void PrepareForRendering(vtkRenderer* renderer, vtkActor* actor);

  // Description:
  // The actual rendering happens here. This method is called only when
  // SupportedPrimitive is present in typeflags when Render() is invoked.
  virtual int RenderPrimitive(unsigned long flags, vtkDataArray* n,
    vtkUnsignedCharArray* c, vtkDataArray* t, vtkRenderer* ren);

  vtkSmartPointer<vtkDataSetToPiston>            DataSetToPiston;
  vtkSmartPointer<vtkTwoScalarsToColorsPainter>  ScalarsToColors;
  vtkCamera                                     *Camera;
  vtkActor                                      *Actor;
  char         *OpacityArrayName;
  int           EnableOpacity;
  int           Direction;

protected:
  // Description:
  // Allocates buffers that are shared between CUDA and GL
  void PrepareDirectRenderBuffers(int nPoints, int nCells);

  static bool CudaGLInitted;

  class InternalInfo;
  InternalInfo *Internal;

private:
  vtkPistonPolygonsPainter(const vtkPistonPolygonsPainter&); // Not implemented.
  void operator=(const vtkPistonPolygonsPainter&); // Not implemented.

};


#endif

