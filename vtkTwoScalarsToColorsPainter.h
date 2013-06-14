/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTwoScalarsToColorsPainter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// .NAME vtkTwoScalarsToColorsPainter
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
// .SECTION Description
// this painter blends two arrays in the color.
// It sends the two array through two different lookup tables, and blends the result
// two blending modes are possible : Multiply and Add.
// in the Multiply mode, the two colors are scaled to 0-1, multiplied, and rescaled to 0-255.
// in the Add mode, the two components are simply added and clamped to 0-255.

#ifndef __vtkTwoScalarsToColorsPainter_h__
#define __vtkTwoScalarsToColorsPainter_h__

#include "vtkOpenGLScalarsToColorsPainter.h"
#include <vector>

class VTK_EXPORT vtkTwoScalarsToColorsPainter : public vtkOpenGLScalarsToColorsPainter
{
public :
  static vtkTwoScalarsToColorsPainter* New();
  vtkTypeMacro(vtkTwoScalarsToColorsPainter, vtkOpenGLScalarsToColorsPainter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set/Get the name of the second array to blend with.
  vtkSetStringMacro(OpacityArrayName);
  vtkGetStringMacro(OpacityArrayName);

  // Description:
  // Enable/disble this painter
  vtkSetMacro(EnableOpacity, int);
  vtkGetMacro(EnableOpacity, int);

  // Description:
  // if ScalarVisibility is disabled, then this filter creates an array
  // of constant colors and blend it with the opacity array.
  // This parameter says if this should be done on point/cell or field data.
  vtkSetMacro(OpacityScalarMode, int);
  vtkGetMacro(OpacityScalarMode, int);

  const char *GetArrayName() { return this->ArrayName; }
  // Description:
  // For alpha blending, we sometime premultiply the colors
  // with alpha and change the alpha blending function.
  // This call returns whether we are premultiplying or using
  // the default blending function.
  // Currently this checks if the actor has a texture, if not
  // it returns true.
  // TODO: It is possible to make this decision
  // dependent on key passed down from a painter upstream
  // that makes a more informed decision for alpha blending
  // depending on extensions available, for example.
  virtual int GetPremultiplyColorsWithAlpha(vtkActor* actor);

  //BTX
  std::vector<float>* ComputeScalarsColorsf();
  void GetScalarRange(double range[2]);
  //ETX

  unsigned char *GetRGBPointer();

protected:
  vtkTwoScalarsToColorsPainter();
  virtual ~vtkTwoScalarsToColorsPainter();

  // Description:
  // Create texture coordinates for the output assuming a texture for the
  // lookuptable has already been created correctly.
  // this->LookupTable is the lookuptable used.
  void MapScalarsToTexture(vtkDataSet* output,
    vtkDataArray* scalars, vtkDataSet* input);

  // Description:
  // Called just before RenderInternal(). We build the Color array here.
  virtual void PrepareForRendering(vtkRenderer* renderer, vtkActor* actor);

  // Description:
  // Generates the colors, if needed.
  // If multiply_with_alpha is set, the colors are multiplied with
  // alpha.
  virtual void MapScalars(vtkDataSet* output,
    double alpha, int multiply_with_alpha,
    vtkDataSet* input);

  // Description:
  // Called before RenderInternal() if the Information has been changed
  // since the last time this method was called.
  virtual void ProcessInformation(vtkInformation*);

  // Description:
  // Returns if we can use texture maps for scalar coloring. Note this doesn't
  // say we "will" use scalar coloring. It says, if we do use scalar coloring,
  // we will use a 1D texture.
  // When rendering multiblock datasets, if any 2 blocks provide different
  // lookup tables for the scalars, then also we cannot use textures. This case
  // can be handled if required.
  int CanUseTextureMapForColoring(vtkDataObject* input);

protected :
  char          *OpacityArrayName;
  int            EnableOpacity;
  int            OpacityScalarMode;
  vtkTimeStamp   BlendTime;
  std::vector<unsigned char> RGBTable;

  //BTX
  std::vector<float> GPUScalarsColors;
  vtkTimeStamp       GPUColourTime;
  //ETX

private:
  vtkTwoScalarsToColorsPainter(const vtkTwoScalarsToColorsPainter&); // Not implemented.
  void operator=(const vtkTwoScalarsToColorsPainter&); // Not implemented.
};

#endif
