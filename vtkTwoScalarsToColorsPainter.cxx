/*=========================================================================

Program:   Visualization Toolkit
Module:    vtkTwoScalarsToColorsPainter.cxx

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

#include "vtkTwoScalarsToColorsPainter.h"

#include "vtkObjectFactory.h"
#include "vtkDataObject.h"
#include "vtkDataSet.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkCompositeDataSet.h"
#include "vtkCompositeDataIterator.h"
#include "vtkLookupTable.h"
#include "vtkAbstractMapper.h"
#include "vtkFloatArray.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkImageData.h"
#include "vtkOpenGL.h"
#include "vtkMapper.h"
//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkTwoScalarsToColorsPainter)
//-----------------------------------------------------------------------------
vtkTwoScalarsToColorsPainter::vtkTwoScalarsToColorsPainter()
{
  this->OpacityArrayName  = NULL;
  this->EnableOpacity     = false;
  this->OpacityScalarMode = VTK_SCALAR_MODE_USE_POINT_FIELD_DATA;
  this->InterpolateScalarsBeforeMapping = 0;
}
//-----------------------------------------------------------------------------
vtkTwoScalarsToColorsPainter::~vtkTwoScalarsToColorsPainter()
{
  this->SetOpacityArrayName(NULL);
}
//-----------------------------------------------------------------------------
int vtkTwoScalarsToColorsPainter::GetPremultiplyColorsWithAlpha(vtkActor* actor)
{
  if (actor && (actor->GetTexture() || actor->GetProperty()->GetNumberOfTextures() > 0)) {
    return 0;
  }
  //  usually we premultiply the colours with the alpha to simplify blend operation
  if (this->EnableOpacity) return 1;
  // let base class decide
  return this->Superclass::GetPremultiplyColorsWithAlpha(actor);
}
//-----------------------------------------------------------------------------
void vtkTwoScalarsToColorsPainter::ProcessInformation(vtkInformation* info)
{
  this->Superclass::ProcessInformation(info);

  // we will not use textures for colours because we add independent opacity 
  // therefore textures would need to be 2D (and potentially very large too)
  if (this->EnableOpacity) {
    this->SetInterpolateScalarsBeforeMapping(false);
  }
}
//-----------------------------------------------------------------------------
void vtkTwoScalarsToColorsPainter::PrepareForRendering(
  vtkRenderer* renderer, vtkActor* actor)
{
  if (!this->EnableOpacity) {
    this->Superclass::PrepareForRendering(renderer, actor);
    return;
  }

  vtkDataObject* input = this->GetInput();
  if (!input)
  {
    vtkErrorMacro("No input present.");
    return;
  }

  // If the input has changed, the output should also reflect
  if (!this->OutputData ||
    !this->OutputData->IsA(input->GetClassName()) ||
    this->OutputUpdateTime < this->MTime ||
    this->OutputUpdateTime < this->GetInput()->GetMTime())
  {
    if (this->OutputData)
    {
      this->OutputData->FastDelete();
      this->OutputData = 0;
    }
    // Create a shallow-copied clone with no output scalars.
    this->OutputData = this->NewClone(input);
    this->OutputUpdateTime.Modified();
  }

  if (!this->ScalarVisibility && !this->EnableOpacity)
  {
    // Nothing to do here.
    this->ColorTextureMap = 0;
    return;
  }

  // Build the colors.
  // As per the vtkOpenGLPolyDataMapper's claim, this
  // it not a very expensive task, as the colors are cached
  // and hence we do this always.

  // Determine if we are going to use a texture for coloring or use vertex
  // colors. This need to be determine before we iterate over all the blocks in
  // the composite dataset to ensure that we employ the technique for all the
  // blocks.
  this->ScalarsLookupTable = 0;
  int useTexture = 0;
  // Remove texture map if present.
  this->ColorTextureMap = 0;
  this->UsingScalarColoring = 0;

  // Now if we have composite data, we need to MapScalars for all leaves.
  if (input->IsA("vtkCompositeDataSet"))
  {
    vtkCompositeDataSet* cdInput = vtkCompositeDataSet::SafeDownCast(input);
    vtkCompositeDataSet* cdOutput = vtkCompositeDataSet::SafeDownCast(this->OutputData);
    vtkCompositeDataIterator* iter = cdInput->NewIterator();
    for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
    {
      vtkDataSet* pdInput = vtkDataSet::SafeDownCast(iter->GetCurrentDataObject());
      vtkDataSet* pdOutput = vtkDataSet::SafeDownCast(cdOutput->GetDataSet(iter));
      if (pdInput && pdOutput) {
        this->MapScalars(pdOutput, actor->GetProperty()->GetOpacity(),
          this->GetPremultiplyColorsWithAlpha(actor), pdInput);
      }
    }

    iter->FastDelete();
  }
  else
  {
    this->MapScalars(vtkDataSet::SafeDownCast(this->OutputData),
      actor->GetProperty()->GetOpacity(),
      this->GetPremultiplyColorsWithAlpha(actor), vtkDataSet::SafeDownCast(
      input));
  }
  this->LastUsedAlpha = actor->GetProperty()->GetOpacity();
  this->GetLookupTable()->SetAlpha(this->LastUsedAlpha);
  this->LastUsedMultiplyWithAlpha = this->GetPremultiplyColorsWithAlpha(actor);
}
//-----------------------------------------------------------------------------
//
//
//
//-----------------------------------------------------------------------------
static inline void vtkMultiplyColorsWithOpacity(vtkDataArray* array, vtkDataArray *opacity, bool multiply_with_alpha)
{
  vtkUnsignedCharArray* colors = vtkUnsignedCharArray::SafeDownCast(array);
  if (!opacity || !colors || colors->GetNumberOfComponents() != 4)
  {
    return;
  }
  unsigned char* ptr = colors->GetPointer(0);

  vtkIdType numValues = colors->GetNumberOfTuples() * colors->GetNumberOfComponents();
  if (numValues < 4) {
    return;
  }

  vtkIdType tuple = 0;
  for (vtkIdType cc=0; cc<numValues; cc+=4, ptr+=4, tuple++) {
    // the alpha component set by actor opacity
    double alpha = (0x0ff & ptr[3])/255.0;

    // @TODO : assume values are 0-1 for now
    alpha = alpha * opacity->GetTuple1(tuple);
    if (multiply_with_alpha) {
      ptr[0] = static_cast<unsigned char>(0x0ff & static_cast<int>((0x0ff&ptr[0])*alpha + 0.5));
      ptr[1] = static_cast<unsigned char>(0x0ff & static_cast<int>((0x0ff&ptr[1])*alpha + 0.5));
      ptr[2] = static_cast<unsigned char>(0x0ff & static_cast<int>((0x0ff&ptr[2])*alpha + 0.5));
    }
    ptr[3]   = static_cast<unsigned char>(0x0ff & static_cast<int>(0x0ff*alpha + 0.5));
  }
}

//-----------------------------------------------------------------------------
// This method is copied directly from vtkScalarsToColors, with the addition that
// we blend in our per vertex opacity array
// we do not create any textures or use them anywhere.
void vtkTwoScalarsToColorsPainter::MapScalars(vtkDataSet* output, double alpha, 
  int multiply_with_alpha, vtkDataSet* input)
{
  if (!this->EnableOpacity) {
    this->Superclass::MapScalars(output, alpha, multiply_with_alpha, input);
    return;
  }

  int cellFlag;
  double orig_alpha;
  vtkDataArray* scalars = vtkAbstractMapper::GetScalars(input,
    this->ScalarMode, this->ArrayAccessMode, this->ArrayId,
    this->ArrayName, cellFlag);

  // We want the opacity array to be the same type as the color array (cell/point). 
  // if no colours ... then a blank LUT with just opacity per vertex/cell
  vtkDataArray* opacity;
  if (this->ScalarVisibility)
  {
    // if we map scalars to colors, then the opacity array has to
    // have the same scalarmode.
    opacity = vtkAbstractMapper::GetScalars(input, this->ScalarMode,
      VTK_GET_ARRAY_BY_NAME, -1, this->OpacityArrayName, cellFlag);

  }
  else
  { // no scalar color array, let us build one with constant color
    opacity = vtkAbstractMapper::GetScalars(input, this->OpacityScalarMode,
      VTK_GET_ARRAY_BY_NAME, -1, this->OpacityArrayName, cellFlag);

  }

  vtkPointData* oppd = output->GetPointData();
  vtkCellData*  opcd = output->GetCellData();
  vtkFieldData* opfd = output->GetFieldData();

  int arraycomponent = this->ArrayComponent;
  // This is for a legacy feature: selection of the array component to color by
  // from the mapper.  It is now in the lookuptable.  When this feature
  // is removed, we can remove this condition.
  if (scalars == 0 || scalars->GetNumberOfComponents() <= this->ArrayComponent)
  {
    arraycomponent = 0;
  }

  if (!this->ScalarVisibility || scalars == 0 || input == 0)
  {
    return;
  }

  // Let subclasses know that scalar coloring was employed in the current pass.
  // it is used in opengl scalars to colours as follows :
  // "if we are doing vertex colors then set lmcolor to adjust
  // the current materials ambient and diffuse values using
  // vertex color commands otherwise tell it not to".

  this->UsingScalarColoring = 1;

  vtkScalarsToColors* lut = 0;
  // Get the lookup table.
  if (scalars->GetLookupTable())
  {
    lut = scalars->GetLookupTable();
  }
  else
  {
    lut = this->GetLookupTable();
    lut->Build();
  }

  if (!this->UseLookupTableScalarRange)
  {
    lut->SetRange(this->ScalarRange);
  }

  // Try to reuse the old colors.
  vtkDataArray* colors;
  if (cellFlag == 0)
  {
    colors = oppd->GetScalars();
  }
  else if (cellFlag == 1)
  {
    colors = opcd->GetScalars();
  }
  else
  {
    colors = opfd->GetArray("Color");
  }

  // The LastUsedAlpha checks ensures that opacity changes are reflected
  // correctly when this->MapScalars(..) is called when iterating over a
  // composite dataset.
  if (colors && opacity && 
    this->LastUsedAlpha == alpha &&
    this->LastUsedMultiplyWithAlpha == multiply_with_alpha)
  {
    if (this->GetMTime() < colors->GetMTime() &&
      input->GetMTime() < colors->GetMTime() &&
      lut->GetMTime() < colors->GetMTime() &&
      opacity->GetMTime() < colors->GetMTime())
    {
      // using old colors.
      return;
    }
  }

  // Get rid of old colors.
  colors = 0;
  orig_alpha = lut->GetAlpha();
  lut->SetAlpha(alpha);
  colors = lut->MapScalars(scalars, this->ColorMode, arraycomponent);
  lut->SetAlpha(orig_alpha);

  // It is possible that the LUT simply returns the scalars as the
  // colors. In which case, we allocate a new array to ensure
  // that we don't modify the array in the input.
  if (scalars == colors)
  {
    // Since we will be changing the colors array
    // we create a copy.
    colors->FastDelete();
    colors = scalars->NewInstance();
    colors->DeepCopy(scalars);
    colors->SetName("Color");
  }
  else {
    colors->SetName("Color");
  }
  vtkMultiplyColorsWithOpacity(colors, opacity, multiply_with_alpha);

  if (cellFlag == 0)
  {
    oppd->SetScalars(colors);
  }
  else if (cellFlag == 1)
  {
    opcd->SetScalars(colors);
  }
  else
  {
    // Typically, when a name is assigned of the scalars array in PointData or CellData
    // it implies 3 component colors. This implication does not hold for FieldData.
    // For colors in field data, we use the component count of the color array
    // to decide if the colors are opaque colors.
    // These colors are nothing but cell colors,
    // except when rendering TStrips, in which case they represent
    // the triange colors.
    colors->SetName("Color");
    opfd->AddArray(colors);
  }
  colors->FastDelete();
}
//-----------------------------------------------------------------------------
void vtkTwoScalarsToColorsPainter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
//-----------------------------------------------------------------------------
unsigned char *vtkTwoScalarsToColorsPainter::GetRGBAPointer()
{
  if (!this->LookupTable) {
    vtkErrorMacro(<< "Invalid look up table");
    return NULL;
  }
  if (!this->UseLookupTableScalarRange) {
    this->LookupTable->SetRange(this->ScalarRange);
  }
  this->LookupTable->Build();
  if (vtkLookupTable::SafeDownCast(this->LookupTable)) {
    return vtkLookupTable::SafeDownCast(this->LookupTable)->GetPointer(0);
  }
  if (!(this->LookupTable->GetMTime() > this->GetMTime() || this->GPUColourTime.GetMTime() < this->GetMTime())) {
    return &this->RGBTable[0];
  }
  this->RGBTable.resize(256*4);
  
  //
  // generate 256 values of scalar ramp from range[0] to range[1]
  //
  double range[2];
  std::vector<float> values(256);
  this->GetScalarRange(range);
  for (int i=0; i<256; ++i) {
    values[i] = range[0] + i * ((range[1] - range[0]) / 256.0);
  }

  //
  // generate RGB values for the scalar ramp 
  //
  if (!this->UseLookupTableScalarRange) {
    this->LookupTable->SetRange(this->ScalarRange);
  }
  this->LookupTable->Build();
  this->LookupTable->MapScalarsThroughTable(&values[0], &this->RGBTable[0], VTK_FLOAT, 256, 1, 4);

  // touch build time (should be done last)
  this->GPUColourTime.Modified();
  return &this->RGBTable[0];
}
//-----------------------------------------------------------------------------
void vtkTwoScalarsToColorsPainter::GetScalarRange(double range[2])
{
  if (!this->UseLookupTableScalarRange) {
    range[0] = this->ScalarRange[0];
    range[1] = this->ScalarRange[1];
    return;
  }
  range[0] = this->LookupTable->GetRange()[0];
  range[1] = this->LookupTable->GetRange()[1];
}

