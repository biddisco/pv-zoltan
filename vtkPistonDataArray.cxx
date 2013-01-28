/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPistonIntArray.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkPistonDataArray.h"
//
#include "vtkObjectFactory.h"
#include "vtkPistonReference.h"

//----------------------------------------------------------------------------
namespace vtkpiston {
  //forward declarations of methods defined in the cuda implementation
  void CopyDataArrayToGPU(vtkIntArray   *da, vtkPistonIntArray *pa);
  void CopyDataArrayToGPU(vtkFloatArray *da, vtkPistonFloatArray *pa);
  //forward declarations of methods defined in the cuda implementation
  void CopyDataArrayFromGPU(vtkPistonIntArray *pa, vtkIntArray *da);
  void CopyArrayFromGPU(vtkPistonFloatArray *pa, vtkFloatArray *da);
}

//----------------------------------------------------------------------------

vtkStandardNewMacro(vtkPistonIntArray);
vtkStandardNewMacro(vtkPistonFloatArray);

//----------------------------------------------------------------------------
vtkPistonIntArray::vtkPistonIntArray()
{
  this->Reference = new vtkPistonReference;
  this->OwnReference = true;
}

//----------------------------------------------------------------------------
vtkPistonIntArray::~vtkPistonIntArray()
{
  if (this->OwnReference)
    {
    delete this->Reference;
    }
}

//----------------------------------------------------------------------------
void vtkPistonIntArray::DeepCopy(vtkPistonIntArray* src)
{
  if (vtkPistonIntArray* const pdo = vtkPistonIntArray::SafeDownCast(src))
    {
    if (this->OwnReference)
    {
      delete this->Reference;
    }
    this->Reference = new vtkPistonReference(pdo->Reference);
    this->OwnReference = true;
    this->Modified();
    }

  this->Superclass::DeepCopy(src);
}

//----------------------------------------------------------------------------
int vtkPistonIntArray::GetReferredType()
{
  return this->Reference->type;
}

//----------------------------------------------------------------------------
void * vtkPistonIntArray::GetReferredData()
{
  return this->Reference->data;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkPistonFloatArray::vtkPistonFloatArray()
{
  this->Reference = new vtkPistonReference;
  this->OwnReference = true;
}

//----------------------------------------------------------------------------
vtkPistonFloatArray::~vtkPistonFloatArray()
{
  if (this->OwnReference)
    {
    delete this->Reference;
    }
}

//----------------------------------------------------------------------------
void vtkPistonFloatArray::DeepCopy(vtkPistonFloatArray* src)
{
  if (vtkPistonFloatArray* const pdo = vtkPistonFloatArray::SafeDownCast(src))
    {
    if (this->OwnReference)
    {
      delete this->Reference;
    }
    this->Reference = new vtkPistonReference(pdo->Reference);
    this->OwnReference = true;
    this->Modified();
    }

  this->Superclass::DeepCopy(src);
}

//----------------------------------------------------------------------------
int vtkPistonFloatArray::GetReferredType()
{
  return this->Reference->type;
}

//----------------------------------------------------------------------------
void * vtkPistonFloatArray::GetReferredData()
{
  return this->Reference->data;
}

