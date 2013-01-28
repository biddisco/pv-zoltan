/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPistonIntArray.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPistonIntArray - A GPU resident data array.
//
// .SECTION Description
// vtkPistonIntArray is a basic data structure for storing datasets on
// GPU. This class provides the infrastructure for the VTK pipeline to
// work with the data as it does the rest of the vtkDataObjects.
// The GPU side structures are managed through the internal
// vtkPistonReference instance to keep the GPU/CPU code conceptually
// distinct.
//
// .SECTION See Also
// vtkPistonReference

#ifndef __vtkPistonIntArray_h
#define __vtkPistonIntArray_h

#include "vtkAcceleratorsPistonModule.h" // For export macro
#include "vtkIntArray.h"
#include "vtkFloatArray.h"

class vtkInformation;
class vtkInformationVector;
class vtkPistonReference;
class vtkTimeStamp;
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
class VTK_EXPORT vtkPistonIntArray : public vtkIntArray
{
public:
  static vtkPistonIntArray* New();
  vtkTypeMacro(vtkPistonIntArray, vtkIntArray);

  // Description:
  // A convenience handle to get type of what is stored in the reference.
  int GetReferredType();

  // Description:
  // A convenience handle to get whatever is actually stored in the reference.
  void * GetReferredData();

//BTX
  // Description:
  // GPU level representation and storage this manages.
  vtkPistonReference *GetReference() { return this->Reference; };
//ETX
  // Description:
  // deep copy the data from src into this object.
  virtual void DeepCopy(vtkPistonIntArray* src);

protected:
   vtkPistonIntArray();
  ~vtkPistonIntArray();

  vtkPistonReference *Reference;
  bool OwnReference;

private:
  vtkPistonIntArray(const vtkPistonIntArray&); // Not implemented
  void operator=(const vtkPistonIntArray&); // Not implemented
};
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
class VTK_EXPORT vtkPistonFloatArray : public vtkFloatArray
{
public:
  static vtkPistonFloatArray* New();
  vtkTypeMacro(vtkPistonFloatArray, vtkFloatArray);

  // Description:
  // A convenience handle to get type of what is stored in the reference.
  int GetReferredType();

  // Description:
  // A convenience handle to get whatever is actually stored in the reference.
  void * GetReferredData();

//BTX
  // Description:
  // GPU level representation and storage this manages.
  vtkPistonReference *GetReference() { return this->Reference; };
//ETX
  // Description:
  // deep copy the data from src into this object.
  virtual void DeepCopy(vtkPistonFloatArray* src);

protected:
   vtkPistonFloatArray();
  ~vtkPistonFloatArray();

  vtkPistonReference *Reference;
  bool OwnReference;

private:
  vtkPistonFloatArray(const vtkPistonFloatArray&); // Not implemented
  void operator=(const vtkPistonFloatArray&); // Not implemented
};
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#endif
