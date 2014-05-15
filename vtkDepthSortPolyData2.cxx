/*=========================================================================

Program:   Visualization Toolkit
Module:    vtkDepthSortPolyData2.cxx

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkDepthSortPolyData2.h"

#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkPolyData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"

#include <algorithm>
#include <functional>

//
// depth, cellId, pts_offset
//
#ifdef PV_ZOLTAN_USE_BOOST_TUPLE
  #include <boost/tuple/tuple.hpp>
  typedef boost::tuple<double, vtkIdType, vtkIdType> depthInfo;
  #define tuple_namespace boost
#else
  #include <tuple>
  typedef std::tuple<double, vtkIdType, vtkIdType> depthInfo;
  #define tuple_namespace std
#endif
typedef std::vector<depthInfo> depthList;

//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkDepthSortPolyData2);
//-----------------------------------------------------------------------------
vtkDepthSortPolyData2::vtkDepthSortPolyData2()
{
  this->Direction          = VTK_DIRECTION_BACK_TO_FRONT;
  this->DepthSortMode      = VTK_SORT_FIRST_POINT;
  this->FastPolygonMode    = 1;
  this->UseCachedSortOrder = 1;
  this->SortingList        = new depthList;
}
//-----------------------------------------------------------------------------
vtkDepthSortPolyData2::~vtkDepthSortPolyData2()
{
  delete reinterpret_cast<depthList*>(this->SortingList);
}
//-----------------------------------------------------------------------------
void FloatOrDoubleArrayPointer(vtkDataArray *dataarray, float *&F, double *&D) {
  if (dataarray && vtkFloatArray::SafeDownCast(dataarray)) {
    F = vtkFloatArray::SafeDownCast(dataarray)->GetPointer(0);
    D = NULL;
  }
  if (dataarray && vtkDoubleArray::SafeDownCast(dataarray)) {
    D = vtkDoubleArray::SafeDownCast(dataarray)->GetPointer(0);
    F = NULL;
  }
  //
  if (dataarray && !F && !D) {
    vtkGenericWarningMacro(<< dataarray->GetName() << "must be float or double");
  }
}
//-----------------------------------------------------------------------------
template<typename T> 
void CentreBoundsFromPtIds(vtkIdType *pts, vtkIdType npts, T *points, T result[3]) 
{
  T bounds[6];
  T *pt = &points[pts[0]*3];
  bounds[0] = bounds[1] = pt[0];
  bounds[2] = bounds[3] = pt[1];
  bounds[4] = bounds[5] = pt[2];
  for (vtkIdType i=1; i<npts; i++) {
    T *pt = &points[pts[0]*3];
    bounds[0] = std::min(pt[0],bounds[0]);
    bounds[1] = std::max(pt[0],bounds[1]);
    bounds[2] = std::min(pt[1],bounds[2]);
    bounds[3] = std::max(pt[1],bounds[3]);
    bounds[4] = std::min(pt[2],bounds[4]);
    bounds[5] = std::max(pt[2],bounds[5]);
  }
  result[0] = (bounds[0]+bounds[1])/2.0;
  result[1] = (bounds[2]+bounds[3])/2.0;
  result[2] = (bounds[4]+bounds[5])/2.0;
}
//-----------------------------------------------------------------------------
template<typename T>
void insertionSort(T *data, vtkIdType N)
{
  vtkIdType i;
  T key;
  for (vtkIdType j=1; j<N; j++) {
    key = data[j];
    for (i=j-1; (i>=0) && (data[i]<key); i--) {
      data[i+1] = data[i];
    }
    data[i+1] = key;
  }
}
//-----------------------------------------------------------------------------
int vtkDepthSortPolyData2::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  // if the user has not requested fast mode, default to the standard sort
  if (!this->FastPolygonMode) {
    return this->vtkDepthSortPolyData::RequestData(request, inputVector, outputVector);
  }

  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Compute the sort vector from camera position
  double vectorD[3], origin[3];
  if ( this->Camera == NULL) {
    vtkErrorMacro(<<"Need a camera to sort");
    return 0;
  }
  this->ComputeProjectionVector(vectorD, origin);
  float vectorF[3] = { (float)vectorD[0], (float)vectorD[1], (float)vectorD[2] };
  //
  vtkCellArray *polys = input->GetPolys();
  vtkPoints   *points = input->GetPoints();
  vtkIdType  numCells = polys->GetNumberOfCells();
  //
  // Points are always float or double, optimize access by avoiding copies
  // we copy the pointer and use the data directly.
  //
  float     *pointsF = NULL;
  double    *pointsD = NULL;
  if (numCells>0) {
    FloatOrDoubleArrayPointer(points->GetData(), pointsF, pointsD);
  }

  vtkIdType *cellArrayData = polys->GetPointer();

  bool usingCachedSortOrder = (this->UseCachedSortOrder 
    && polys->GetMTime()<this->LastSortTime 
    && points->GetMTime()<this->LastSortTime 
    && this->DepthOrder);

  // allocate space for depth values and sorted order
  depthList *ListToSort = static_cast<depthList*>(this->SortingList);
  ListToSort->resize(numCells);
  vtkIdType *pts;
  vtkIdType  npts;
  // traverse polygon list and compute depth
  polys->InitTraversal();
  for (vtkIdType index=0; index<numCells; index++) {
    if (!usingCachedSortOrder) {
      // get N and pointer to point Ids for this cell
      polys->GetNextCell(npts, pts);
      // set the cell ID in 1st tuple entry
      tuple_namespace::get<1>((*ListToSort)[index]) = index;
      // set the offset to the cell {N,ptIds} in 2nd tuple entry
      tuple_namespace::get<2>((*ListToSort)[index]) = static_cast<vtkIdType>(pts-cellArrayData-1);
    }
    else {
      // get N and pointer to point Ids for this cell from last cached iteration
      vtkIdType val = tuple_namespace::get<2>((*ListToSort)[index]);
      npts =  cellArrayData[val];
      pts  = &cellArrayData[val+1];
    }

    if (this->DepthSortMode == VTK_SORT_FIRST_POINT) {
      // set depth using float/double operation
      if (pointsF) {
        float *x = &pointsF[pts[0]*3];
        tuple_namespace::get<0>(ListToSort->operator[](index)) = vtkMath::Dot(x,vectorF);
      }
      else {
        double *x = &pointsD[pts[0]*3];
        tuple_namespace::get<0>(ListToSort->operator[](index)) = vtkMath::Dot(x,vectorD);
      }
    }
    else // if (this->DepthSortMode == VTK_SORT_BOUNDS_CENTER)
    {
      // compute depth and store it in 0th entry of tuple
      if (pointsF) {
        float x[3];
        CentreBoundsFromPtIds<float>(pts, npts, pointsF, x);
        tuple_namespace::get<0>(ListToSort->operator[](index)) = vtkMath::Dot(x,vectorF);
      }
      else {
        double x[3];
        CentreBoundsFromPtIds<double>(pts, npts, pointsD, x);
        tuple_namespace::get<0>(ListToSort->operator[](index)) = vtkMath::Dot(x,vectorD);
      }
    }
//    else // VTK_SORT_PARAMETRIC_CENTER )
//    {
//    }
  }
  this->LastSortTime.Modified();
  //
  timer->StopTimer();
  double buildtime = timer->GetElapsedTime();
  //
  vtkSmartPointer<vtkTimerLog> timer2 = vtkSmartPointer<vtkTimerLog>::New();
  timer2->StartTimer();
  // Sort the tuples, using std::sort (quicksort?) if not cached, shellsort if cached
  if (!usingCachedSortOrder) {
    if (this->Direction==VTK_DIRECTION_BACK_TO_FRONT) {
      std::sort(ListToSort->begin(), ListToSort->end(), std::greater<depthInfo>());
    }
    else {
      std::sort(ListToSort->begin(), ListToSort->end(), std::less<depthInfo>());
    }
  }
  else {
    //    insertionSort<depthInfo>(&ListToSort->operator[](0), ListToSort->size());
    //    shellsort<depthInfo>(&ListToSort->operator[](0), ListToSort->size());
    if (this->Direction==VTK_DIRECTION_BACK_TO_FRONT) {
      std::stable_sort(ListToSort->begin(), ListToSort->end(), std::greater<depthInfo>());
    }
    else {
      std::stable_sort(ListToSort->begin(), ListToSort->end(), std::less<depthInfo>());
    }
  }
  timer2->StopTimer();
  double sorttime = timer2->GetElapsedTime();

  vtkSmartPointer<vtkTimerLog> timer3 = vtkSmartPointer<vtkTimerLog>::New();
  timer3->StartTimer();
  //
  if (!this->DepthOrder) {
    this->DepthOrder = vtkSmartPointer<vtkIdTypeArray>::New();
    this->DepthOrder->SetName("DepthOrder");
    this->DepthOrder->SetNumberOfComponents(2);
  }
  this->DepthOrder->SetNumberOfTuples(numCells);
  vtkIdType *DepthData = this->DepthOrder->GetPointer(0);
  for (vtkIdType cellId=0; cellId<numCells; cellId++) {
    // tuple consists of "sorted cell Id", "Offset into cellArray list for {N,PtIds}"
    DepthData[cellId*2+0] = tuple_namespace::get<1>(ListToSort->operator[](cellId));
    DepthData[cellId*2+1] = tuple_namespace::get<2>(ListToSort->operator[](cellId));
  }

  // Points are left alone
  output->ShallowCopy(input);
  output->GetCellData()->AddArray(this->DepthOrder);

  timer3->StopTimer();
  std::cout << setprecision(6) << "BuildTime : << " <<  buildtime << " Cached " << usingCachedSortOrder << " SortTime " << sorttime << " Finalize = " << timer2->GetElapsedTime() << std::endl;

  return 1;
}
//-----------------------------------------------------------------------------
