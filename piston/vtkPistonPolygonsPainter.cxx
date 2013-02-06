/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPistonPolygonsPainter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkPistonPolygonsPainter.h"

#include "vtkgl.h"
#include "vtkOpenGLExtensionManager.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkTransform.h"
#include "vtkPainterDeviceAdapter.h"
//
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
//
#include "vtkObjectFactory.h"
#include "vtkTimerLog.h"
#include "vtkCompositeDataSet.h"
#include "vtkCompositeDataIterator.h"
//
#include "vtkPistonDataObject.h"
//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkPistonPolygonsPainter);
//-----------------------------------------------------------------------------
struct  cudaGraphicsResource; //keeps vtkpiston namespace decl from claiming it
#define NUM_INTEROP_BUFFERS 4
bool vtkPistonPolygonsPainter::CudaGLInitted = false;
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
namespace vtkpiston {
  // Forward declarations of methods defined in the cuda implementation
  int  GetCudaDeviceCount();
  void CudaGLInit(int device);
  int  QueryNumVerts(vtkPistonDataObject *id);
  int  QueryVertsPer(vtkPistonDataObject *id);
  int  QueryNumCells(vtkPistonDataObject *id);
  void CudaRegisterBuffer(struct cudaGraphicsResource **vboResource, GLuint vboBuffer);
  void CudaTransferToGL(vtkPistonDataObject *id, unsigned long dataObjectMTimeCache, vtkTwoScalarsToColorsPainter *psc,
       struct cudaGraphicsResource **vboResources, double alpha, bool &hasNormals, bool &hasColors, bool &useindexbuffers);
  bool AlmostEqualRelativeAndAbs(float A, float B, float maxDiff, float maxRelDiff);
  //
  void DepthSortPolygons(vtkPistonDataObject *id, double *cameravec, int direction);
}
//-----------------------------------------------------------------------------
class vtkPistonPolygonsPainter::InternalInfo {
public:
  InternalInfo() {
    this->BufferSize = 0;
    this->CellCount = 0;
    this->DataObjectMTimeCache = 0;
    this->clearBuffers();
  }
  ~InternalInfo() {
    this->clearBuffers();
  }
  void clearBuffers() {
    if (this->BufferSize != 0 && this->vboBuffers[0] != -1) {
      vtkgl::DeleteBuffers(NUM_INTEROP_BUFFERS, this->vboBuffers);
    }
    for (int i=0; i<NUM_INTEROP_BUFFERS; i++) {
      vboBuffers[i] = -1;
      vboResources[i] = NULL;
    }
  }

  int     BufferSize;
  int     CellCount;
  GLuint  vboBuffers[NUM_INTEROP_BUFFERS];
  struct  cudaGraphicsResource* vboResources[NUM_INTEROP_BUFFERS];
  //
  unsigned long DataObjectMTimeCache;
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
vtkPistonPolygonsPainter::vtkPistonPolygonsPainter()
{
  this->SetSupportedPrimitive(vtkPainter::POLYS);
  this->DataSetToPiston  = vtkSmartPointer<vtkDataSetToPiston>::New();
  this->OpacityArrayName = NULL;
  this->EnableOpacity    = 0;
  this->Direction        = VTK_DIRECTION_BACK_TO_FRONT;
  //
  this->Internal = new vtkPistonPolygonsPainter::InternalInfo();
}
//-----------------------------------------------------------------------------
vtkPistonPolygonsPainter::~vtkPistonPolygonsPainter()
{
  this->PrepareDirectRenderBuffers(0, 0);
  delete this->Internal;
  delete [] this->OpacityArrayName;
}
//-----------------------------------------------------------------------------
void vtkPistonPolygonsPainter::SetScalarsToColors(vtkTwoScalarsToColorsPainter *painter)
{
  this->ScalarsToColors = painter;
}
//-----------------------------------------------------------------------------
int device_binding(int mpi_rank)
{
  int local_rank = mpi_rank;
  int dev_count, use_dev_count, my_dev_id;
  char *str;

  if ((str = getenv ("MV2_COMM_WORLD_LOCAL_RANK")) != NULL)
  {
    local_rank = atoi (str);
    printf ("MV2_COMM_WORLD_LOCAL_RANK %s\n", str);
  }

  if ((str = getenv ("MPISPAWN_LOCAL_NPROCS")) != NULL)
  {
    //num_local_procs = atoi (str);
    printf ("MPISPAWN_LOCAL_NPROCS %s\n", str);
  }

  dev_count = vtkpiston::GetCudaDeviceCount();
  if ((str = getenv ("NUM_GPU_DEVICES")) != NULL)
  {
    use_dev_count = atoi (str);
    printf ("NUM_GPU_DEVICES %s\n", str);
  }
  else
  {
    use_dev_count = dev_count;
  }

  my_dev_id = local_rank % use_dev_count;
  printf ("local rank = %d dev id = %d\n", local_rank, my_dev_id);
  return my_dev_id;
}
//-----------------------------------------------------------------------------
void vtkPistonPolygonsPainter::InitCudaGL(vtkRenderWindow *rw, int rank, int displayId)
{
  if (!vtkPistonPolygonsPainter::CudaGLInitted)
  {
    vtkOpenGLExtensionManager *em = vtkOpenGLExtensionManager::New();
    em->SetRenderWindow(rw);
    if (!em->LoadSupportedExtension("GL_VERSION_1_5"))
    {
      cerr << "WARNING: Can not use direct piston rendering, reverting to CPU rendering path." << endl;
      em->FastDelete();
      return;
    }
    em->FastDelete();
    if (displayId<0 || displayId>=vtkpiston::GetCudaDeviceCount()) {
      // try another method to get the device ID
      displayId = device_binding(rank);
    }
    vtkPistonPolygonsPainter::CudaGLInitted = true;
    vtkpiston::CudaGLInit(displayId);
  }
}
//-----------------------------------------------------------------------------
void vtkPistonPolygonsPainter::PrepareDirectRenderBuffers(int nPoints, int nCells)
{
  if (nPoints==this->Internal->BufferSize && nCells==this->Internal->CellCount) {
    return;
  }
  if (this->Internal->BufferSize != 0) {
    this->Internal->clearBuffers();
  }

  this->Internal->BufferSize = nPoints;
  this->Internal->CellCount  = nCells;
  if (this->Internal->BufferSize == 0) {
    return;
  }

  // Prep shared mem buffer between gl and cuda
  vtkgl::GenBuffers(NUM_INTEROP_BUFFERS, this->Internal->vboBuffers);

  // points 3*n float {x,y,z}
  vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER,
    this->Internal->vboBuffers[0]);
  vtkgl::BufferData(vtkgl::ARRAY_BUFFER,
    this->Internal->BufferSize*3*sizeof(float), 0,
    vtkgl::DYNAMIC_DRAW);

  // normals 3*n float {n0,n1,n2}
  vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER,
    this->Internal->vboBuffers[1]);
  vtkgl::BufferData(vtkgl::ARRAY_BUFFER,
    this->Internal->BufferSize*3*sizeof(float), 0,
    vtkgl::DYNAMIC_DRAW);

  // colors 3*n float {R,G,B} 
  vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER,
    this->Internal->vboBuffers[2]);
  vtkgl::BufferData(vtkgl::ARRAY_BUFFER,
    this->Internal->BufferSize*4*sizeof(float), 0,
    vtkgl::DYNAMIC_DRAW);

  // indexes 3*nCells int : triangles assumed {a,b,c} 
  vtkgl::BindBuffer(vtkgl::ELEMENT_ARRAY_BUFFER,
    this->Internal->vboBuffers[3]);
  vtkgl::BufferData(vtkgl::ELEMENT_ARRAY_BUFFER,
    this->Internal->CellCount*3*sizeof(int), 0,
    vtkgl::DYNAMIC_DRAW);

  vtkpiston::CudaRegisterBuffer(&this->Internal->vboResources[0],
    this->Internal->vboBuffers[0]);
  vtkpiston::CudaRegisterBuffer(&this->Internal->vboResources[1],
    this->Internal->vboBuffers[1]);
  vtkpiston::CudaRegisterBuffer(&this->Internal->vboResources[2],
    this->Internal->vboBuffers[2]);
  vtkpiston::CudaRegisterBuffer(&this->Internal->vboResources[3],
    this->Internal->vboBuffers[3]);
}
//-----------------------------------------------------------------------------
void vtkPistonPolygonsPainter::RenderOnGPU(vtkCamera *cam, vtkActor *act)
{
  vtkPistonDataObject *id = vtkPistonDataObject::SafeDownCast(this->DataSetToPiston->GetOutputDataObject(0));
  if (!id) {
    return;
  }

  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  int nPoints = vtkpiston::QueryNumVerts(id);
  int nCells = vtkpiston::QueryNumCells(id);
  this->PrepareDirectRenderBuffers(nPoints, nCells);

  // Transfer what is in tdo to buffer and render it directly on card
  bool hasNormals = false;
  bool hasColors = false;
  bool useindexbuffers = false;

  double cameravec[3], origin[3];
  this->ComputeProjectionVector(cam, act, cameravec, origin);
  if (this->Direction>=0) {
    vtkpiston::DepthSortPolygons(id, cameravec, this->Direction);
  }
  vtkpiston::CudaTransferToGL(id, this->Internal->DataObjectMTimeCache,
    this->ScalarsToColors,
    this->Internal->vboResources, 
    act->GetProperty()->GetOpacity(),
    hasNormals, hasColors, useindexbuffers);

  // Draw the result
  glEnableClientState(GL_VERTEX_ARRAY);
  vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER, this->Internal->vboBuffers[0]);
  glVertexPointer(3, GL_FLOAT, 0, 0);

  if (hasNormals)
  {
    glEnableClientState(GL_NORMAL_ARRAY);
    vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER, this->Internal->vboBuffers[1]);
    glNormalPointer(GL_FLOAT, 0, 0);
  }

  if (hasColors)
  {
//    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
//    glEnable(GL_BLEND);
//    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
//    glEnable(GL_COLOR_MATERIAL);
    /*
    // because we have used premultiple_with_alpha
    // otherwise we'd use glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);

    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    glDepthFunc( GL_LEQUAL );
    glEnable( GL_DEPTH_TEST );
*/
    glEnableClientState(GL_COLOR_ARRAY);
    vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER, this->Internal->vboBuffers[2]);
    glColorPointer(4, GL_FLOAT, 0, 0);
  }
  else {
    //    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    //    glEnable(GL_COLOR_MATERIAL);
    //    glDepthFunc( GL_LEQUAL );
    //    glEnable( GL_DEPTH_TEST );
    //    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    //    glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glEnable(GL_BLEND);
    //glEnableClientState(GL_COLOR_ARRAY);
    //vtkgl::BindBuffer(vtkgl::ARRAY_BUFFER, this->Internal->vboBuffers[2]);
    //glColorPointer(4, GL_FLOAT, 0, 0);
  }

  if (useindexbuffers) {
    //
    int vertsPer = vtkpiston::QueryVertsPer(id);
    vtkgl::BindBuffer(vtkgl::ELEMENT_ARRAY_BUFFER, this->Internal->vboBuffers[3]);
    switch (vertsPer) {
    case 4:
      glDrawElements(GL_QUADS, nCells*4, GL_UNSIGNED_INT, (GLvoid*)0);
      break;
    case 3:
      glDrawElements(GL_TRIANGLES, nCells*3, GL_UNSIGNED_INT, (GLvoid*)0);
      break;
    default:
      glDrawElements(GL_POINTS, nCells*1, GL_UNSIGNED_INT, (GLvoid*)0);
    }
  }
  else {
    int vertsPer = vtkpiston::QueryVertsPer(id);
    switch (vertsPer) {
    case 4:
      glDrawArrays(GL_QUADS, 0, nPoints);
      break;
    case 3:
      glDrawArrays(GL_TRIANGLES, 0, nPoints);
      break;
    default:
      glDrawArrays(GL_POINTS, 0, nPoints);
    }
  }

  glDisableClientState(GL_VERTEX_ARRAY);
  if (hasNormals) glDisableClientState(GL_NORMAL_ARRAY);
  if (hasColors) glDisableClientState(GL_COLOR_ARRAY);

  // Update object modified time
  this->Internal->DataObjectMTimeCache = id->GetMTime();

  timer->StopTimer();
  double rendertime = timer->GetElapsedTime();
  //  std::cout << setprecision(6) << "RenderTime : << " <<  rendertime << std::endl;
}
//-----------------------------------------------------------------------------
void vtkPistonPolygonsPainter::ComputeProjectionVector(
  vtkCamera *cam, vtkActor *act, 
  double vector[3], double origin[3])
{
  double *focalPoint = cam->GetFocalPoint();
  double   *position = cam->GetPosition();
  double focalPt[3], pos[3];

  vtkSmartPointer<vtkTransform> Transform = vtkSmartPointer<vtkTransform>::New();
  Transform->SetMatrix(act->GetMatrix());
  Transform->Inverse();

  Transform->TransformPoint(focalPoint,focalPt);
  Transform->TransformPoint(position,pos);

  for (int i=0; i<3; i++) {
    vector[i] = focalPt[i] - pos[i];
    origin[i] = pos[i];
  }
}
//-----------------------------------------------------------------------------
void vtkPistonPolygonsPainter::PrepareForRendering(vtkRenderer* renderer, vtkActor* actor)
{
  vtkDataObject* input = this->GetInput();
  if (!input)
  {
    vtkErrorMacro("No input present.");
    return;
  }

  //
  // does a shallow copy of input to output
  //
  this->Superclass::PrepareForRendering(renderer, actor);

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
        // this will break with real composite data
        this->DataSetToPiston->SetInputData(input);
        if (this->EnableOpacity) {
          this->DataSetToPiston->SetOpacityArrayName(this->OpacityArrayName);
        }
        else {
          this->DataSetToPiston->SetOpacityArrayName(NULL);
        }
        this->DataSetToPiston->SetScalarArrayName(this->ScalarsToColors->GetArrayName());
        this->DataSetToPiston->Update();
        break;
      }
    }

    iter->FastDelete();
  }
  else
  {
    this->DataSetToPiston->SetInputData(input);
    if (this->EnableOpacity) {
      this->DataSetToPiston->SetOpacityArrayName(this->OpacityArrayName);
    }
    else {
      this->DataSetToPiston->SetOpacityArrayName(NULL);
    }
    this->DataSetToPiston->SetScalarArrayName(this->ScalarsToColors->GetArrayName());
    this->DataSetToPiston->Update();
  }
  //
  this->Camera = renderer->GetActiveCamera();
  this->Actor  = actor;
}
//-----------------------------------------------------------------------------
int vtkPistonPolygonsPainter::RenderPrimitive(unsigned long idx, vtkDataArray* n,
    vtkUnsignedCharArray* c, vtkDataArray* t, vtkRenderer* ren)
{
  this->RenderOnGPU(this->Camera, this->Actor);
  return 1;
}

