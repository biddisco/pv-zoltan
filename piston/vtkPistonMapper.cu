#ifdef _WIN32
 #include <windows.h>
#endif

#include <iostream>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <thrust/functional.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>

#include "vtkgl.h"
#include <cuda_gl_interop.h>

#include "piston/piston_math.h"
//
#include "vtkPistonDataObject.h"
#include "vtkPistonDataWrangling.h"
#include "vtkPistonReference.h"

#include "../vtkTwoScalarsToColorsPainter.h"

namespace vtkpiston {

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
typedef thrust::device_vector<float>::iterator FloatIterator;
typedef thrust::tuple<FloatIterator, FloatIterator> FloatIteratorTuple;
typedef thrust::tuple<float&, float&> FloatTuple;
typedef thrust::zip_iterator<FloatIteratorTuple> Float2Iterator;
//------------------------------------------------------------------------------
// The colour map struct is templated (T) over the iterator type that we will use
// in our transform function.
//
// Because we use complicated zip/tuple iterators which are partially dereferenced
// by the tuple/iterator code inside the transform call, we separately template
// the operator() over the final dereferenced (element) type that will be passed
// to the operator().
//
// Two specializations are provided for color tables with/without vertex opacity
// 1) scalars        : T=float(iterator), element=float 
// 2) scalar/opacity : T=FloatIteratorTuple(iterator), element=FloatTuple 

template<typename T>
struct color_map : public thrust::unary_function<T, float4>
{
  const float min;
  const float max;
  const int   size;
  float      *table;
  float       alpha;
  float      *opacity;

  color_map(float *table, int arrSize, float rMin, float rMax, double a, float *opacityarray) :
    min(rMin),
    max(rMax),
    size((arrSize / 3) - 1),
    table(table),
    alpha(a),
    opacity(opacityarray)
    {
    }
  
  // the internal calculation which is independent of template types
  __host__ __device__ inline float4 calc(float val, float opac) 
  { 
    // convert val to lookuptable index
    int index = 0;
    if ((max - min) > 0.0) {
      index = ( (val - min) / (max - min) ) * size;
    }
    if (index < 0)    index = 0;
    if (index > size) index = size;
    // convert to RGB tuple index
    index *= 3; 
    return make_float4(table[index], table[index + 1], table[index + 2], opac);
  };

  // declare an empty general templated functor operator which we will specialize later
  template<typename element>
  __host__ __device__ float4 operator()(element t) { 
    // should throw ("Error - this function must be specialized for the type used");
    return make_float4(1,1,1,1);
  }
};
//------------------------------------------------------------------------------
// doubly templated specialization
// overload the colormap operator() for tuple iterators which will hold <scalar, opacity>
template <> template<>
__host__ __device__ float4 color_map<FloatIteratorTuple>::operator()<FloatTuple>(FloatTuple t)
{
  float val  = thrust::get<0>(t);
  float opac = thrust::get<1>(t)*alpha;
  return calc(val, opac);
}
//------------------------------------------------------------------------------
// doubly templated specialization
// overload the colormap operator() for single color array
template <> template<>
__host__ __device__ float4 color_map<float>::operator()<float>(float t)
{
  float val = t;
  return calc(val, alpha);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void CudaGLInit()
{
  cudaDeviceProp prop;
  int dev;

  // Fill it with zeros
  memset(&prop,0,sizeof(cudaDeviceProp));

  // Pick a GPU capable of 1.0 or better
  prop.major=1; prop.minor=0;
  cudaChooseDevice(&dev,&prop);

  // Set OpenGL device
  cudaError_t res = cudaGLSetGLDevice(dev);

  if (res != cudaSuccess)
    {
    std::cerr << "Set device failed  ... " << cudaGetErrorString(res) << endl;
    return;
    }
}

//------------------------------------------------------------------------------
void CudaRegisterBuffer(struct cudaGraphicsResource **vboResource,
                        GLuint vboBuffer)
{
  cudaError_t res =
    cudaGraphicsGLRegisterBuffer(vboResource, vboBuffer,
                                cudaGraphicsMapFlagsWriteDiscard);
  if (res != cudaSuccess)
  {
    std::cerr << "Register buffer failed ... " << cudaGetErrorString(res) << endl;
    return;
  }
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
struct distance_functor 
{
  float3 cameravector;

  // construct with a constant camera vector
  __host__ __device__ distance_functor(float3 &cam) : cameravector(cam) {}

  template <typename Tuple>
  __host__ __device__
  void operator()(Tuple t)
  {
    thrust::get<1>(t) = dot(thrust::get<0>(t), cameravector);  
  }
};
//------------------------------------------------------------------------------
struct celldistance_functor 
{
  const float *vertex_distances;
  
  // construct with a precomputed distance vector for every vertex
  __host__ __device__ celldistance_functor(float *v) : vertex_distances(v) {}
  
  template <typename Tuple>
  __host__ __device__
  void operator()(Tuple t)
  { 
    thrust::get<1>(t) = (vertex_distances[thrust::get<0>(t).x] + 
                         vertex_distances[thrust::get<0>(t).y] +
                         vertex_distances[thrust::get<0>(t).z])/3.0;
  }
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void DepthSortPolygons(vtkPistonDataObject *id, double *cameravec, int direction)
{
  vtkPistonReference *tr = id->GetReference();
  if (tr->type != VTK_POLY_DATA || tr->data == NULL) {
    // Type mismatch, don't bother trying
    return;
  }
  vtk_polydata *pD = (vtk_polydata *)tr->data;

  //
  // we need to compute the distance to the camera for each cell.
  // Perform a dot product of each vertex with the supplied camera vector
  //

  // prepare an array for the distances
  thrust::device_vector<float> distances(pD->points->size());

  // initialize our functor which will compute distance and store in a vector
  float3 cam = make_float3(cameravec[0], cameravec[1], cameravec[2]);
  distance_functor distance(cam);

  // apply distance functor using input and output arrays using zip_iterator
  thrust::for_each(
    thrust::make_zip_iterator(thrust::make_tuple(pD->points->begin(), distances.begin())),
    thrust::make_zip_iterator(thrust::make_tuple(pD->points->end(),   distances.end())),
    distance);

  //
  // To compute the average distance for each cell, we must
  // sum/gather 3 distances (one for each vertex) for every cell by
  // looking up the vertex indices from the cell array tuples
  //

  // prepare an array for the distances
  thrust::device_vector<float> cell_distances(pD->nCells);

  celldistance_functor celldist(thrust::raw_pointer_cast(distances.data()));

  thrust::for_each(
    thrust::make_zip_iterator(thrust::make_tuple(pD->cells->begin(), cell_distances.begin())),
    thrust::make_zip_iterator(thrust::make_tuple(pD->cells->end(),   cell_distances.end())),
    celldist);

  //
  // now we want to sort the cells using the average distance
  // we must copy the cell vertex index tuple during the sort
  //
  if (direction==0) {
    thrust::sort_by_key(cell_distances.begin(), cell_distances.end(), pD->cells->begin(), 
      thrust::greater<float>());
  }
  else {
    thrust::sort_by_key(cell_distances.begin(), cell_distances.end(), pD->cells->begin(), 
      thrust::less<float>());
  }
}

//------------------------------------------------------------------------------
void CudaTransferToGL(vtkPistonDataObject *id, unsigned long dataObjectMTimeCache,
                      vtkTwoScalarsToColorsPainter *psc,
                      cudaGraphicsResource **vboResources,
                      double alpha,
                      bool &hasNormals, bool &hasColors, 
                      bool &useindexbuffers)
{
  vtkPistonReference *tr = id->GetReference();
  if (tr->type != VTK_POLY_DATA || tr->data == NULL)
    {
    // Type mismatch, don't bother trying
    return;
    }

  vtk_polydata *pD = (vtk_polydata *)tr->data;

  // Claim access to buffer for cuda
  cudaError_t res;
  res = cudaGraphicsMapResources(4, vboResources, 0);
  if (res != cudaSuccess)
  {
    cerr << "Claim for CUDA failed ... " << cudaGetErrorString(res) << endl;
    return;
  }

  size_t num_bytes;
  float3 *vertexBufferData;
  uint3  *cellsBufferData;
  float  *normalsBufferData;
  float4 *colorsBufferData; 

  res = cudaGraphicsResourceGetMappedPointer
      ((void **)&vertexBufferData, &num_bytes, vboResources[0]);
  if(res != cudaSuccess) {
    cerr << "Get mappedpointer for vertices failed ... "
         << cudaGetErrorString(res) << endl;
    return;
  }
  res = cudaGraphicsResourceGetMappedPointer
      ((void **)&normalsBufferData, &num_bytes, vboResources[1]);
  if(res != cudaSuccess) {
    cerr << "Get mappedpointer for normals failed ... "
         << cudaGetErrorString(res) << endl;
    return;
  }
  res = cudaGraphicsResourceGetMappedPointer
      ((void **)&colorsBufferData, &num_bytes, vboResources[2]);
  if(res != cudaSuccess)
  {
    cerr << "Get mappedpointer for colors failed ... "
         << cudaGetErrorString(res) << endl;
    return;
  }

  res = cudaGraphicsResourceGetMappedPointer
      ((void **)&cellsBufferData, &num_bytes, vboResources[3]);
  if(res != cudaSuccess)
  {
    std::string errormsg = cudaGetErrorString(res);
    cerr << "Get mappedpointer for cell indices failed ... "
         << cudaGetErrorString(res) << endl;
    return;
  }

  // Copy on card verts to the shared on card gl buffer
  thrust::copy(pD->points->begin(), pD->points->end(),
               thrust::device_ptr<float3>(vertexBufferData));

  // Copy on card cell indices to the shared on card gl buffer
  if (pD->cells) {
    useindexbuffers = true;
    thrust::copy(pD->cells->begin(), pD->cells->end(),
                 thrust::device_ptr<uint3>(cellsBufferData));
  }

  hasNormals = false;
  if (pD->normals)
    {
    hasNormals = true;

    // Copy on card verts to the shared on card gl buffer
    thrust::copy(pD->normals->begin(), pD->normals->end(),
                 thrust::device_ptr<float>(normalsBufferData));
    }
  hasColors = false;


  if (0 && pD->colors)
  {
    thrust::copy(pD->colors->begin(), pD->colors->end(), 
      thrust::device_ptr<float4>(colorsBufferData));
  }
  else if (pD->scalars)
  {
/*
    double scalarRange[2];
    id->GetScalarsRange(scalarRange);

    hasColors = true;

    if(id->GetMTime() > dataObjectMTimeCache)
      {
      vtkPiston::minmax_pair<float> result = vtkPiston::find_min_max(
                                              pD->scalars);

      scalarRange[0] = static_cast<double>(result.min_val);
      scalarRange[1] = static_cast<double>(result.max_val);

      // Set parameters to compute scalars colors
      const int numvalues = 256;
      id->SetScalarsRange(scalarRange);
      psc->SetTableRange(scalarRange[0], scalarRange[1]);
      psc->SetNumberOfValues(numvalues);
      }
*/

    std::vector<float> *colors = psc->ComputeScalarsColorsf();

    // Copy to GPU
    thrust::device_vector<float> onGPU(colors->begin(), colors->end());
    float *raw_ptr = thrust::raw_pointer_cast(&onGPU[0]);

    // Now run each scalar data through the map to choose a color for it

    double scalarRange[2];
    psc->GetScalarRange(scalarRange);

    float * opacitydata = pD->opacities ? thrust::raw_pointer_cast(pD->opacities->data()) : NULL;

    if (opacitydata) {
      // Now we'll create some zip_iterators for A and B
      Float2Iterator _first = thrust::make_zip_iterator(thrust::make_tuple(pD->scalars->begin(), pD->opacities->begin()));; 
      Float2Iterator  _last = thrust::make_zip_iterator(thrust::make_tuple(pD->scalars->end(),   pD->opacities->end()));

      color_map<FloatIteratorTuple> colorMap(raw_ptr, onGPU.size(), scalarRange[0], scalarRange[1], alpha, opacitydata);
      thrust::copy(thrust::make_transform_iterator(_first, colorMap),
                   thrust::make_transform_iterator(_last,  colorMap),
                   thrust::device_ptr<float4>(colorsBufferData));
    }
    else {
      color_map<float> colorMap(raw_ptr, onGPU.size(), scalarRange[0], scalarRange[1], alpha, opacitydata);
      thrust::copy(thrust::make_transform_iterator(pD->scalars->begin(), colorMap),
                   thrust::make_transform_iterator(pD->scalars->end(),   colorMap),
                   thrust::device_ptr<float4>(colorsBufferData));
    }
  }

  // Allow GL to access again
  res = cudaGraphicsUnmapResources(4, vboResources, 0);
  if (res != cudaSuccess)
  {
    cerr << "Release from CUDA failed ... " << cudaGetErrorString(res) << endl;
    return;
  }

  return;
}
//------------------------------------------------------------------------------
} //namespace
