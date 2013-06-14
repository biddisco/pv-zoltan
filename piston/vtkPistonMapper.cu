#ifdef _WIN32
 #include <windows.h>
#endif

#include <iostream>

// include before gl_interop
#include "vtkgl.h"
//
// CUDA
//
#include <cuda_gl_interop.h>
#include <vector_types.h>
#include <vector_functions.h>

//
// Thrust
//
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>

//
// Piston
//
//#include "piston/piston_math.h"
//
#include "vtkPistonDataObject.h"
#include "vtkPistonDataWrangling.h"
#include "vtkPistonReference.h"

#include "../vtkTwoScalarsToColorsPainter.h"

#define USE_FLOAT_FOR_DEPTH_SORT

namespace vtkpiston {

inline __host__ __device__ float dot(float3 a, float3 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline __host__ __device__ float dot(double3 a, double3 b)
{ 
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
typedef thrust::device_vector<float>::iterator FloatIterator;
typedef thrust::tuple<FloatIterator, FloatIterator> FloatIteratorTuple;
typedef thrust::tuple<float&, float&> FloatTuple;
typedef thrust::zip_iterator<FloatIteratorTuple> Float2Iterator;
//typedef thrust::detail::tuple_of_iterator_references<float &, float &, thrust::null_type, thrust::null_type, thrust::null_type, thrust::null_type, thrust::null_type, thrust::null_type, thrust::null_type, thrust::null_type> TupleDef;

static __inline__ __host__ __device__ double3 make_double3(const float3 &f3)
{
  double3 t; t.x = f3.x; t.y = f3.y; t.z = f3.z; return t;
}

static __inline__ __host__ __device__ double3 make_double3(const float  &x, const float  &y, const float  &z)
{
  double3 t; t.x = x; t.y = y; t.z = z; return t;
}

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

struct color_map 
{
  const float min;
  const float max;
  const int   size;
  unsigned char *table;
  float       alpha;
  float      *opacity;
  float      *scalars;

  color_map(unsigned char *RGBtable, size_t arrSize, float rMin, float rMax, double a, float *scalararray, float *opacityarray) :
    min(rMin),
    max(rMax),
    size((arrSize / 3) - 1),
    table(RGBtable),
    alpha(a),
    scalars(scalararray),
    opacity(opacityarray)
    {
    }
//  const float inv = 0.0039215686274509803921568627451;
  // the internal calculation which is independent of template types
  __host__ __device__ inline uchar4 calc(float val, float opac) 
  { 
    // convert val to lookuptable index
    int index = 0;
    if ((max - min) > 0.0) {
      index = ( (val - min) / (max - min) ) * size;
    }
    if (index < 0)    index = 0;
    if (index > size) index = size;
    // convert to RGB tuple index
    // 1/255 = 0.0039215686274509803921568627451
    index *= 3; 
    return make_uchar4(table[index]*opac, table[index + 1]*opac, table[index + 2]*opac, 255.0*opac);
  };

  __host__ __device__ uchar4 color_map::operator()(float t)
  {
    float val = t;
    return calc(val, alpha);
  }

//  __host__ __device__ float4 color_map::operator()(const TupleDef &t)
//  {
//    float val  = thrust::get<0>(t);
//    float opac = thrust::get<1>(t)*alpha;
//    return calc(val, opac);
//  }

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
int GetCudaDeviceCount()
{
  int num;
  cudaGetDeviceCount(&num);
  return num;
}
//------------------------------------------------------------------------------
void CudaGLInit(int device)
{
//  cudaDeviceProp prop;
    // Fill it with zeros
//  memset(&prop,0,sizeof(cudaDeviceProp));
    // Pick a GPU capable of 1.0 or better
//  prop.major=1; prop.minor=0;
//  cudaChooseDevice(&dev,&prop);
 
  cudaSetDevice(device);

  // Set OpenGL device
  cudaError_t res = cudaGLSetGLDevice(device);

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
// Compute the distance from the camera to a single point
//------------------------------------------------------------------------------
template <typename T>
struct distance_functor 
{
  T cameravector;

  // construct with a constant camera vector
  __host__ __device__ distance_functor(T &cam) : cameravector(cam) {}

  template <typename Tuple>
  __host__ __device__
  void operator()(Tuple t)
  {
    thrust::get<1>(t) = dot(thrust::get<0>(t), cameravector);  
  }
};
//------------------------------------------------------------------------------
// Takes 3 point distances and computes the mean (i.e. triangular cell) distance
//------------------------------------------------------------------------------
template <typename T>
struct celldistance_functor 
{
  const T *vertex_distances;
  
  // construct with a precomputed distance vector for every vertex
  __host__ __device__ celldistance_functor(T *v) : vertex_distances(v) {}
  
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

#ifdef USE_FLOAT_FOR_DEPTH_SORT
  // prepare an array for the distances
  pD->distances.resize(pD->points->size());
  // initialize our functor which will compute distance and store in a vector
  float3 cam = make_float3(cameravec[0], cameravec[1], cameravec[2]);

  distance_functor<float3> distance(cam);
  // apply distance functor using input and output arrays using zip_iterator
  thrust::for_each(
    thrust::make_zip_iterator(thrust::make_tuple(pD->points->begin(), pD->distances.begin())),
    thrust::make_zip_iterator(thrust::make_tuple(pD->points->end(),   pD->distances.end())),
    distance);
  //
  // To compute the average distance for each cell, we must
  // sum/gather 3 distances (one for each vertex) for every cell by
  // looking up the vertex indices from the cell array tuples
  //
  // prepare an array for the distances
  thrust::device_vector<float> cell_distances(pD->nCells);

  celldistance_functor<float> celldist(thrust::raw_pointer_cast(pD->distances.data()));

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

#else
  // prepare an array for the distances
  thrust::device_vector<double> distances(pD->points->size());
  // initialize our functor which will compute distance and store in a vector
  double3 cam = make_double3(cameravec[0], cameravec[1], cameravec[2]);

  distance_functor<double3> distance(cam);
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
  thrust::device_vector<double> cell_distances(pD->nCells);

  celldistance_functor<double> celldist(thrust::raw_pointer_cast(distances.data()));

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
      thrust::greater_equal<double>());
  }
  else {
    thrust::sort_by_key(cell_distances.begin(), cell_distances.end(), pD->cells->begin(), 
      thrust::less<double>());
  }
#endif

}

//------------------------------------------------------------------------------
void CudaTransferToGL(vtkPistonDataObject *id, unsigned long dataObjectMTimeCache,
                      cudaGraphicsResource **vboResources,
                      unsigned char *colorptr,
                      double alpha,
                      bool &hasNormals, bool &hasColors, 
                      bool &useindexbuffers)
{
  vtkPistonReference *tr = id->GetReference();
  if (tr->type != VTK_POLY_DATA || tr->data == NULL) {
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
  uchar4 *colorsBufferData; 

  res = cudaGraphicsResourceGetMappedPointer((void **)&vertexBufferData, &num_bytes, vboResources[0]);
  if (res != cudaSuccess) {
    cerr << "Get mappedpointer for vertices failed ... "
         << cudaGetErrorString(res) << endl;
    return;
  }
  res = cudaGraphicsResourceGetMappedPointer((void **)&normalsBufferData, &num_bytes, vboResources[1]);
  if (res != cudaSuccess) {
    cerr << "Get mappedpointer for normals failed ... "
         << cudaGetErrorString(res) << endl;
    return;
  }
  res = cudaGraphicsResourceGetMappedPointer((void **)&colorsBufferData, &num_bytes, vboResources[2]);
  if (res != cudaSuccess)
  {
    cerr << "Get mappedpointer for colors failed ... "
         << cudaGetErrorString(res) << endl;
    return;
  }

  res = cudaGraphicsResourceGetMappedPointer((void **)&cellsBufferData, &num_bytes, vboResources[3]);
  if (res != cudaSuccess)
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


  if (pD->colors)
  {
    thrust::copy(pD->colors->begin(), pD->colors->end(), 
      thrust::device_ptr<uchar4>(colorsBufferData));
    hasColors = true;
  }
  else if (pD->distances.size()>0) {
    pD->colors->resize(pD->distances.size());
    float *raw_scalar = thrust::raw_pointer_cast(&pD->distances[0]);
    // Copy RGB values to GPU
    thrust::device_vector<unsigned char> onGPU(colorptr, colorptr+(256*4));
    unsigned char *raw_RGBA = thrust::raw_pointer_cast(&onGPU[0]);
    // create a lookuptable
    color_map colorMap(raw_RGBA, 256, 0, 2, 0.5, raw_scalar, NULL);
    // map all scalars through LUT
    thrust::transform(pD->distances.begin(), pD->distances.end(), pD->colors->begin(), colorMap);
  }
  else if (pD->scalars)
  {
/*
    std::vector<float> *colors = psc->ComputeScalarsColorsf();

    // Copy to GPU
    thrust::device_vector<float> onGPU(colors->begin(), colors->end());
    float *raw_ptr = thrust::raw_pointer_cast(&onGPU[0]);

    // Now run each scalar data through the map to choose a color for it

    double scalarRange[2];
    psc->GetScalarRange(scalarRange);

    float * opacitydata = pD->opacities ? thrust::raw_pointer_cast(pD->opacities->data()) : NULL;


  if (opacitydata==0) {
    color_map colorMap(raw_RGBA, rgba.size(), scalarRange[0], scalarRange[1], alpha, raw_scalar, raw_opacity);

    thrust::transform(scalardata.begin(), scalardata.end(), colours.begin(), colorMap);
  }
  else {
    // Now we'll create some zip_iterators for A and B
    auto first = thrust::make_zip_iterator(thrust::make_tuple(scalardata.begin(), opacitydata.begin()));
    auto last  = thrust::make_zip_iterator(thrust::make_tuple(scalardata.end(),   opacitydata.end()));

    color_map colorMap(raw_RGBA, rgba.size(), scalarRange[0], scalarRange[1], alpha, raw_scalar, raw_opacity);

    thrust::transform(first, last, colours.begin(), colorMap);
  }

    if (opacitydata) {
      // Now we'll create some zip_iterators for A and B
      Float2Iterator _first = thrust::make_zip_iterator(thrust::make_tuple(pD->scalars->begin(), pD->opacities->begin()));; 
      Float2Iterator  _last = thrust::make_zip_iterator(thrust::make_tuple(pD->scalars->end(),   pD->opacities->end()));

      color_map<FloatTuple> colorMap(raw_ptr, onGPU.size(), scalarRange[0], scalarRange[1], alpha, opacitydata);
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
*/
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
