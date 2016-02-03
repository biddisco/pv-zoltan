/*=========================================================================

  Module                  : vtkPartitionOutline.h

  Copyright (C) CSCS - Swiss National Supercomputing Centre.
  You may use modify and and distribute this code freely providing
  1) This copyright notice appears on all copies of source code
  2) An acknowledgment appears with any substantial usage of the code
  3) If this code is contributed to any other open source project, it
  must not be reformatted such that the indentation, bracketing or
  overall style is modified significantly.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/
//
// .NAME vtkZoltanV2PartitionFilter Efficiently distribute datasets in parallel
// .SECTION Description
// vtkZoltanV2PartitionFilter is the abstract base class used for parallel 
// load balancing/partitioning using the Zoltan library from Trilinos.
//
// .SECTION See Also
// vtkParticlePartitionFilter, vtkMeshPartitionFilter
//
#ifndef __vtkZoltanV2PartitionFilter_h
#define __vtkZoltanV2PartitionFilter_h
//
#include <vector>                // std used throughout
#include <map>                   // std used throughout
#include <string>                // std used throughout
//
#include "vtkDataSetAlgorithm.h" // superclass
#include "vtkBoundingBox.h"      // used as parameter
#include "vtkSmartPointer.h"     // for memory safety
//
#include "vtkZoltanBasePartitionFilter.h"
#include "zoltan.h"

#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>


// standard vtk classes
class  vtkMultiProcessController;
class  vtkPoints;
class  vtkIdTypeArray;
class  vtkIntArray;
class  vtkPointSet;
class  vtkDataSetAttributes;
class  vtkCellArray;
class  vtkTimerLog;
class  vtkFieldData;
class  vtkPKdTree;
// our special extent translator
class  vtkBoundsExtentTranslator;
class vtkInformationDataObjectMetaDataKey;
class vtkInformationIntegerRequestKey;
class vtkInformationIntegerKey;


//----------------------------------------------------------------------------
//#define JB_DEBUG__
//#define EXTRA_ZOLTAN_DEBUG 1
//----------------------------------------------------------------------------
#ifdef EXTRA_ZOLTAN_DEBUG
  #define INC_PACK_COUNT pack_count++;
  #define INC_UNPACK_COUNT unpack_count++;
  #define INC_SIZE_COUNT size_count++;
  #define CLEAR_ZOLTAN_DEBUG pack_count = 0; size_count = 0; unpack_count = 0;
#else
  #define INC_PACK_COUNT 
  #define INC_UNPACK_COUNT 
  #define INC_SIZE_COUNT 
  #define CLEAR_ZOLTAN_DEBUG 
#endif
//----------------------------------------------------------------------------
#if defined JB_DEBUG__ && !defined VTK_WRAPPING_CXX
#define OUTPUTTEXT(a) std::cout <<(a); std::cout.flush();

  #undef vtkDebugMacro
  #define vtkDebugMacro(a)  \
  { \
    if (this->UpdatePiece>=0) { \
      vtkOStreamWrapper::EndlType endl; \
      vtkOStreamWrapper::UseEndl(endl); \
      vtkOStrStreamWrapper vtkmsg; \
      vtkmsg << "P(" << this->UpdatePiece << "): " a << "\n"; \
      OUTPUTTEXT(vtkmsg.str()); \
      vtkmsg.rdbuf()->freeze(0); \
    } \
  }

  #undef  vtkErrorMacro
  #define vtkErrorMacro(a) vtkDebugMacro(a)  
#endif
//----------------------------------------------------------------------------
//
// GCC has trouble resolving some templated function pointers, 
// we explicitly declare the types and then cast them as args where needed
//    
typedef int  (*zsize_fn) (void *, int , int , ZOLTAN_ID_PTR , ZOLTAN_ID_PTR , int *);
typedef void (*zpack_fn) (void *, int , int , ZOLTAN_ID_PTR , ZOLTAN_ID_PTR , int , int , char *, int *);
typedef void (*zupack_fn)(void *, int , ZOLTAN_ID_PTR , int , char *, int *);
typedef void (*zprem_fn) (void *, int , int , int , ZOLTAN_ID_PTR , ZOLTAN_ID_PTR , int *, int *, int , ZOLTAN_ID_PTR , ZOLTAN_ID_PTR , int *, int *, int *);

// Zoltan 2 typedefs
typedef int localId_t;
#ifdef HAVE_ZOLTAN2_LONG_LONG_INT
typedef long long globalId_t;
#else
typedef int globalId_t;
#endif
//----------------------------------------------------------------------------
class vtkInformationDoubleKey;
class vtkInformationDoubleVectorKey;
class vtkInformationIntegerKey;
//----------------------------------------------------------------------------
class VTK_EXPORT vtkZoltanV2PartitionFilter : public vtkZoltanBasePartitionFilter
{
  public:
    static vtkZoltanV2PartitionFilter *New();
    vtkTypeMacro(vtkZoltanV2PartitionFilter,vtkZoltanBasePartitionFilter);

  protected:
     vtkZoltanV2PartitionFilter();
    ~vtkZoltanV2PartitionFilter();

    virtual void InitializeZoltanLoadBalance();
    virtual void ExecuteZoltanPartition(vtkPointSet *output, vtkPointSet *input);
    virtual void GetZoltanBoundingBoxes(vtkBoundingBox &globalBounds);

//BTX
    template<typename U>
    friend struct vtkZoltan2Helper;
//ETX
    Teuchos::ParameterList ZoltanParams;

  private:
    vtkZoltanV2PartitionFilter(const vtkZoltanV2PartitionFilter&);  // Not implemented.
    void operator=(const vtkZoltanV2PartitionFilter&);  // Not implemented.
};

#endif
