/*=========================================================================

  Module                  : vtkZoltanV1PartitionFilter.h
  Revision of last commit : $Rev: 884 $
  Author of last commit   : $Author: biddisco $
  Date of last commit     : $Date:: 2010-04-06 12:03:55 +0200 #$

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
// .NAME vtkZoltanV1PartitionFilter distribute datasets in parallel
// .SECTION Description
// vtkZoltanV1PartitionFilter is a parallel load balancing/partitioning 
// filter for datasets. It uses the Zoltan library from the Trilinos 
// package to perform the redistribution.

#ifndef __vtkZoltanV1PartitionFilter_h
#define __vtkZoltanV1PartitionFilter_h
//
#include <vector>
#include <map>
//
#include "vtkDataSetAlgorithm.h" // superclass
#include "vtkBoundingBox.h"
#include "vtkSmartPointer.h"
//
#include "zoltan.h"
// standard vtk classes
class  vtkMultiProcessController;
class  vtkPoints;
class  vtkIdTypeArray;
class  vtkIntArray;
class  vtkPointSet;
class  vtkDataSetAttributes;
class  vtkCellArray;
class  vtkTimerLog;
// our special extent translator
class  vtkBoundsExtentTranslator;

//----------------------------------------------------------------------------
// #define JB_DEBUG__
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
#define JB_DEBUG__

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
//----------------------------------------------------------------------------

class VTK_EXPORT vtkZoltanV1PartitionFilter : public vtkDataSetAlgorithm
{
  public:
    static vtkZoltanV1PartitionFilter *New();
    vtkTypeMacro(vtkZoltanV1PartitionFilter,vtkDataSetAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    // Description:
    // By default this filter uses the global controller,
    // but this method can be used to set another instead.
    virtual void SetController(vtkMultiProcessController*);

    // Description:
    // Specify the name of a scalar array which will be used to fetch
    // the index of each point. This is necessary only if the particles
    // change position (Id order) on each time step. The Id can be used
    // to identify particles at each step and hence track them properly.
    // If this array is NULL, the global point ids are used.  If an Id
    // array cannot otherwise be found, the point index is used as the ID.
    vtkSetStringMacro(IdChannelArray);
    vtkGetStringMacro(IdChannelArray);
    
    // Description:
    // Set the maximum aspect ratio allowed for any pair of axes when subdividing
    vtkSetMacro(MaxAspectRatio, double);
    vtkGetMacro(MaxAspectRatio, double);
  
//BTX
    // Description:
    // Return the Bounding Box for a partition
    vtkBoundingBox *GetPartitionBoundingBox(int partition);

    // Description:
    // Return the Bounding Box of all the data (uses MPI to collect all boxes)
    vtkBoundingBox GetGlobalBounds(vtkDataSet *input);
//ETX

    //----------------------------------------------------------------------------
    // Structure to hold all the dataset/mesh/points related data we pass to
    // and from zoltan during the calbacks
    //----------------------------------------------------------------------------
    typedef struct CallbackData {
      vtkZoltanV1PartitionFilter   *self;
      int                           ProcessRank;
      vtkPointSet                  *Input;
      vtkPointSet                  *Output;
      int                           PointType;
      std::vector<int>              ProcessOffsetsPointId;      // offsets into Ids for each process {0, N1, N1+N2, N1+N2+N3...}
      std::vector<int>              ProcessOffsetsCellId;       // offsets into Ids for each process {0, N1, N1+N2, N1+N2+N3...}
      vtkIdType                     InputNumberOfLocalPoints;
      vtkIdType                     InputNumberOfLocalCells;
      vtkIdType                     OutputNumberOfLocalPoints;
      vtkIdType                     OutputNumberOfLocalCells;
      vtkIdType                     OutputNumberOfPointsWithHalo;
      vtkPoints                    *OutputPoints; 
      void                         *InputPointsData;  // float/double
      void                         *OutputPointsData; // float/double
      int                           NumberOfFields;
      int                           MaxCellSize;
      vtkIdType                     OutPointCount;
      vtkIdType                     OutCellCount;
      vtkSmartPointer<vtkCellArray> OutputCellArray;
      std::vector<vtkIdType>        LocalToLocalIdMap;
      std::map<vtkIdType,vtkIdType> ReceivedGlobalToLocalIdMap;
      //
      // The variables below are used twice, once for points, then again for cells
      // but the information is updated in between
      //
      std::vector<void*>            InputArrayPointers;
      std::vector<void*>            OutputArrayPointers;
      std::vector<int>              MemoryPerTuple;
      int                           TotalSizePerId;
    } CallbackData;

    //----------------------------------------------------------------------------
    // Structure holding pointers that zoltan returns for lists of assignments
    // to processors etc
    //----------------------------------------------------------------------------
    typedef struct ZoltanLoadBalanceData {
      int changes, numGidEntries, numLidEntries, numImport, numExport;
      ZOLTAN_ID_PTR importGlobalGids;
      ZOLTAN_ID_PTR importLocalGids; 
      ZOLTAN_ID_PTR exportGlobalGids;
      ZOLTAN_ID_PTR exportLocalGids;
      int *importProcs;
      int *importToPart;
      int *exportProcs;
      int *exportToPart;
    } ZoltanLoadBalanceData;

    //----------------------------------------------------------------------------
    // Structure we use as a temp storage for 
    // to processors etc
    //----------------------------------------------------------------------------
    typedef struct {
      std::vector<ZOLTAN_ID_TYPE> GlobalIds;
      std::vector<ZOLTAN_ID_TYPE> LocalIds;
      std::vector<int> Procs;
    } PartitionInfo;

//BTX
    // Description:
    // zoltan callback to return number of points participating in load/balance
    static int get_number_of_objects_points(void *data, int *ierr);

    // Description:
    // Zoltan callback which fills the Ids for each point in the exchange
    static void get_object_list_points(void *data, int sizeGID, int sizeLID,
      ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr);

    // Description:
    // Zoltan callback which returns the dimension of geometry (3D for us)
    static int get_num_geometry(void *data, int *ierr);

    // Description:
    // Zoltan callback which returns coordinate geometry data (points)
    // templated here to alow float/double instances in our implementation
    template<typename T>
    static void get_geometry_list(
      void *data, int sizeGID, int sizeLID, int num_obj, 
      ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
      int num_dim, double *geom_vec, int *ierr);

    // Description:
    // A ZOLTAN_OBJ_SIZE_FN query function returns the size (in bytes) of the data buffer 
    // that is needed to pack all of a single object's data.
    // Here we add up the size of all the field arrays for points + the geometry itself
    template<typename T>
    static int zoltan_obj_size_func_points(void *data, 
      int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, 
      ZOLTAN_ID_PTR local_id, int *ierr);

    // Description:
    // Zoltan callback to pack all the data for one point into a buffer
    template<typename T>
    static void zoltan_pack_obj_func_points(void *data, int num_gid_entries, int num_lid_entries,
      ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size, char *buf, int *ierr);

    // Description:
    // Zoltan callback to unpack all the data for one point from a buffer
    template<typename T>
    static void zoltan_unpack_obj_func_points(void *data, int num_gid_entries,
      ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr);

    // Description:
    // Zoltan callback for Pre migration setup/initialization
    template<typename T>
    static void zoltan_pre_migrate_func_points(void *data, int num_gid_entries, int num_lid_entries,
      int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
      int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
      ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr);
//ETX

    void        InitializeZoltanLoadBalance();
    static void add_Id_to_interval_map(CallbackData *data, vtkIdType GID, vtkIdType LID);
    vtkIdType   global_to_local_Id(vtkIdType GID);

  protected:
     vtkZoltanV1PartitionFilter();
    ~vtkZoltanV1PartitionFilter();

    int  GatherDataTypeInfo(vtkPoints *points);
    bool GatherDataArrayInfo(vtkDataArray *data, int &datatype, std::string &dataname, int &numComponents);
    void SetupFieldArrayPointers(vtkDataSetAttributes *fields);

    // Override to specify support for vtkPointSet input type.
    virtual int FillInputPortInformation(int port, vtkInformation* info);

    // Override to specify different type of output
    virtual int FillOutputPortInformation(int vtkNotUsed(port), 
      vtkInformation* info);

    // Description:
    virtual int RequestInformation(vtkInformation*,
                            vtkInformationVector**,
                            vtkInformationVector*);
    // Description:
    virtual int RequestUpdateExtent(vtkInformation*, 
                                   vtkInformationVector**, 
                                   vtkInformationVector*);
    
    // Description:
    virtual int RequestData(vtkInformation*,
                            vtkInformationVector**,
                            vtkInformationVector*);
//BTX
    vtkSmartPointer<vtkIdTypeArray> GenerateGlobalIds(vtkIdType Npoints, vtkIdType Ncells, const char *ptidname, vtkIdTypeArray *ptIds);

    vtkSmartPointer<vtkIntArray> BuildCellToProcessList(
      vtkDataSet *data, PartitionInfo &partitioninfo, 
      std::vector<int> &ProcessOffsetsCellId,
      int numExport,
      ZOLTAN_ID_PTR exportGlobalGids ,
      ZOLTAN_ID_PTR exportLocalGids,
      int *exportProcs);

    MPI_Comm GetMPIComm();
    int PartitionPoints(vtkInformation* info, vtkInformationVector** inputVector, vtkInformationVector* outputVector);

    //
    vtkBoundingBox             *LocalBox;
    std::vector<vtkBoundingBox> BoxList;
//ETX
    //
    vtkMultiProcessController   *Controller;
     //
    int                          UpdatePiece;
    int                          UpdateNumPieces;
    char                        *IdChannelArray;
    double                       MaxAspectRatio;
    vtkBoundsExtentTranslator   *ExtentTranslator;
    vtkBoundsExtentTranslator   *InputExtentTranslator;
    vtkSmartPointer<vtkTimerLog> Timer;
    std::string                  IdsName;
    //
    struct Zoltan_Struct       *ZoltanData;
    CallbackData                ZoltanCallbackData;
    ZoltanLoadBalanceData       LoadBalanceData;

    //
    // For debugging
    //
    int pack_count;
    int unpack_count;
    int size_count;

  private:
    vtkZoltanV1PartitionFilter(const vtkZoltanV1PartitionFilter&);  // Not implemented.
    void operator=(const vtkZoltanV1PartitionFilter&);  // Not implemented.
};

#endif
