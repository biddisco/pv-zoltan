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

#include <vector>
#include <map>
//
#include "vtkDataSetAlgorithm.h" // superclass
#include "vtkBoundingBox.h"
#include "zoltan.h"
#include "vtkSmartPointer.h"

class  vtkMultiProcessController;
class  vtkPoints;
class  vtkIdTypeArray;
class  vtkIntArray;
class  vtkBoundsExtentTranslator;
class  vtkPointSet;
class  vtkDataSetAttributes;
class  vtkCellArray;
struct ProcessExchangeVariables;

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
  
    // Description:
    // The thickness of the region between each partition that is used for 
    // ghost cell exchanges. Any particles within this overlap region of another
    // processor will be duplicated on neighbouring processors (possibly multiple times
    // at corner region overlaps)
    vtkSetMacro(GhostCellOverlap, double);
    vtkGetMacro(GhostCellOverlap, double);
        
//BTX
    // Description:
    // Return the Bounding Box for a partition
    vtkBoundingBox *GetPartitionBoundingBox(int partition);

    // Description:
    // Return the Bounding Box for a partition plus the extended region
    // all around the box where ghost cells might be required/present
    vtkBoundingBox *GetPartitionBoundingBoxWithHalo(int partition);
//ETX

    //----------------------------------------------------------------------------
    // Structure to hold mesh data 
    //----------------------------------------------------------------------------
    typedef struct ProcessExchangeVariables {
      int                           ProcessRank;
      vtkPointSet                  *Input;
      vtkPointSet                  *Output;
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
      int                           NumberOfPointFields;
      int                           NumberOfCellFields;
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
    } ProcessExchangeVariables;

  protected:
     vtkZoltanV1PartitionFilter();
    ~vtkZoltanV1PartitionFilter();

    int  GatherDataTypeInfo(vtkPoints *points);
    bool GatherDataArrayInfo(vtkDataArray *data, int &datatype, std::string &dataname, int &numComponents);
    void InitBoundingBoxes(vtkDataSet *input, vtkBoundingBox &box);
    void SetupFieldArrayPointers(vtkDataSetAttributes *fields, ProcessExchangeVariables &mesh, int &NumberOfFields);

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
    vtkSmartPointer<vtkIdTypeArray> GenerateGlobalIds(vtkIdType Npoints, vtkIdType Ncells, const char *ptidname, ProcessExchangeVariables *mesh);

    typedef struct {
      std::vector<ZOLTAN_ID_TYPE> GlobalIds;
      std::vector<ZOLTAN_ID_TYPE> LocalIds;
      std::vector<int> Procs;
    } PartitionInfo;

    void FindPointsInHaloRegions(vtkPoints *pts, vtkIdTypeArray *IdArray, PartitionInfo &ghostinfo);

    vtkSmartPointer<vtkIntArray> BuildCellToProcessList(
      vtkDataSet *data, PartitionInfo &partitioninfo, 
      std::vector<int> &ProcessOffsetsCellId,
      int numExport,
      ZOLTAN_ID_PTR exportGlobalGids ,
      ZOLTAN_ID_PTR exportLocalGids,
      int *exportProcs);

    //
    vtkBoundingBox             *LocalBoxHalo;
    vtkBoundingBox             *LocalBox;
    std::vector<vtkBoundingBox> BoxList;
    std::vector<vtkBoundingBox> BoxListWithHalo;
//ETX
    //
    vtkMultiProcessController  *Controller;
     //
    int                         UpdatePiece;
    int                         UpdateNumPieces;
    char                       *IdChannelArray;
    double                      GhostCellOverlap;
    double                      MaxAspectRatio;
    vtkBoundsExtentTranslator  *ExtentTranslator;
    int                         ExchangePoints;
    int                         ExchangeCells;
    int                         ExchangeHaloPoints;
    //
  private:
    vtkZoltanV1PartitionFilter(const vtkZoltanV1PartitionFilter&);  // Not implemented.
    void operator=(const vtkZoltanV1PartitionFilter&);  // Not implemented.
};

#endif
