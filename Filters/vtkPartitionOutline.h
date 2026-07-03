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
// .NAME vtkPartitionOutline - Generate outline boxes from partitioned regions of data
//
// .SECTION Description
//
// .SECTION See Also
// 

#ifndef vtkPartitionOutline_h
#define vtkPartitionOutline_h

#include "vtkPVZoltanFiltersModule.h" // for export macro
#include "vtkPolyDataAlgorithm.h"
class vtkPKdTree;

class VTKPVZOLTANFILTERS_EXPORT vtkPartitionOutline : public vtkPolyDataAlgorithm {
  public:
    // Description:
    // Standard Type-Macro
    vtkTypeMacro(vtkPartitionOutline,vtkPolyDataAlgorithm);

    // Description:
    // Create an instance of vtkPartitionOutline
    static vtkPartitionOutline *New();

    // Description:
    // By default, each process generates only the box for its own rank
    // when this option is enabled, each node produces all N boxes which may be
    // useful for some debugging purposes
    vtkSetMacro(AllBoxesOnAllProcesses,int);
    vtkGetMacro(AllBoxesOnAllProcesses,int);
    vtkBooleanMacro(AllBoxesOnAllProcesses,int);

    // Description:
    // Specify a factor (default=1) to increase (>1) or decrease (<1)
    // the generated boxes by.
    vtkSetMacro(InflateFactor,double);
    vtkGetMacro(InflateFactor,double);

    vtkSetMacro(ShowKdTreeBounds,int);
    vtkGetMacro(ShowKdTreeBounds,int);
    vtkBooleanMacro(ShowKdTreeBounds,int);
    void ShowPKdTree(vtkPKdTree *tree, vtkPolyData *output);

  protected:
     vtkPartitionOutline(void);
    ~vtkPartitionOutline() override;
    //
    int RequestInformation (vtkInformation *,
                            vtkInformationVector **,
                            vtkInformationVector *) override;
    //
    int RequestData(vtkInformation *request,
                    vtkInformationVector** inputVector,
                    vtkInformationVector* outputVector) override;
    //
    int FillInputPortInformation(int port, vtkInformation* info) override;

    int    AllBoxesOnAllProcesses;
    int    ShowKdTreeBounds;
    double InflateFactor;
private:
  vtkPartitionOutline(const vtkPartitionOutline&) = delete;
  void operator=(const vtkPartitionOutline&) = delete;
};

#endif

