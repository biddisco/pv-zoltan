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

#ifndef _vtkPartitionOutline_h
#define _vtkPartitionOutline_h

#include "vtkPolyDataAlgorithm.h"
class vtkPKdTree;

class VTK_EXPORT vtkPartitionOutline : public vtkPolyDataAlgorithm {
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

    void ShowPKdTree(vtkPKdTree *tree, vtkPolyData *output);

  protected:
     vtkPartitionOutline(void);
    ~vtkPartitionOutline();
    //
    virtual int RequestInformation (vtkInformation *,
                                    vtkInformationVector **,
                                    vtkInformationVector *);
    //
    virtual int RequestData(vtkInformation *request,
                            vtkInformationVector** inputVector,
                            vtkInformationVector* outputVector);
    //
    virtual int FillInputPortInformation(int port, vtkInformation* info);

    int    AllBoxesOnAllProcesses;
    double InflateFactor;
private:
  // Not implemented.
  vtkPartitionOutline(const vtkPartitionOutline&);
  void operator=(const vtkPartitionOutline&);
};

#endif

