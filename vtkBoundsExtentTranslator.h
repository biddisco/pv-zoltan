/*=========================================================================

  Module                  : vtkBoundsExtentTranslator.h

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
// .NAME vtkBoundsExtentTranslator - Extent translation using bounding boxes of data partitions
// .SECTION Description
// vtkBoundsExtentTranslator provides a vtkExtentTranslator that is
// compute the structured extents of data using the bounding boxes stored
// for each process piece. When using a partitioner such as Zoltan, each piece
// had an independent size - based on the image resolution, the extents are computed
// from the bounding boxes. The translator may also be used as a container
// for bounding information for pieces.

#ifndef __vtkBoundsExtentTranslator_h
#define __vtkBoundsExtentTranslator_h

#include "vtkInformation.h"
#include "vtkInformationObjectMetaDataKey.h"
#include "vtkExtentTranslator.h"
#include "vtkSmartPointer.h"
#include <vector>

class vtkMultiProcessController;
class vtkPKdTree;
class vtkBoundingBox;

class VTK_EXPORT vtkBoundsExtentTranslator : public vtkExtentTranslator
{
public:
  vtkTypeMacro(vtkBoundsExtentTranslator,vtkExtentTranslator);
  
  static vtkInformationObjectMetaDataKey* META_DATA();
  static vtkBoundsExtentTranslator* New();

  virtual void ShallowCopy(vtkBoundsExtentTranslator *trans);

  // Description:
  // Set the bounding box to be used for a piece.
  virtual void SetKdTree(vtkPKdTree *cuts);
  virtual vtkPKdTree *GetKdTree();  

  // Description:
  // Set the number of pieces that will be stored.
  // Executives set this and the value must match that set by the
  // PartitionFilter that generated this data, otherwise failure
  // is probable
  virtual void SetNumberOfPieces(int pieces);

  // Description:
  // Used by executives to compute a structured extent.
  // Only valid when a default Spacing has been set, otherwise
  // all pieces return the WholeExtent
  virtual int PieceToExtent();  
  virtual int PieceToExtentByPoints();
  virtual int BoundsToExtentThreadSafe(
      double *bounds, int *wholeExtent, int *resultExtent);
  virtual int PieceToExtentThreadSafe(int piece, int numPieces, 
                              int ghostLevel, int *wholeExtent, 
                              int *resultExtent, int splitMode, 
                              int byPoints);

  // Description:
  // Set the bounding box to be used for a piece.
  virtual void SetBoundsForPiece(int piece, double* bounds);
  
  // Description:
  // Set the bounding box to be used for a piece.
  virtual void SetBoundsForPiece(int piece, vtkBoundingBox &box);
  
  // Description:
  // Set the bounding box + halo region to be used for a piece.
  virtual void SetBoundsHaloForPiece(int piece, double* bounds);

  // Description:
  // Set the bounding box + halo region to be used for a piece.
  virtual void SetBoundsHaloForPiece(int piece, vtkBoundingBox &box);

  // Description:
  // If a modified bounding box table is required, the bounds
  // must be set for each piece and UserBoundsEnabled must be set
  void SetUserBoundsForPiece(int piece, double* bounds);
  void SetUserBoundsForPiece(int piece, vtkBoundingBox &box);

  // Description:
  // Set the maximum ghost overlap region that is required 
  vtkSetMacro(UserBoundsEnabled, int);
  vtkGetMacro(UserBoundsEnabled, int); 
  vtkBooleanMacro(UserBoundsEnabled, int); 

  // Description:
  // Set the maximum ghost overlap region that is required 
  vtkSetMacro(BoundsHalosEnabled, int);
  vtkGetMacro(BoundsHalosEnabled, int); 
  vtkBooleanMacro(BoundsHalosEnabled, int); 

  // Description:  
  // Get the bounds table entry for the given piece.  
  // Structured Extents should be calculated using PieceToExtent
  // (after a valid spacing has been set).
  virtual void    GetBoundsForPiece(int piece, double *bounds);
  virtual double *GetBoundsForPiece(int piece);
  virtual void    GetBoundsHaloForPiece(int piece, double *bounds);
  virtual double *GetBoundsHaloForPiece(int piece);
  
  // Description:
  // Set the maximum ghost overlap region that is required 
  vtkSetMacro(MaximumGhostDistance, double);
  vtkGetMacro(MaximumGhostDistance, double);
  
  // Description:
  // To convert a bounding box to a StructuredExtent, a spacing
  // (between samples) is necessary - and the global bounds.
  vtkSetVector3Macro(Spacing, double);
  vtkGetVector3Macro(Spacing, double);

  // Description:
  // After setting all bounds, call this to compute the global/whole bounds
  virtual void InitWholeBounds();

  // Description:
  // Returns the union of all bounds
  virtual double *GetWholeBounds();

  // Description:
  // if you only know the local bounds and not the bounds of all processes, 
  // call ExchangeBoundsForAllProcesses passing in a bounds object for the local
  // piece and a controller. The class will initialize the correct number of pieces and
  // bounds for all pieces.
  void ExchangeBoundsForAllProcesses(vtkMultiProcessController *controller, double localbounds[6]);

protected:
   vtkBoundsExtentTranslator();
  ~vtkBoundsExtentTranslator();
  
  // Store the extent table in a single array.  Every 6 values form an extent.
  std::vector<double>         BoundsTable;
  int                         BoundsHalosEnabled;
  std::vector<double>         BoundsTableHalo;
  int                         UserBoundsEnabled;
  std::vector<double>         BoundsTableUser;
  double                      WholeBounds[6];
  double                      Spacing[3];
  double                      MaximumGhostDistance;
  vtkSmartPointer<vtkPKdTree> KdTree;
   
private:
  vtkBoundsExtentTranslator(const vtkBoundsExtentTranslator&);  // Not implemented.
  void operator=(const vtkBoundsExtentTranslator&);  // Not implemented.
};

#endif
