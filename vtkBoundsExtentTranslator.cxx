/*=========================================================================

  Module                  : vtkBoundsExtentTranslator.cxx

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
#include "vtkBoundsExtentTranslator.h"
#include "vtkObjectFactory.h"
#include "vtkMultiProcessController.h"
#include "vtkMath.h"
#include "vtkBoundingBox.h"
#include "vtkPKdTree.h"

vtkStandardNewMacro(vtkBoundsExtentTranslator);

//----------------------------------------------------------------------------
vtkBoundsExtentTranslator::vtkBoundsExtentTranslator()
{
  this->MaximumGhostDistance = 0;
  this->BoundsHalosEnabled   = 0;
  this->UserBoundsEnabled    = 0;
  this->KdTree               = NULL;
  this->WholeBounds[0] = this->WholeBounds[2] = this->WholeBounds[4] = 0;
  this->WholeBounds[1] = this->WholeBounds[3] = this->WholeBounds[5] = -1;
  this->Spacing[0]     = this->Spacing[1]     = this->Spacing[2]     = 0;
}

//----------------------------------------------------------------------------
vtkBoundsExtentTranslator::~vtkBoundsExtentTranslator()
{
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::ShallowCopy(vtkBoundsExtentTranslator *trans)
{
  // inherited members we copy just in case
  this->Piece          = trans->Piece;
  this->NumberOfPieces = trans->NumberOfPieces;
  this->GhostLevel     = trans->GhostLevel;
  this->SplitMode      = trans->SplitMode;
  memcpy(this->Extent,      trans->Extent,      sizeof(int)*6);
  memcpy(this->WholeExtent, trans->WholeExtent, sizeof(int)*6);
  //
  this->MaximumGhostDistance = trans->MaximumGhostDistance;
  this->BoundsHalosEnabled   = trans->BoundsHalosEnabled;
  this->UserBoundsEnabled    = trans->UserBoundsEnabled;
  this->KdTree               = trans->KdTree; // @TODO warning, smartpointer copy, might need attention
  this->BoundsTable          = trans->BoundsTable;
  this->BoundsTableHalo      = trans->BoundsTableHalo;
  this->BoundsTableUser      = trans->BoundsTableUser;
  memcpy(this->WholeBounds, trans->WholeBounds, sizeof(double)*6);
  memcpy(this->Spacing,     trans->Spacing,     sizeof(double)*3);
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::SetNumberOfPieces(int pieces)
{
  // Allocate a table for this number of pieces.
  this->BoundsTable.resize(pieces*6);
  this->BoundsTableHalo.resize(pieces*6);
  this->BoundsTableUser.resize(pieces*6);
  this->Superclass::SetNumberOfPieces(pieces);
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::ExchangeBoundsForAllProcesses(vtkMultiProcessController *controller, double localbounds[6])
{
  // init bounds storage arrays
  int piece = controller->GetLocalProcessId();
  this->SetNumberOfPieces(controller->GetNumberOfProcesses());
  //
  controller->AllGather(localbounds, &this->BoundsTable[0], 6);
  //
  for (int i=0; i<controller->GetNumberOfProcesses(); i++) {
    this->SetBoundsForPiece(i, &this->BoundsTable[i*6]);
  }
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::SetBoundsForPiece(int piece, double* bounds)
{
  if ((piece*6)>this->BoundsTable.size() || (piece < 0))
    {
    vtkErrorMacro("Piece " << piece << " does not exist.  "
                  "GetNumberOfPieces() is " << this->GetNumberOfPieces());
    return;
    }
  memcpy(&this->BoundsTable[piece*6], bounds, sizeof(double)*6);
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::SetBoundsForPiece(int piece, vtkBoundingBox &box)
{
  double bounds[6];
  box.GetBounds(bounds);
  this->SetBoundsForPiece(piece, bounds);
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::SetBoundsHaloForPiece(int piece, double* bounds)
{
  if ((piece*6)>this->BoundsTableHalo.size() || (piece < 0))
    {
    vtkErrorMacro("Piece " << piece << " does not exist.  "
                  "GetNumberOfPieces() is " << this->GetNumberOfPieces());
    return;
    }
  memcpy(&this->BoundsTableHalo[piece*6], bounds, sizeof(double)*6);
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::SetBoundsHaloForPiece(int piece, vtkBoundingBox &box)
{
  double bounds[6];
  box.GetBounds(bounds);
  this->SetBoundsHaloForPiece(piece, bounds);
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::SetUserBoundsForPiece(int piece, double* bounds)
{
  if ((piece*6)>this->BoundsTableUser.size() || (piece < 0))
    {
    vtkErrorMacro("Piece " << piece << " does not exist.  "
                  "GetNumberOfPieces() is " << this->GetNumberOfPieces());
    return;
    }
  memcpy(&this->BoundsTableUser[piece*6], bounds, sizeof(double)*6);
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::SetUserBoundsForPiece(int piece, vtkBoundingBox &box)
{
  double bounds[6];
  box.GetBounds(bounds);
  this->SetUserBoundsForPiece(piece, bounds);
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::InitWholeBounds()
{
  vtkBoundingBox wholeBounds(&this->BoundsTable[0]);
  for (int p=1; p<this->GetNumberOfPieces(); p++) {
    wholeBounds.AddBounds(&this->BoundsTable[p*6]);
  }
  wholeBounds.GetBounds(WholeBounds);
}

//----------------------------------------------------------------------------
double *vtkBoundsExtentTranslator::GetWholeBounds()
{
  return this->WholeBounds;
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::GetBoundsForPiece(int piece, double* bounds)
{
  if ((piece*6)>this->BoundsTable.size() || (piece < 0))
    {
    vtkErrorMacro("Piece " << piece << " does not exist.  "
                  "GetNumberOfPieces() is " << this->GetNumberOfPieces());
    return;
    }
  memcpy(bounds, &this->BoundsTable[piece*6], sizeof(double)*6);
}

//----------------------------------------------------------------------------
double* vtkBoundsExtentTranslator::GetBoundsForPiece(int piece)
{
  static double emptyBounds[6] = {0,-1,0,-1,0,-1};
  if ((piece*6)>this->BoundsTable.size() || (piece < 0))
    {
    vtkErrorMacro("Piece " << piece << " does not exist.  "
                  "GetNumberOfPieces() is " << this->GetNumberOfPieces());
    return emptyBounds;
    }
  return &this->BoundsTable[piece*6];
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::GetBoundsHaloForPiece(int piece, double* bounds)
{
  if ((piece*6)>this->BoundsTableHalo.size() || (piece < 0))
    {
    vtkErrorMacro("Piece " << piece << " does not exist.  "
                  "GetNumberOfPieces() is " << this->GetNumberOfPieces());
    return;
    }
  memcpy(bounds, &this->BoundsTableHalo[piece*6], sizeof(double)*6);
}

//----------------------------------------------------------------------------
double* vtkBoundsExtentTranslator::GetBoundsHaloForPiece(int piece)
{
  static double emptyBounds[6] = {0,-1,0,-1,0,-1};
  if ((piece*6)>this->BoundsTableHalo.size() || (piece < 0))
    {
    vtkErrorMacro("Piece " << piece << " does not exist.  "
                  "GetNumberOfPieces() is " << this->GetNumberOfPieces());
    return emptyBounds;
    }
  return &this->BoundsTableHalo[piece*6];
}

//----------------------------------------------------------------------------
void vtkBoundsExtentTranslator::SetKdTree(vtkPKdTree *cuts)
{
  this->KdTree = cuts;
}

//----------------------------------------------------------------------------
vtkPKdTree *vtkBoundsExtentTranslator::GetKdTree()
{
  return this->KdTree;
}

//----------------------------------------------------------------------------
// Make sure these inherited methods report an error is anyone is calling them
//----------------------------------------------------------------------------
int vtkBoundsExtentTranslator::PieceToExtentByPoints()
{
  vtkErrorMacro("PieceToExtentByPoints not supported.");
  return 0;
}

//----------------------------------------------------------------------------
int vtkBoundsExtentTranslator::PieceToExtent()
{
  return 
    this->PieceToExtentThreadSafe(this->Piece, this->NumberOfPieces,
                                  this->GhostLevel, this->WholeExtent,
                                  this->Extent, this->SplitMode, 0);
}

//----------------------------------------------------------------------------
#define BOUNDSEXTENT_ISAXPY(a,s,l,y,z) \
  a[0] = s[0]>0 ? (z[0] - y[0])/s[0] : 0; \
  a[1] = s[1]>0 ? (z[1] - y[1])/s[1] : 0; \
  a[2] = s[2]>0 ? (z[2] - y[2])/s[2] : 0;
//----------------------------------------------------------------------------
int vtkBoundsExtentTranslator::BoundsToExtentThreadSafe(
  double *bounds, int *wholeExtent, int *resultExtent)
{
  double scaling[3], lengths[3];
  vtkBoundingBox wholebounds(this->GetWholeBounds());
  vtkBoundingBox piecebounds(bounds);
  wholebounds.GetLengths(lengths);
  const double *origin = wholebounds.GetMinPoint();
  double dims[3] = {
    (double)(wholeExtent[1]-wholeExtent[0]), 
    (double)(wholeExtent[3]-wholeExtent[2]), 
    (double)(wholeExtent[5]-wholeExtent[4]) 
  };

  for (int i=0; i<3; i++) {
    scaling[i] = (dims[i]>0.0) ? lengths[i]/dims[i] : 0.0;
  }

  const double *minpt = piecebounds.GetMinPoint();
  const double *maxpt = piecebounds.GetMaxPoint();

  for (int i=0; i<3; i++) {
    resultExtent[i*2+0] = scaling[i]>0 ? static_cast<int>(0.5+(minpt[i]-origin[i])/scaling[i]) : 0;
    resultExtent[i*2+1] = scaling[i]>0 ? static_cast<int>(0.5+(maxpt[i]-origin[i])/scaling[i]) : 0;
  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkBoundsExtentTranslator::PieceToExtentThreadSafe(
  int piece, int numPieces, 
  int ghostLevel, int *wholeExtent, 
  int *resultExtent, int splitMode, 
  int byPoints)
{
  return this->BoundsToExtentThreadSafe(this->GetBoundsForPiece(piece), wholeExtent, resultExtent);
}


