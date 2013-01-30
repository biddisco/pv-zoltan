/*=========================================================================

  Module                  : vtkPartitionOutline.cxx

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
#include "vtkPartitionOutline.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkIntArray.h"
#include "vtkPolyData.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkBoundingBox.h"
//
#include "vtkBoundsExtentTranslator.h"
#include "vtkAppendPolyData.h"
#include "vtkOutlineSource.h"
#include "vtkPKdTree2.h"
//
#include <cmath>
//---------------------------------------------------------------------------
vtkStandardNewMacro(vtkPartitionOutline);
//---------------------------------------------------------------------------
vtkPartitionOutline::vtkPartitionOutline(void)
{
  this->AllBoxesOnAllProcesses = 0;
  this->ShowKdTreeBounds       = 0;
  this->InflateFactor          = 1.0;
}
//---------------------------------------------------------------------------
vtkPartitionOutline::~vtkPartitionOutline(void) {
}
//----------------------------------------------------------------------------
int vtkPartitionOutline::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkPartitionOutline::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);
  return 1;
}
//----------------------------------------------------------------------------
int vtkPartitionOutline::RequestData(vtkInformation *request,
                                       vtkInformationVector** inputVector,
                                       vtkInformationVector* outputVector)
{
  // get the info objects
  vtkInformation *inInfo  = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  //
  vtkPointSet *input  = vtkPointSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  //
  vtkExtentTranslator *translator = inInfo ? vtkExtentTranslator::SafeDownCast(
    inInfo->Get(vtkStreamingDemandDrivenPipeline::EXTENT_TRANSLATOR())) : NULL;
  vtkBoundsExtentTranslator *bet = vtkBoundsExtentTranslator::SafeDownCast(translator);
  // if the extent translator has not been initialized well - don't use it
  if (bet && bet->GetNumberOfPieces()==0) {
    bet = NULL;
  }

  if (bet && this->ShowKdTreeBounds) {
    this->ShowPKdTree(bet->GetKdTree(), output);
    return 1;
  }


  //
  vtkSmartPointer<vtkAppendPolyData> polys = vtkSmartPointer<vtkAppendPolyData>::New();
  int piece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  size_t boxes = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());;
  if (bet) {
    boxes = bet->GetNumberOfPieces();  
  }
  vtkSmartPointer<vtkIntArray> processIds = vtkSmartPointer<vtkIntArray>::New();
  processIds->SetName("ProcessId");
  vtkSmartPointer<vtkIdTypeArray> quantity = vtkSmartPointer<vtkIdTypeArray>::New();
  quantity->SetName("Occupation");
  //
  vtkBoundingBox box;
  for (int i=0; i<boxes; i++) {
    for (int h=0; h<2; h++) {
      bool add = false;
      if (this->AllBoxesOnAllProcesses || i==piece) {
        if (h==0) {
          if (bet) box.SetBounds(bet->GetBoundsForPiece(i));
          else box.SetBounds(input->GetBounds());
          add = true;
        }
        else if (h==1 && bet && bet->GetBoundsHalosPresent()) {
          box.SetBounds(bet->GetBoundsHaloForPiece(i));
          add = true;
        }
      }
      if (add) {
        double p1[3],p2[3];
        box.GetMinPoint(p1[0], p1[1], p1[2]);
        box.GetMaxPoint(p2[0], p2[1], p2[2]);
        for (int j=0; j<3; j++) {
          double l = box.GetLength(j);
          double d = (l-l*this->InflateFactor);
          p1[j] += d/2.0; 
          p2[j] -= d/2.0; 
        }
        vtkSmartPointer<vtkOutlineSource> cube = vtkSmartPointer<vtkOutlineSource>::New();
        cube->SetBounds(p1[0],p2[0],p1[1],p2[1],p1[2],p2[2]);
        cube->Update();
        polys->AddInputData(cube->GetOutput());
        for (int p=0; p<8; p++) processIds->InsertNextValue(i);
        vtkIdType number = input->GetNumberOfPoints();
        for (int p=0; p<8; p++) quantity->InsertNextValue(number);
      }
    }
  }
  polys->Update();
  output->SetPoints(polys->GetOutput()->GetPoints());
  output->SetLines(polys->GetOutput()->GetLines());
  output->GetPointData()->AddArray(processIds);
  output->GetPointData()->AddArray(quantity);

  return 1;
}
//---------------------------------------------------------------------------
void vtkPartitionOutline::ShowPKdTree(vtkPKdTree *tree, vtkPolyData *output)
{
  if (vtkPKdTree2::SafeDownCast(tree)) {
    vtkPKdTree2::SafeDownCast(tree)->SetInflateFactor(this->InflateFactor);
  }
  tree->SetGenerateRepresentationUsingDataBounds(1);
  tree->GenerateRepresentation(-1, output);
}
//----------------------------------------------------------------------------

