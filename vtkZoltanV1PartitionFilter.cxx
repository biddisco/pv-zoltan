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
#include "vtkZoltanV1PartitionFilter.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetAttributes.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"
#include "vtkIdTypeArray.h"
#include "vtkBoundingBox.h"
#include "vtkMath.h"
#include "vtkPointLocator.h"
#include "vtkCubeSource.h"
#include "vtkInformationDoubleKey.h"
#include "vtkInformationDoubleVectorKey.h"
#include "vtkInformationIntegerKey.h"
//
// For PARAVIEW_USE_MPI
#include "vtkPVConfig.h"
#ifdef PARAVIEW_USE_MPI
  #include "vtkMPI.h"
  #include "vtkMPIController.h"
  #include "vtkMPICommunicator.h"
#endif
#include "vtkDummyController.h"
//
#include "vtkNew.h"
#include "vtkHexahedron.h"
#include "vtkPKdTree.h"
#include "vtkBSPCuts.h"
#include "vtkKdTreeGenerator.h"
//
#include "vtkPKdTree2.h"
#include "vtkBoundsExtentTranslator.h"
#include "vtkZoltanV1PartitionFilter.h"
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>
#include <algorithm>
#include <map>
#include <stack>
#include <iostream>
#include <ostream>
#include <sstream>
#include <iterator>
//
#include "zz_const.h"
#include "rcb.h"

#include "vtkZoltanBasePartitionFilter.txx"

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkZoltanV1PartitionFilter);
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// vtkZoltanV1PartitionFilter :: implementation 
//----------------------------------------------------------------------------
vtkZoltanV1PartitionFilter::vtkZoltanV1PartitionFilter()
{
}
//----------------------------------------------------------------------------
vtkZoltanV1PartitionFilter::~vtkZoltanV1PartitionFilter()
{
}

//----------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::ExecuteZoltanPartition(
    vtkPointSet *output,
    vtkPointSet *input)
{
  //
  // Zoltan can now partition our points.
  // After this returns, we have redistributed points and the Output holds
  // the list of correct points/fields etc for each process
  //
  int zoltan_error = Zoltan_LB_Partition(
      this->ZoltanData,                            // input (all remaining fields are output)
      &this->LoadBalanceData.changes,              // 1 if partitioning was changed, 0 otherwise
      &this->LoadBalanceData.numGidEntries,        // Number of integers used for a global ID
      &this->LoadBalanceData.numLidEntries,        // Number of integers used for a local ID
      &this->LoadBalanceData.numImport,            // Number of vertices to be sent to me
      &this->LoadBalanceData.importGlobalGids,     // Global IDs of vertices to be sent to me
      &this->LoadBalanceData.importLocalGids,      // Local IDs of vertices to be sent to me
      &this->LoadBalanceData.importProcs,          // Process rank for source of each incoming vertex
      &this->LoadBalanceData.importToPart,         // New partition for each incoming vertex
      &this->LoadBalanceData.numExport,            // Number of vertices I must send to other processes
      &this->LoadBalanceData.exportGlobalGids,     // Global IDs of the vertices I must send
      &this->LoadBalanceData.exportLocalGids,      // Local IDs of the vertices I must send
      &this->LoadBalanceData.exportProcs,          // Process to which I send each of the vertices
      &this->LoadBalanceData.exportToPart);        // Partition to which each vertex will belong

  if (zoltan_error != ZOLTAN_OK) {
    printf("Zoltan_LB_Partition NOT OK...\n");
    MPI_Finalize();
    Zoltan_Destroy(&this->ZoltanData);
    exit(0);
  }
}

//----------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::GetZoltanBoundingBoxes(vtkBoundingBox &globalBounds)
{
  //
  // Get bounding boxes from zoltan and set them in the ExtentTranslator
  //
  this->BoxList.clear();
  for (int p = 0; p<this->UpdateNumPieces; p++) {
    double bounds[6];
    int ndim;
    if (ZOLTAN_OK == Zoltan_RCB_Box(this->ZoltanData, p, &ndim, &bounds[0], &bounds[2], &bounds[4], &bounds[1], &bounds[3], &bounds[5])) {
      if (bounds[0] == -DBL_MAX) { bounds[0] = globalBounds.GetMinPoint()[0]; }
      if (bounds[1] == DBL_MAX) { bounds[1] = globalBounds.GetMaxPoint()[0]; }
      if (bounds[2] == -DBL_MAX) { bounds[2] = globalBounds.GetMinPoint()[1]; }
      if (bounds[3] == DBL_MAX) { bounds[3] = globalBounds.GetMaxPoint()[1]; }
      if (bounds[4] == -DBL_MAX) { bounds[4] = globalBounds.GetMinPoint()[2]; }
      if (bounds[5] == DBL_MAX) { bounds[5] = globalBounds.GetMaxPoint()[2]; }
      vtkBoundingBox box(bounds);
      this->BoxList.push_back(box);
      this->ExtentTranslator->SetBoundsForPiece(p, bounds);
    }
  }
  this->ExtentTranslator->InitWholeBounds();
}