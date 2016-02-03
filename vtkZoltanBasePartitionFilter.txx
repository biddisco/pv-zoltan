/*=========================================================================

  Module : vtkZoltanBasePartitionFilter

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

// templated functions that need to be instantiated 

//----------------------------------------------------------------------------
// Zoltan callback which returns coordinate geometry data (points)
// templated here to alow float/double instances in our implementation
//----------------------------------------------------------------------------
template<typename T>
void vtkZoltanBasePartitionFilter::get_geometry_list(
  void *data, int sizeGID, int sizeLID, int num_obj, 
  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
  int num_dim, double *geom_vec, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  for (int i=0;  i<num_obj; i++){
    geom_vec[3*i]   = ((T*)(callbackdata->InputPointsData))[3*i+0];
    geom_vec[3*i+1] = ((T*)(callbackdata->InputPointsData))[3*i+1];
    geom_vec[3*i+2] = ((T*)(callbackdata->InputPointsData))[3*i+2];
  }
  *ierr = ZOLTAN_OK;
  return;
}

//----------------------------------------------------------------------------
// A ZOLTAN_OBJ_SIZE_FN query function returns the size (in bytes) of the data buffer 
// that is needed to pack all of a single object's data.
//
// Here we add up the size of all the field arrays for points + the geometry itself
//  
// Function Type:   ZOLTAN_OBJ_SIZE_FN_TYPE
// Arguments:   
//  data             Pointer to user-defined data.
//  num_gid_entries  The number of array entries used to describe a single global ID.  
//  num_lid_entries  The number of array entries used to describe a single local ID.  
//  global_id        Pointer to the global ID of the object.
//  local_id         Pointer to the local ID of the object.
//  ierr             Error code to be set by function.
// Returned Value:   
//  int              The size (in bytes) of the required data buffer (per object).
//----------------------------------------------------------------------------
template<typename T>
int vtkZoltanBasePartitionFilter::zoltan_obj_size_function_pointdata(void *data, 
  int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, 
  ZOLTAN_ID_PTR local_id, int *ierr)
{
  INC_SIZE_COUNT
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  *ierr = ZOLTAN_OK;
  return callbackdata->TotalSizePerId + sizeof(T)*3;
}

//----------------------------------------------------------------------------
// Zoltan callback to pack all the data for one point into a buffer
//----------------------------------------------------------------------------
template<typename T>
void vtkZoltanBasePartitionFilter::zoltan_pack_obj_function_pointdata(void *data, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int dest, int size, char *buf, int *ierr)
{
  INC_PACK_COUNT
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  vtkIdType GID = *global_id;
  vtkIdType LID = GID - callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];
  //
  for (int i=0; i<callbackdata->NumberOfFields; i++) {
    int asize = callbackdata->MemoryPerTuple[i];
    char *dataptr = (char*)(callbackdata->InputArrayPointers[i]) + asize*LID;
    memcpy(buf, dataptr, asize);
    buf += asize;
  }
  memcpy(buf, &((T*)(callbackdata->InputPointsData))[LID*3], sizeof(T)*3);  
  *ierr = ZOLTAN_OK;
  return;
}

//----------------------------------------------------------------------------
// Zoltan callback to unpack all the data for one point from a buffer
//----------------------------------------------------------------------------
template<typename T>
void vtkZoltanBasePartitionFilter::zoltan_unpack_obj_function_pointdata(void *data, int num_gid_entries,
  ZOLTAN_ID_PTR global_id, int size, char *buf, int *ierr)
{
  INC_UNPACK_COUNT
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  vtkIdType GID = *global_id;
  //
  for (int i=0; i<callbackdata->NumberOfFields; i++) {
    int asize = callbackdata->MemoryPerTuple[i];
    char *dataptr = (char*)(callbackdata->OutputArrayPointers[i]) + asize*(callbackdata->OutPointCount);
    memcpy(dataptr, buf, asize);
    buf += asize;
  }
//  if (callbackdata->self->UpdatePiece==2 && GID <100) { std::cout <<"Received " << GID << std::endl; }
  add_Id_to_interval_map(callbackdata, GID, callbackdata->OutPointCount);
  memcpy(&((T*)(callbackdata->OutputPointsData))[callbackdata->OutPointCount*3], buf, sizeof(T)*3);  
  callbackdata->OutPointCount++;
  *ierr = ZOLTAN_OK;
  return;
}

//----------------------------------------------------------------------------
// Function Type: Pre migration callback
// Arguments:   
//  data              Pointer to user-defined data.
//  num_gid_entries   The number of array entries used to describe a single global ID.  
//  num_lid_entries   The number of array entries used to describe a single local ID.  
//  num_import        The number of objects that will be received by this processor.
//  import_global_ids An array of num_import global IDs of objects to be received by this processor. 
//                    may be NULL, as the processor does not necessarily need to know which objects it will receive.
//  import_local_ids  An array of num_import local IDs of objects to be received by this processor. 
//                    may be NULL, as the processor does not necessarily need to know which objects it will receive.
//  import_procs      An array of size num_import listing the processor IDs of the source processors. 
//                    may be NULL, as the processor does not necessarily need to know which objects is will receive.
//  import_to_part    An array of size num_import listing the parts to which objects will be imported. 
//                    may be NULL, as the processor does not necessarily need to know which objects it will receive.
//  num_export        The number of objects that will be sent from this processor to other processors.
//  export_global_ids An array of num_export global IDs of objects to be sent from this processor.
//  export_local_ids  An array of num_export local IDs of objects to be sent from this processor.
//  export_procs      An array of size num_export listing the processor IDs of the destination processors.
//  export_to_part    An array of size num_export listing the parts to which objects will be sent.
//  ierr              Error code to be set by function.
template<typename T>
void vtkZoltanBasePartitionFilter::zoltan_pre_migrate_function_points(
  void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  // newTotal = original points - sent away + received
  vtkIdType N  = callbackdata->Input->GetNumberOfPoints();
  vtkIdType N2 = N + num_import - num_export;
  callbackdata->Output->GetPoints()->SetNumberOfPoints(N2);
  callbackdata->OutputPointsData = callbackdata->Output->GetPoints()->GetData()->GetVoidPointer(0);
  vtkPointData    *inPD  = callbackdata->Input->GetPointData();
  vtkPointData    *outPD = callbackdata->Output->GetPointData();
  outPD->CopyAllocate(inPD, N2);
  //
  callbackdata->self->InitializeFieldDataArrayPointers(callbackdata, inPD, outPD, N2);

  // some points are being sent away, some will be received, we must copy
  // the ones that are not moving from the input to the output.
  // Mark points so we know which local points will still be local after the exchange
  callbackdata->LocalToLocalIdMap.assign(N, 0);
  for (vtkIdType i=0; i<num_export; i++) {
    vtkIdType GID = export_global_ids[i];
    vtkIdType LID = GID - callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];    
    callbackdata->LocalToLocalIdMap[LID] = -1;    
  }

  // Loop over each local point and copy it to the output.
  // WARNING: point Ids are changing so any cells referencing the points
  // must have their Ids updated to the new index - create an IdMap to hold this info.
  callbackdata->OutPointCount = 0;
  for (vtkIdType i=0; i<N; i++) {
    if (callbackdata->LocalToLocalIdMap[i]==0) {
      outPD->CopyData(inPD, i, callbackdata->OutPointCount);
      memcpy(&((T*)(callbackdata->OutputPointsData))[callbackdata->OutPointCount*3], &((T*)(callbackdata->InputPointsData))[i*3], sizeof(T)*3);
      callbackdata->LocalToLocalIdMap[i] = callbackdata->OutPointCount;
      callbackdata->OutPointCount++;
    }
  }
  *ierr = ZOLTAN_OK;
}

//----------------------------------------------------------------------------
template<typename T>
void vtkZoltanBasePartitionFilter::CopyPointsToSelf(
  std::vector<vtkIdType> &LocalPointsToKeep, vtkIdType num_reserved,
  void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  // newTotal = original points - sent away + received
  vtkIdType N  = callbackdata->Input->GetNumberOfPoints();
//    if (this->UpdatePiece==0)
//    for (int i=0 ; i<num_export; i++) {
//        cout<<i<<"...."<<export_global_ids[i]<<"\t"<<export_procs[i]<<endl;
//    }
//    cout<<"Smooth"<<endl;
  //
  // our size estimates of the final number of points can be messed up because the list of points being sent away
  // contains some points which are sent to multiple remote processes. We can't therefore use the size of this list
  // to subtract away from the original points to get final point list size.
  // We therefore first need to count the number of unique points being sent away
  // We will do that whilst we ...
  // Mark points to be moved with -1 so we know which local points will still be local after the exchange
  callbackdata->LocalToLocalIdMap.assign(N, 0);
  vtkIdType uniqueSends = 0;
  for (vtkIdType i=0; i<num_export; i++) {
    vtkIdType GID = export_global_ids[i];
    vtkIdType LID = GID - callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];
    if (callbackdata->LocalToLocalIdMap[LID]==0) {
      callbackdata->LocalToLocalIdMap[LID] = -1;
      uniqueSends++;
    }
  }


  // now compute the final number of points we'll have
  vtkIdType N2 = N + num_reserved + num_import - (uniqueSends - LocalPointsToKeep.size());
    this->Controller->Barrier();
    
//    cout<<this->UpdatePiece
//    <<"\tN2:"<<N2
//    <<"\tN:"<<N
//    <<"\tnum_reserved:"<<num_reserved
//    <<"\tnum_import:"<<num_import
//    <<"\tuniqueSends:"<<uniqueSends
//    <<"\tLocalPointsToKeep:"<<LocalPointsToKeep.size()
//    <<endl;
  callbackdata->Output->GetPoints()->SetNumberOfPoints(N2);
  callbackdata->OutputPointsData = callbackdata->Output->GetPoints()->GetData()->GetVoidPointer(0);
  vtkPointData    *inPD  = callbackdata->Input->GetPointData();
  vtkPointData    *outPD = callbackdata->Output->GetPointData();
  outPD->CopyAllOn();
  outPD->CopyAllocate(inPD, N2);
  //
  // prepare for copying data by setting up pointers to field arrays
  //
  callbackdata->self->InitializeFieldDataArrayPointers(callbackdata, inPD, outPD, N2);
  // Loop over each local point and copy it to the output.
  // WARNING: point Ids are changing so any cells referencing the points
  // must have their Ids updated to the new index - create an IdMap to hold this info.
  callbackdata->OutPointCount = 0;
  for (vtkIdType i=0; i<N; i++) {
    // for each point that is staying on this process
    if (callbackdata->LocalToLocalIdMap[i]==0) { 
      outPD->CopyData(inPD, i, callbackdata->OutPointCount);
      memcpy(
        &((T*)(callbackdata->OutputPointsData))[callbackdata->OutPointCount*3],
        &((T*)(callbackdata->InputPointsData))[i*3],
        sizeof(T)*3);
      callbackdata->LocalToLocalIdMap[i] = callbackdata->OutPointCount;
      callbackdata->OutPointCount++;
    }
  }
  vtkIdType maxInitialId = callbackdata->OutPointCount;
  // ...and the points marked for sending which we also need a local copy of 
  for (vtkIdType i=0; i<LocalPointsToKeep.size(); i++) {
    vtkIdType LID = LocalPointsToKeep[i];
    if (callbackdata->LocalToLocalIdMap[LID]==-1) { // the point was marked as moving, but we need it here too
      outPD->CopyData(inPD, LID, callbackdata->OutPointCount);
      memcpy(
        &((T*)(callbackdata->OutputPointsData))[callbackdata->OutPointCount*3],
        &((T*)(callbackdata->InputPointsData))[LID*3],
        sizeof(T)*3);
      callbackdata->LocalToLocalIdMap[LID] = callbackdata->OutPointCount;
      callbackdata->OutPointCount++;
    }
    else if (callbackdata->LocalToLocalIdMap[LID]>=maxInitialId) {
      vtkErrorMacro(<<callbackdata->ProcessRank << " Point already mapped from " << LID << " to " << callbackdata->LocalToLocalIdMap[LID]);
    }
    else {
      vtkErrorMacro(<<callbackdata->ProcessRank << " Serious Error : Point already mapped from " << LID << " to " << callbackdata->LocalToLocalIdMap[LID]);
    }
  }
  // callbackdata->OutPointCount<N2 is allowed as we may receive points, but > is forbidden
  if (callbackdata->OutPointCount>N2) {
    vtkErrorMacro(<<"Serious Error : Point allocation N2 " << N2 << " greater than " << callbackdata->OutPointCount);
  }
}
//----------------------------------------------------------------------------
