//----------------------------------------------------------------------------
//template<typename T>
//std::ostream& PrintVector(std::ostream& out, int width, std::vector<T> &vec) {
//  out << "[" << std::setprecision(1) << std::fixed;
//  for (auto & x : vec) out << std::setw(width) << x << ", ";
//  return out << "]";
//}

//----------------------------------------------------------------------------
// Zoltan callback which fills the Ids for each object in the exchange
//----------------------------------------------------------------------------
void vtkZoltanV1PartitionFilter::get_object_list_points(void *data, int sizeGID, int sizeLID,
  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID, int wgt_dim, float *obj_wgts, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);
  //
  // Return the IDs of our objects, but no weights.
  // Zoltan will assume equally weighted objects.
  //
  vtkIdType N = callbackdata->Input->GetNumberOfPoints();
  for (int i=0; i<N; i++){
    globalID[i] = i + callbackdata->ProcessOffsetsPointId[callbackdata->ProcessRank];
  }
  *ierr = ZOLTAN_OK;
}
//----------------------------------------------------------------------------
// Function Type: Pre migration callback for halo/other particle exchange
// The difference between this and the standard pre_migrate function for points 
// is that we only add new points and do not remove any.
//----------------------------------------------------------------------------
template <typename T>
void vtkZoltanV1PartitionFilter::zoltan_pre_migrate_function_points_add(
  void *data, int num_gid_entries, int num_lid_entries,
  int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
  int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr)
{
  CallbackData *callbackdata = static_cast<CallbackData*>(data);

  // this is the current actual size of the output point list
  vtkIdType OutputNumberOfPoints = callbackdata->Output->GetNumberOfPoints();
  callbackdata->OutPointCount = OutputNumberOfPoints;
  // resize points to accept ghost cell additions
  OutputNumberOfPoints = OutputNumberOfPoints + num_import;
  callbackdata->Output->GetPoints()->GetData()->Resize(OutputNumberOfPoints);
  callbackdata->Output->GetPoints()->SetNumberOfPoints(OutputNumberOfPoints);
  callbackdata->OutputPointsData = (T*)(callbackdata->Output->GetPoints()->GetData()->GetVoidPointer(0));
  // 
  // The migration taking place might be using the original points (standard send/receive)
  // or the already migrated points (like a halo exchange). When doing a halo exchange, we might
  // be resending points we have just received, which are in our output list and not our
  // input list, so change the pointers accordingly
  //
  vtkPointData *inPD;
  callbackdata->InputPointsData = callbackdata->Input->GetPoints() ? (T*)(callbackdata->Input->GetPoints()->GetData()->GetVoidPointer(0)) : NULL;
  inPD  = callbackdata->Input->GetPointData();

  vtkPointData *outPD = callbackdata->Output->GetPointData();
  callbackdata->self->InitializeFieldDataArrayPointers(callbackdata, inPD, outPD, OutputNumberOfPoints);
}
