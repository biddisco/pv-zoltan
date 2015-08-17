#define VTK_ZOLTAN2_PARTITION_FILTER

#ifdef VTK_ZOLTAN2_PARTITION_FILTER
  #include "vtkZoltanV2PartitionFilter.h" // superclass
  #define VTK_ZOLTAN_PARTITION_FILTER vtkZoltanV2PartitionFilter
#else
  #include "vtkZoltanV1PartitionFilter.h" // superclass
  #define VTK_ZOLTAN_PARTITION_FILTER vtkZoltanV1PartitionFilter
#endif
