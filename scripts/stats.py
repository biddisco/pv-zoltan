#### import the simple module from the paraview
from paraview.simple import *
import paraview.benchmark

def dump_stats():
  ###
  ### memory
  ###
  memory = []
  memuse = paraview.benchmark.get_memuse()[1:]
  for s in memuse:
    localmem = s.split()[1]
    # print (s + ' ' + localmem)
    memory.append(int(localmem))
  
  average_mem = sum(memory) / (float(len(memory))*1024.0*1024.0)
  # print(memory)
  print("Average memory " + str(average_mem))
  print("Total memory " + str(average_mem*len(memory)))
  
  ###
  ### timing logs
  ###
  # print("\nMemory parse_logs "),
  #vtkParticlePartitionFilter
  logs = paraview.benchmark.get_logs()
  print('Logs are ', logs)

  paraview.benchmark.parse_logs(show_parse=True, tabular=True)

  # find max time for filter on all nodes
  #
  # grep "Execute vtkParticlePartitionFilter" slurm.out.log | awk '{print $(NF-1)}' | sort -nr | head -n1
  #
  print("grep \"Execute vtkParticlePartitionFilter\" slurm.out | awk '{print $(NF-1)}' | sort -nr | head -n1")
  print("grep \"Execute vtkDistributedDataFilter\" slurm.out | awk '{print $(NF-1)}' | sort -nr | head -n1")
  
