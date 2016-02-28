#### import the simple module from the paraview
from paraview.simple import *
import paraview.benchmark
import socket, os, sys, re, getopt
import setup_plugins

try:
  opts, args = getopt.getopt(sys.argv[1:],"r:",["resample="])
except getopt.GetoptError:
  print 'test.py -r <reduction factor>'
  sys.exit(2)

paths = setup_plugins.load_plugins()
data_path = paths[0]
output_path = paths[1]

resample = 1
for o, a in opts:
    if o == "-r":
        print("Resample factor " + str(a))
        resample = int(a)
    else:
        assert False, "unhandled option" + str(o) + " " + str(a)

# help(servermanager.vtkProcessModule.GetProcessModule())
nranks = servermanager.vtkProcessModule.GetProcessModule().GetNumberOfLocalPartitions()
rank   = servermanager.vtkProcessModule.GetProcessModule().GetPartitionId()
print ('Number of processes ' + str(nranks) + ' this is rank ' + str(rank))

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# set active view
SetActiveView(None)

filelist = []

for filename in os.listdir(data_path):
  filenumber = int(filename.split("dambreak")[1].split(".")[0])
  #print "File name is " + filename + " number is " + str(filenumber)
  filelist.append(filenumber)

filelist = sorted(filelist)
filelist = ["dambreak"+str(x)+".h5part" for x in filelist]
filelist = filelist[-1:]
print("List of files to work with is :", filelist)

do_render = False

# Convenience method to ask paraview to produce logs with lots of space and highest resolution
paraview.benchmark.maximize_logs()

# create an H5Part reader
dambreak1h5part = H5Part(FileName=data_path + '/' + filelist[0], StepName="Fluid")

for f in filelist:
  #raw_input("Press Enter to continue...")

  print("Setting filename to " + data_path + '/' + f)
  dambreak1h5part.FileName = data_path + '/' + f

  # Properties modified on dambreak1h5part
  dambreak1h5part.Xarray = 'X'
  dambreak1h5part.Yarray = 'Y'
  dambreak1h5part.Zarray = 'Z'
  dambreak1h5part.PointArrays = ['DeltaX', 'ID', 'Kind', 'P', 'VX', 'VY', 'VZ', 'Volume']
  dambreak1h5part.UpdatePipeline()
#  numPoints = dambreak1h5part.GetDataInformation().GetNumberOfPoints()
#  resampled = numPoints / resample
#  ratio = 
#  print("Number of points is " + str(numPoints), "Resample ratio is ", resampled)

  if resample>1:  
    # create a new 'Mask Points'
    maskPoints1 = MaskPoints(Input=dambreak1h5part)
    maskPoints1.MaximumNumberofPoints = 2147483647
    maskPoints1.OnRatio = resample
    maskPoints1.RandomSampling = 0
    maskPoints1.GenerateVertices = 1
    maskPoints1.SingleVertexPerCell = 1
    maskPoints1.ProportionallyDistributeMaximumNumberOfPoints = 1
    maskPoints1.UpdatePipeline()
    # create a new 'H5PartWriter'
    h5PartWriter1 = H5PartWriter(Input=maskPoints1)
  else:
    # create a new 'H5PartWriter'
    h5PartWriter1 = H5PartWriter(Input=dambreak1h5part)

  print("Setting output to " + output_path + '/' + str(resample) + '_' + f)
  h5PartWriter1.StepName = "Step"
  h5PartWriter1.FileName = output_path + '/' + str(resample) + '_' + f;
  h5PartWriter1.UpdatePipeline()

  memory = []
  print("\nMemory use ")
  memuse = paraview.benchmark.get_memuse()
  for s in memuse:
    localmem = s.split()[1]
    print (s + ' ' + localmem)
    memory.append(int(localmem))

  average_mem = sum(memory) / (float(len(memory))*1024.0*1024.0)
  print(memory, " Average memory ", average_mem)

  print("\nMemory parse_logs "),
  paraview.benchmark.parse_logs(show_parse=True, tabular=True)
