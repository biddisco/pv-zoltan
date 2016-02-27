#### import the simple module from the paraview
from paraview.simple import *
import paraview.benchmark
import socket, os, sys, re

# ------------------------------
# setup machine dependent things
# ------------------------------
hostname = socket.gethostname()
if "carona" in hostname:
  print(hostname + " contains Carona, laptop usage")
  plugin_name = '/Users/biddisco/build/egpgv/bin/libpv_zoltan.dylib'
  plugin_name = '/Users/biddisco/build/egpgv/bin/libpv_meshless.dylib'
#  data_path = '/Users/biddisco/data/sphflow/0100millions'
  data_path = '/Users/biddisco/data/sphflow/0001millions/hdf5'
else:
  print("Running on some other machine")
  plugin_name = '/Users/biddisco/build/egpgv/bin/libpv_zoltan.dylib'

### Load pv-zoltan plugin
paraview.servermanager.LoadPlugin(plugin_name)

# help(servermanager.vtkProcessModule.GetProcessModule())
nranks = servermanager.vtkProcessModule.GetProcessModule().GetNumberOfLocalPartitions()
rank   = servermanager.vtkProcessModule.GetProcessModule().GetPartitionId()
print ('Number of processes ' + str(nranks) + ' this is rank ' + str(rank))

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# set active view
SetActiveView(None)

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [640, 480]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.StereoType = 0
renderView1.Background = [0.0, 0.0, 0.0]

# create a new 'H5Part'
filelist = []

for filename in os.listdir(data_path):
  filenumber = int(filename.split("dambreak")[1].split(".")[0])
  #print "File name is " + filename + " number is " + str(filenumber)
  filelist.append(filenumber)

filelist = sorted(filelist)
filelist = ["dambreak"+str(x)+".h5part" for x in filelist]
del filelist[1:]
# print(filelist)

do_render = False

# Convenience method to ask paraview to produce logs with lots of space and highest resolution
paraview.benchmark.maximize_logs()

# create an H5Part reader
dambreak1h5part = H5Part(FileName=data_path + '/' + filelist[0])

for f in filelist:
#  raw_input("Press Enter to continue...")

  print("Setting filename to " + data_path + '/' + f)
  dambreak1h5part.FileName = data_path + '/' + f

  # Properties modified on dambreak1h5part
  dambreak1h5part.Xarray = 'X'
  dambreak1h5part.Yarray = 'Y'
  dambreak1h5part.Zarray = 'Z'
  dambreak1h5part.PointArrays = ['P']

  # create a new 'Mask Points'
  maskPoints1 = MaskPoints(Input=dambreak1h5part)

  # Properties modified on maskPoints1
  maskPoints1.MaximumNumberofPoints = 2147483647
  maskPoints1.GenerateVertices = 1
  maskPoints1.SingleVertexPerCell = 1
  maskPoints1.ProportionallyDistributeMaximumNumberOfPoints = 1
  maskPoints1.UpdatePipeline()

  print("Setting output to " + data_path + '/' + 'resampled_' + f)
  # create a new 'H5PartWriter'
  h5PartWriter1 = H5PartWriter(Input=maskPoints1)
  h5PartWriter1.FileName = data_path + '/' + 'resampled_' + f;
  h5PartWriter1.UpdatePipeline()

  memory = []
  print("\nMemory use ")
  memuse = paraview.benchmark.get_memuse()
  for s in memuse:
    localmem = s.split()[1]
    print (s + ' ' + localmem)
    memory.append(int(localmem))

  average_mem = sum(memory) / (float(len(memory))*1024.0*1024.0)

  print(memory, average_mem)
  print("\nMemory parse_logs "),
  paraview.benchmark.parse_logs(show_parse=True, tabular=True)
