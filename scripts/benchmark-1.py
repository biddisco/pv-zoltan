#### import the simple module from the paraview
from paraview.simple import *
import paraview.benchmark
import socket, os, sys, re
import psutil

# ------------------------------
# setup machine dependent things
# ------------------------------
hostname = socket.gethostname()
if "carona" in hostname:
  print(hostname + " contains Carona, laptop usage")
  plugin_name = '/Users/biddisco/build/egpgv/bin/libpv_zoltan.dylib'
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
#  dambreak1h5part.UpdatePipeline()

  # show data in view
#  damBreak_Display = Show(dambreak1h5part, renderView1)

  # reset view to fit data
#  renderView1.ResetCamera()
#  RenderAllViews()

#  Hide(damBreak_Display, renderView1)
  # show data in view
  #dambreak1h5partDisplay = Show(dambreak1h5part, renderView1)

  # create a new 'Particle Partition Filter'
  particlePartitionFilter1 = ParticlePartitionFilter(Input=dambreak1h5part)
   # Properties modified on particlePartitionFilter1
  particlePartitionFilter1.WeightsScalarArray = ''
  particlePartitionFilter1.KeepInversePointLists = 0
  particlePartitionFilter1.Maxaspectratiobetweenboundingboxaxes = 0.005

  # show data in view
  particlePartitionFilter1Display = Show(particlePartitionFilter1, renderView1)
  particlePartitionFilter1Display.ColorArrayName = ('POINT_DATA', 'vtkGhostLevels')

  if (do_render):
    # reset view to fit data
    renderView1.ResetCamera()
    
    #### saving camera placements for all active views

    # current camera placement for renderView1
    renderView1.CameraPosition = [2.605978786945343, 0.0, 3.511761331785175]
    renderView1.CameraFocalPoint = [2.605978786945343, 0.0, 0.2749998830695404]
    renderView1.CameraParallelScale = 0.8377355073812321

    #### uncomment the following to render all views
    RenderAllViews()

    # alternatively, if you want to write images, you can use SaveScreenshot(...).
    SaveScreenshot('image' + f + '.png')
  else:
    particlePartitionFilter1.UpdatePipeline()

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
#  print("\nMemory print_logs "),
#  paraview.benchmark.print_logs()
