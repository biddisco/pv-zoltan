#### import the simple module from the paraview
from paraview.simple import *
import paraview.benchmark
import socket, os, sys, re, getopt
import setup_plugins, stats

try:
  opts, args = getopt.getopt(sys.argv[1:],"f:g:p:",["filter=","ghostoverlap=","filepath="])
except getopt.GetoptError:
  print 'test.py -g <ghost overlap>, -f <0=PPF, 1=D3>, -p inputfile'
  sys.exit(2)

paths = setup_plugins.load_plugins()
data_path = paths[0]
output_path = paths[1]
image_path = paths[2]

filepath = ""
ghostoverlap = 0
filter = 0
for o, a in opts:
    if o == "-g":
        print("ghostoverlap " + str(a))
        ghostoverlap = float(a)
    elif o == "-p":
        print("filepath " + str(a))
        filepath = a
    elif o == "-f":
        filter = int(a)
        print("filter ",  'PPF' if (filter==0) else 'D3')
    else:
        assert False, "unhandled option" + str(o) + " " + str(a)

# help(servermanager.vtkProcessModule.GetProcessModule())
nranks = servermanager.vtkProcessModule.GetProcessModule().GetNumberOfLocalPartitions()
rank   = servermanager.vtkProcessModule.GetProcessModule().GetPartitionId()
print ('Number of processes ' + str(nranks) + ' this is rank ' + str(rank))

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
obj = servermanager.misc.GlobalMapperProperties()
obj.GlobalImmediateModeRendering = 1

# set active view
SetActiveView(None)

# create a new 'H5Part'
filelist = []

do_render = True

# Convenience method to ask paraview to produce logs with lots of space and highest resolution
paraview.benchmark.maximize_logs()

# create an H5Part reader
dambreak1h5part = H5Part(FileName=filepath,StepName="Step")

#raw_input("Press Enter to continue...")

# Properties modified on dambreak1h5part
dambreak1h5part.Xarray = 'Coords_0'
dambreak1h5part.Yarray = 'Coords_1'
dambreak1h5part.Zarray = 'Coords_2'
dambreak1h5part.PointArrays = ['DeltaX', 'ID', 'Kind', 'P', 'VX', 'VY', 'VZ', 'Volume']

# create a new 'Particle Partition Filter'
if (filter==0):
  partitionFilter1 = ParticlePartitionFilter(Input=tetrahedralize1)
  partitionFilter1.WeightsScalarArray = ''
  partitionFilter1.KeepInversePointLists = 0
  partitionFilter1.MaxAspectRatio = 5
  partitionFilter1.GhostHaloSize = ghostoverlap
else:
  partitionFilter1 = D3(Input=tetrahedralize1)
  partitionFilter1.MinimalMemory = 1
  partitionFilter1.BoundaryMode = 'Assign cells uniquely'
  
if (do_render):
  # Create a new 'Render View'
  renderView1 = CreateView('RenderView')
  renderView1.ViewSize = [640, 480]
  renderView1.AxesGrid = 'GridAxes3DActor'
  renderView1.StereoType = 0
  renderView1.Background = [0.0, 0.0, 0.0]

  # show data in view
  partitionFilter1Display = Show(partitionFilter1, renderView1)
  partitionFilter1Display.ColorArrayName = ('POINT_DATA', 'vtkGhostLevels')

  # reset view to fit data
  #renderView1.ResetCamera()
  
  # current camera placement for renderView1
  renderView1.CameraPosition = [1.6058090848410411, -7.784277150391815, 1.031715159929654]
  renderView1.CameraFocalPoint = [1.6058090848410411, 0.007768508687273745, 1.031715159929654]
  renderView1.CameraViewUp = [0.0, 0.0, 1.0]
  renderView1.CameraParallelScale = 2.0167298168780916

  #### uncomment the following to render all views
  #RenderAllViews()

  # alternatively, if you want to write images, you can use SaveScreenshot(...).
  imagename = image_path + '/' + os.path.basename(filepath) + '.png'
  print('writing image ', imagename)
  SaveScreenshot(imagename)
  
else:
  # Update
  partitionFilter1.UpdatePipeline()
  pass

numPoints_1 = dambreak1h5part.GetDataInformation().GetNumberOfPoints()
numPoints_2 = partitionFilter1.GetDataInformation().GetNumberOfPoints()
print("Num Dambreak Points ", numPoints_1, "Num PPF Points ", numPoints_2)
print("Num Ghost Points ", numPoints_2-numPoints_1, " ")

stats.dump_stats()
