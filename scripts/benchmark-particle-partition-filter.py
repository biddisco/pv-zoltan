#### import the simple module from the paraview
from paraview.simple import *
import paraview.benchmark
import socket, os, sys, re, getopt
import setup_plugins

try:
  opts, args = getopt.getopt(sys.argv[1:],"g:p:",["ghostoverlap=","filepath="])
except getopt.GetoptError:
  print 'test.py -g <ghost overlap>'
  sys.exit(2)

paths = setup_plugins.load_plugins()
data_path = paths[0]
output_path = paths[1]

filepath = ""
ghostoverlap = 0
for o, a in opts:
    if o == "-g":
        print("ghostoverlap " + str(a))
        ghostoverlap = float(a)
    elif o == "-p":
        print("filepath " + str(a))
        filepath = a
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

# create a new 'H5Part'
filelist = []

do_render = False

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
dambreak1h5part.UpdatePipeline()

# create a new 'Particle Partition Filter'
particlePartitionFilter1 = ParticlePartitionFilter(Input=dambreak1h5part)
 # Properties modified on particlePartitionFilter1
particlePartitionFilter1.WeightsScalarArray = ''
particlePartitionFilter1.KeepInversePointLists = 0
particlePartitionFilter1.Maxaspectratiobetweenboundingboxaxes = 5

if (do_render):
  # Create a new 'Render View'
  renderView1 = CreateView('RenderView')
  renderView1.ViewSize = [640, 480]
  renderView1.AxesGrid = 'GridAxes3DActor'
  renderView1.StereoType = 0
  renderView1.Background = [0.0, 0.0, 0.0]

# show data in view
#  damBreak_Display = Show(dambreak1h5part, renderView1)

# reset view to fit data
#  renderView1.ResetCamera()
#  RenderAllViews()

#  Hide(damBreak_Display, renderView1)
# show data in view

  # reset view to fit data
  renderView1.ResetCamera()
  
# show data in view
#particlePartitionFilter1Display = Show(particlePartitionFilter1, renderView1)
#particlePartitionFilter1Display.ColorArrayName = ('POINT_DATA', 'vtkGhostLevels')

  # current camera placement for renderView1
  renderView1.CameraPosition = [2.605978786945343, 0.0, 3.511761331785175]
  renderView1.CameraFocalPoint = [2.605978786945343, 0.0, 0.2749998830695404]
  renderView1.CameraParallelScale = 0.8377355073812321

  #### uncomment the following to render all views
  RenderAllViews()

  # alternatively, if you want to write images, you can use SaveScreenshot(...).
  SaveScreenshot('image' + f + '.png')
  
else:
  #particlePartitionFilter1.UpdatePipeline()
  pass

memory = []
print("\nMemory use (dropping client stats ")
memuse = paraview.benchmark.get_memuse()[1:]
for s in memuse:
  localmem = s.split()[1]
  print (s + ' ' + localmem)
  memory.append(int(localmem))

average_mem = sum(memory) / (float(len(memory))*1024.0*1024.0)
print(memory)
print("Average memory ", average_mem, " len ", len(memory))

print("\nMemory parse_logs "),
paraview.benchmark.parse_logs(show_parse=True, tabular=True)
