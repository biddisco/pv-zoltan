#### import the simple module from the paraview
from paraview.simple import *
import paraview.benchmark
import socket, os, sys, re, getopt, math
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

paraview.simple._DisableFirstRenderCameraReset()


##############
# Pipeline
##############

numcells = 150*150*150*nranks
sidelength = 1 + int(math.pow(numcells, 1.0/3.0)/2)
print("Wavelet size {-" + str(sidelength) + "," + str(sidelength) + "}") 
print("Wavelet numcells " + str(8*sidelength*sidelength*sidelength)) 
waveletside = int(sidelength/2)

# create a new 'Wavelet'
wavelet1 = Wavelet()
# Properties modified on wavelet1
wavelet1.WholeExtent = [-waveletside, waveletside, -waveletside, waveletside, -waveletside, waveletside]

# create a new 'Tetrahedralize'
tetrahedralize1 = Tetrahedralize(Input=wavelet1)

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

# create a new 'Shrink'
shrink1 = Shrink(Input=partitionFilter1)
# Properties modified on shrink1
shrink1.ShrinkFactor = 0.9


##################
# colour table
##################

# get color transfer function/color map for 'RTData'
rTDataLUT = GetColorTransferFunction('RTData')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
rTDataLUT.ApplyPreset('Black-Body Radiation', True)

# get opacity transfer function/opacity map for 'RTData'
rTDataPWF = GetOpacityTransferFunction('RTData')
# Properties modified on rTDataPWF
rTDataPWF.Points = [37.35310363769531, 0.0, 0.5, 0.0, 100.0, 0.0, 0.5, 0.0, 276.8288269042969, 1.0, 0.5, 0.0]

# Properties modified on rTDataLUT
rTDataLUT.EnableOpacityMapping = 1


##################
# rendering
##################

renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [640, 480]

# show data in view
shrink1Display = Show(shrink1, renderView1)
# trace defaults for the display properties.
shrink1Display.ColorArrayName = ['POINTS', 'RTData']
shrink1Display.LookupTable = rTDataLUT

# reset view to fit data
#renderView1.ResetCamera()

# show color bar/color legend
shrink1Display.SetScalarBarVisibility(renderView1, True)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-129.07375961872597, 417.62235377155565, 143.9824442336569]
renderView1.CameraViewUp = [0.7306437409080044, 0.4132080413165191, -0.5435244598574407]
renderView1.CameraParallelScale = 121.11167448863597

#### uncomment the following to render all views
#RenderAllViews()
#alternatively, if you want to write images, you can use 
SaveScreenshot("wavelet-1.png")
#
stats.dump_stats()
