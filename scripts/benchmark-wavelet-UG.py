#### import the simple module from the paraview
from paraview.simple import *
import paraview.benchmark
import socket, os, sys, re, getopt, math
import setup_plugins, stats

try:
  opts, args = getopt.getopt(sys.argv[1:],"f:g:r:",["filter=","ghostoverlap=","res"])
except getopt.GetoptError:
  print 'test.py -f <0=PPF, 1=D3>, -g <ghost>, -r <resolution>'
  sys.exit(2)

paths = setup_plugins.load_plugins()
data_path = paths[0]
output_path = paths[1]
image_path = paths[2]

resolution = 100
ghostoverlap = 0
filter = 0
for o, a in opts:
    if o == "-f":
        filter = int(a)
        print("filter ",  'PPF' if (filter==0) else 'D3')
    elif o == "-g":
        print("ghostoverlap " + str(a))
        ghostoverlap = float(a)
    elif o == "-r":
        print("resolution " + str(a))
        resolution = int(a)
    else:
        assert False, "unhandled option" + str(o) + " " + str(a)

nranks = servermanager.vtkProcessModule.GetProcessModule().GetNumberOfLocalPartitions()
rank   = servermanager.vtkProcessModule.GetProcessModule().GetPartitionId()
print ('Number of processes ' + str(nranks) + ' this is rank ' + str(rank))

paraview.simple._DisableFirstRenderCameraReset()


##############
# Pipeline
##############

numcells = resolution*resolution*resolution*nranks
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

filter2 = tetrahedralize1

# create a new 'Particle Partition Filter'
if (filter==0):
  partitionFilter1 = MeshPartitionFilter(Input=filter2)
  partitionFilter1.PointWeightsArrayName = ''
  partitionFilter1.KeepInversePointLists = 0
  partitionFilter1.MaxAspectRatio = 5
  partitionFilter1.GhostHaloSize = ghostoverlap
else:
  partitionFilter1 = D3(Input=filter2)
  partitionFilter1.MinimalMemory = 1
  partitionFilter1.BoundaryMode = 'Assign cells uniquely'
  
filter3 = partitionFilter1

# create a new 'Shrink'
shrink1 = Shrink(Input=filter3)
# Properties modified on shrink1
shrink1.ShrinkFactor = 0.9

filter4 = shrink1

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
renderView1.ViewSize = [1024, 1024]

# show data in view
shrink1Display = Show(filter4, renderView1)
# trace defaults for the display properties.
shrink1Display.ColorArrayName = ['POINTS', 'RTData']
shrink1Display.LookupTable = rTDataLUT

# reset view to fit data
#renderView1.ResetCamera()

# show color bar/color legend
#shrink1Display.SetScalarBarVisibility(renderView1, True)

#### saving camera placements for all active views

# current camera placement for renderView1
distance = sidelength*math.sqrt(math.sqrt(2.0)*10.0)
renderView1.CameraPosition = [distance, distance, distance]
renderView1.CameraViewUp = [0.7306437409080044, 0.4132080413165191, -0.5435244598574407]
renderView1.CameraParallelScale = 121.11167448863597

#### uncomment the following to render all views
#RenderAllViews()
#alternatively, if you want to write images, you can use 
SaveScreenshot("wavelet-filter-"+str(filter)+".png")
#
stats.dump_stats()
