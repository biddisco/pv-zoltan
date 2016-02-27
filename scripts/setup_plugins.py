#### import the simple module from the paraview
from paraview.simple import *
import socket, os, sys, re

# ------------------------------
# setup machine dependent things
# ------------------------------
def load_plugins():
  hostname = socket.gethostname()
  plugins = []
  if "carona" in hostname:
    print(hostname + " contains Carona, laptop usage")
    plugins.append('/Users/biddisco/build/egpgv/bin/libpv_zoltan.dylib')
    plugins.append('/Users/biddisco/build/egpgv/bin/libpv_meshless.dylib')
  #  data_path = '/Users/biddisco/data/sphflow/0100millions'
    data_path = '/Users/biddisco/data/sphflow/0001millions/hdf5'
  elif "daint" in hostname:
    print("Running on some other machine")
    data_path = '/scratch/daint/biddisco/data/sphflow/medium'
    plugins.append('/scratch/daint/biddisco/egpgv/libpv_zoltan.so')
    plugins.append('/scratch/daint/biddisco/egpgv/libpv_meshless.so')
  else:
    # error, we must put some paths in here
    pass

  ### Load pv-zoltan plugin
  for p in plugins:
    print("Loading plugin " + p)
    paraview.servermanager.LoadPlugin(p)

  return data_path
