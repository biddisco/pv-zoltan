#### import the simple module from the paraview
from paraview.simple import *
import socket, os, sys, re

# ------------------------------
# setup machine dependent things
# ------------------------------
def load_plugins():
  hostname = socket.gethostname()
  plugins = []
  if "carona" in hostname or "public" in hostname:
    print(hostname + " contains Carona, laptop usage")
    data_path = '/Users/biddisco/data/sphflow/0100millions'
  #  data_path = '/Users/biddisco/data/sphflow/0001millions/hdf5'
    output_path = '/Users/biddisco/data/sphflow/resampled'
    image_path = '/Users/biddisco/data/sphflow/images'
    plugins.append('/Users/biddisco/build/egpgv/bin/libpv_zoltan.dylib')
    plugins.append('/Users/biddisco/build/egpgv/bin/libpv_meshless.dylib')
  elif "daint" in hostname:
    print("Running on some other machine")
    data_path = '/scratch/daint/biddisco/data/sphflow/big'
    output_path = '/scratch/daint/biddisco/data/sphflow/resampled'
    image_path = '/scratch/daint/biddisco/data/sphflow'
    plugins.append('/scratch/daint/biddisco/egpgv/libpv_zoltan.so')
    plugins.append('/scratch/daint/biddisco/egpgv/libpv_meshless.so')
  else:
    print("Running on some other machine - using daint settings ", hostname)
    data_path = '/scratch/daint/biddisco/data/sphflow/big'
    output_path = '/scratch/daint/biddisco/data/sphflow/resampled'
    image_path = '/scratch/daint/biddisco/data/sphflow'
    plugins.append('/scratch/daint/biddisco/egpgv/libpv_zoltan.so')
    plugins.append('/scratch/daint/biddisco/egpgv/libpv_meshless.so')

  ### Load pv-zoltan plugin
  for p in plugins:
    print("Loading plugin " + p)
    paraview.servermanager.LoadPlugin(p)

  return [data_path,output_path, image_path]
