"""
Main file running the Fastflood model on the python side
"""

import dagger as dag
import numpy as np
import matplotlib.pyplot as plt
from scabbard import io
from scabbard import enumeration as en
import math
import warnings

class FastFlood(object):
	"""docstring for FastFlood"""
	def __init__(self):
		# To provide multiple constructors, I moved all the constructor to dedicated functions at the end of the page
		# One may want to load a file, start with a numpy array, generate random topo, ...
		pass



	def set_precipitations(self, precipitation):
		if(isinstance(precipitation, np.ndarray) == False):
			self.precipitation = np.full_like(self.topography,precipitation)
		else:
			self.precipitation = precipitation.ravel()
		self._.set_Qbase(self.precipitation * self.cellarea)



	def run(self, dt = 1e-3, method = en.Flooder.SFD_STATIC, experimental = False, **kwargs):

		if(isinstance(dt,np.ndarray)):
			self._.spatial_dt(dt.ravel())
		else:
			self._.set_dt(dt)


		if(self.monitor_dhw):
			self.dhw = self.get_hw()

		if(method == en.Flooder.SFD_STATIC):
			self._.run_SFD()

		elif(method == en.Flooder.MFD_STATIC):
			if(experimental == False):
				self._.run_MFD()
			else:
				self._.run_MFD_exp( kwargs["N_MFD_NOGRAPH"])


		elif method == en.Flooder.MFD_DYNAMIC:
			self._.run_MFD_dynamic()

		elif method == en.Flooder.CAESAR_LS:
			self._.caesar_lisflood()
		elif method == en.Flooder.CAESAR_LS_OMP:
			self._.caesar_lisflood_OMP()

		if(self.monitor_dhw):
			self.dhw = self.get_hw() - self.dhw 



	def get_hw(self):
		
		return self._.get_hw().reshape(self.rshp)

	def get_dhw(self):
		if(self.monitor_dhw == False):
			raise AttributeError("Cannot access to dhw if the option monitor_dhw is not active. ")
		return self.dhw.reshape(self.rshp)

	def get_Qwin(self):
		return self._.get_Qwin().reshape(self.rshp)

	def get_Qwout(self):
		return self._.get_Qwout().reshape(self.rshp)

	def get_HS(self):
		if(self.topology == en.Topology.D8):
			return dag.hillshade(self.connector, self._.get_topography()).reshape(self.rshp)
		else:
			warnings.warn("D4 hillshading not available, backporting to d8 (can take time)")
			tdag = dag.D8N(self.nx, self.ny, self.dx, self.dy, self.x_min, self.y_min)
			return dag.hillshade(tdag,self._.get_topography()).reshape(self.rshp)

	def get_topography(self):
		return self._.get_topography().reshape(self.rshp)





















def FF_from_file(file_name, precipitation = 1e-4, topology = en.Topology.D8):

	fastfloodobject = FastFlood()
		
	# Loading DEM data with rasterio
	dem = io.load_raster(file_name)
	
	# Storing topography
	fastfloodobject.topography = dem['array'].ravel()
	# model dimensions
	fastfloodobject.nx = dem['nx']
	fastfloodobject.ny = dem['ny']
	fastfloodobject.nxy = fastfloodobject.nx * fastfloodobject.ny
	fastfloodobject.dx = dem['dx']
	fastfloodobject.dy = dem['dy']
	fastfloodobject.cellarea = fastfloodobject.dx * fastfloodobject.dy
	fastfloodobject.dxy = math.sqrt(fastfloodobject.dx**2 + fastfloodobject.dy**2)
	fastfloodobject.x_min = dem["x_min"]
	fastfloodobject.x_max = dem["x_max"]
	fastfloodobject.y_min = dem["y_min"]
	fastfloodobject.y_max = dem["y_max"]
	fastfloodobject.extent = [fastfloodobject.x_min, fastfloodobject.x_max, fastfloodobject.y_min, fastfloodobject.y_max]

	# The connector
	fastfloodobject.connector = dag.D8N(dem["nx"], dem["ny"], dem["dx"], dem["dy"], dem["x_min"], dem["y_min"])
	fastfloodobject.DAG = dag.graph(fastfloodobject.connector)
	# fastfloodobject.DAG.init_graph(fastfloodobject.connector)
	fastfloodobject.DAG.set_LMR_method(dag.LMR.cordonnier_fill)
	fastfloodobject.rshp = (fastfloodobject.ny,fastfloodobject.nx)

	if(isinstance(precipitation, np.ndarray) == False):
		precipitation = np.full_like(fastfloodobject.topography,precipitation)

	fastfloodobject._ = dag.FF(fastfloodobject.DAG, fastfloodobject.connector, fastfloodobject.topography, precipitation.ravel())
	fastfloodobject.DAG.compute_graph(fastfloodobject.topography, True, True)
	fastfloodobject.topology = topology

	if(topology == en.Topology.D8):
		fastfloodobject._.set_topological_number(0.5)
	elif(topology == en.Topology.D4):
		fastfloodobject._.set_topological_number(1)


	# Monitoring options
	fastfloodobject.monitor_dhw = False
	fastfloodobject.dhw = np.zeros_like(fastfloodobject.topography)

	return fastfloodobject




def FF_from_array(nx,ny,dx,dy,xmin,ymin, array, topology = en.Topology.D8, precipitation = 1e-4, BCs = None):

	fastfloodobject = FastFlood()
		
	# Storing topography
	fastfloodobject.topography = array.ravel()
	# model dimensions
	fastfloodobject.nx = nx
	fastfloodobject.ny = ny
	fastfloodobject.nxy = fastfloodobject.nx * fastfloodobject.ny
	fastfloodobject.dx = dx
	fastfloodobject.dy = dy
	fastfloodobject.cellarea = fastfloodobject.dx * fastfloodobject.dy
	fastfloodobject.dxy = math.sqrt(fastfloodobject.dx**2 + fastfloodobject.dy**2)
	fastfloodobject.x_min = xmin
	fastfloodobject.x_max = xmin + (nx + 1) * dx
	fastfloodobject.y_min = ymin
	fastfloodobject.y_max = ymin + (ny + 1) * dy
	fastfloodobject.extent = [fastfloodobject.x_min, fastfloodobject.x_max, fastfloodobject.y_min, fastfloodobject.y_max]

	# The connector
	fastfloodobject.connector = dag.D8N(fastfloodobject.nx, fastfloodobject.ny, fastfloodobject.dx, fastfloodobject.dy, fastfloodobject.x_min, fastfloodobject.y_min) 
	
	if(BCs is not None):
		fastfloodobject.connector.set_custom_boundaries(BCs)

	# The graph
	fastfloodobject.DAG = dag.graph(fastfloodobject.connector)
	# fastfloodobject.DAG.init_graph(fastfloodobject.connector)
	fastfloodobject.DAG.set_LMR_method(dag.LMR.cordonnier_fill)
	fastfloodobject.rshp = (fastfloodobject.ny,fastfloodobject.nx)

	if(isinstance(precipitation, np.ndarray) == False):
		precipitation = np.full_like(fastfloodobject.topography,precipitation)

	fastfloodobject._ = dag.FF(fastfloodobject.DAG, fastfloodobject.connector, fastfloodobject.topography, precipitation.ravel())
	fastfloodobject.DAG.compute_graph(fastfloodobject.topography, True, True)
	fastfloodobject.topology = topology

	if(topology == en.Topology.D8):
		fastfloodobject._.set_topological_number(0.5)
	elif(topology == en.Topology.D4):
		fastfloodobject._.set_topological_number(1)


	# Monitoring options
	fastfloodobject.monitor_dhw = False
	fastfloodobject.dhw = np.zeros_like(fastfloodobject.topography)

	return fastfloodobject









		#end of class

