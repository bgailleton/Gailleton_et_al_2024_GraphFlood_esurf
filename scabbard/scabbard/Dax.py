import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cmcrameri as cm
import string
import random


class Dax(object):
	
	"""
		docstring for Dax
	"""

	def __init__(self, ax, key = None, zorder = 1, fig = None):
		
		super(Dax, self).__init__()

		self.fig = fig
		
		self.ax = ax

		self.zorder = zorder
		
		# // Generating random key if none is given	
		if(key is not None):
			self.key = key
		else:
			self.key = ''.join(random.choices(string.ascii_letters + string.digits, k=8))


	def __hash__(self):
		return hash(self.key)

	def update(self):
		pass


class callbax_image(object):
	
	def __init__(self, im, callback, rshp, clim = None, callback_params = None):
		self.im = im
		self.callback = callback
		self.callback_params = callback_params
		self.rshp = rshp
		self.clim = clim

	def update(self):
		if self.callback_params is None :
			arr = self.callback().reshape(self.rshp)
		else:
			arr = self.callback(*self.callback_params).reshape(self.rshp)
		self.im.set_data(arr)
		if(self.clim is None):
			self.im.set_clim(np.nanmin(arr), np.nanmax(arr))

class callbax_sline(object):
	
	def __init__(self, ax, sline, callback, axylim = None, axylim_ignore = False):
		self.ax = ax
		self.sline = sline[0]
		self.callback = callback
		self.axylim = axylim
		self.axylim_ignore = axylim_ignore

	def update(self):
		
		arr = self.callback()
		self.sline.set_xdata(arr[0])
		self.sline.set_ydata(arr[1])

		if(self.axylim is None and self.axylim_ignore == False):
			perc = 0.05 * (np.nanmax(arr[1]) - np.nanmin(arr[1]))
			self.ax.set_ylim(np.nanmin(arr[1]) - perc, np.nanmax(arr[1]) + perc)
			perc = 0.05 * (np.nanmax(arr[0]) - np.nanmin(arr[0]))
			self.ax.set_xlim(np.nanmin(arr[0]) - perc, np.nanmax(arr[0]) + perc)




class RGridDax(Dax):
	"""docstring for RGridDax"""
	def __init__(self,Rgrid, ax, key = None, cmap = "gist_earth", hillshade = True, alpha_hillshade = 0.45, clim = None, callback_topo = None):
		
		super(RGridDax, self).__init__(ax, key)
		self.grid = Rgrid
		self.base_ax = self.ax.imshow(self.grid.Z2D, extent = self.grid.extent(), cmap = cmap, zorder =  self.zorder)
		self.hillshade_on = hillshade

		self.clim = clim if clim is not None else [self.grid.min(), self.grid.max()]

		if(hillshade):
			self.hillshade_ax = self.ax.imshow(self.grid.hillshade, extent = self.grid.extent(), cmap = "gray", alpha = alpha_hillshade, zorder = self.zorder)

		self.ax.set_xlabel("X (m)")
		self.ax.set_ylabel("Y (m)")

	def update(self):
		'''
		'''
		self.base_ax.set_data(self.grid.Z2D)
		if(self.hillshade_on):
			self.hillshade_ax.set_data(self.grid.hillshade)


	def drape_on(self, array, cmap = "gray", clim = None, delta_zorder = 1, alpha = 0.5, callback = None, callback_params = None):

		if clim is not None:
			vmin = clim[0]
			vmax = clim[1]
		else:
			vmin = np.nanmin(array)
			vmax = np.nanmax(array)


		ttaaxx = self.ax.imshow(array.reshape(self.grid.rshp), extent = self.grid.extent(), cmap = cmap, alpha = alpha, zorder = self.zorder + delta_zorder, vmin = vmin, vmax = vmax)

		if(callback == None):
			def nonecback():
				pass
			callback = nonecback


		return callbax_image(ttaaxx, callback, self.grid.rshp, clim = clim, callback_params = callback_params)

























































# dskjhksadfhla