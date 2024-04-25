'''
Contains generic widgets customisation for matplotlib qt backend
For example specific types of sliders (float, min max)
'''

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

# from PySide6 import QtWidgets, QtCore
from matplotlib.backends.qt_compat import QtWidgets, QtCore
import matplotlib

import scabbard as scb

# from scabbard.steenbok.gui_colormap_picker import str2cmap, ColorMapWidget

# Ensure using the Qt6Agg backend with PySide6
matplotlib.use('QtAgg')






class MapWidget(QtWidgets.QWidget):
	'''
	MapWidget, plots a topographic map from an environment
	'''

	
	def __init__(self, env):
		super(MapWidget, self).__init__()

		# this widget is made to be connected to an environment, for updates
		self.env = env

		# grid layout is the easiest layout to organise everything by row and col
		self.grid_layout = QtWidgets.QGridLayout(self)

		# Create two matplotlib figures and set their canvases
		self.figure, self.ax1 = plt.subplots()
		self.canvas = FigureCanvas(self.figure)


		# Add canvases and toolbars to the layout
		self.grid_layout.addWidget(NavigationToolbar(self.canvas, self), 0, 0)  # Row 1, Column 0
		self.grid_layout.addWidget(self.canvas, 1, 0)  # Row 0, Column 0

		# Initial plots
		self.imTopo = self.ax1.imshow(self.env.grid.Z2D, cmap='gist_earth')
		self.imHS = self.ax1.imshow(self.env.grid.hillshade, cmap='gray', alpha = 0.6)

		self.drapePlot = None

		self.ax1.set_xlabel("X (m)")
		self.ax1.set_ylabel("Y (m)")

		# Controls layout for the first plot
		self.controls_layout = QtWidgets.QHBoxLayout()
		self.grid_layout.addLayout(self.controls_layout, 2, 0)  # Row 2, Column 0



		# Draw the canvases
		self.canvas.draw()
		# self.canvas2.draw()




def map_widget_from_fname(fname):
	'''
	Loads a map widget from a file
	'''
	
	env = scb.env_from_DEM(fname)

	return MapWidget(env)

def map_widget_from_fname_for_graphflood(fname):
	'''
	Loads a map widget from a file
	'''
	
	env = scb.env_from_DEM(fname)
	mapw = MapWidget(env)
	mapw.imHS.set_alpha(1.)
	mapw.drapePlot = mapw.ax1.imshow(np.zeros_like(env.grid.Z2D), cmap='Blues', alpha = 0.6, vmin = 0, vmax = 0.8)
	plt.colorbar(mapw.drapePlot, label = 'Flow depth (m)')

	return mapw