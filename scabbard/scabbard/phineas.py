# Plotting wizard
import click
import scabbard as scb
import matplotlib.pyplot as plt
import dagger as dag

# ANSI escape sequences for colors
RESET = "\x1b[0m"
RED = "\x1b[31m"
GREEN = "\x1b[32m"
YELLOW = "\x1b[33m"
BLUE = "\x1b[34m"
MAGENTA = "\x1b[35m"
CYAN = "\x1b[36m"


@click.command()
@click.argument('fname', type = str)
def simplemapwizard(fname):
	plt.ioff()
	dem = scb.raster2RGrid(fname)
	atlas = scb.Dplot.basemap(dem)
	atlas.fig.show()
	plt.pause(0.01)
	input("press Enter to continue")










	


@click.command()
@click.option('-c', '--courant', 'courant',  type = float, default = 1e-3)
@click.option('-d', '--dt', 'dt',  type = float, default = None)
@click.option('-P', '--precipitations', 'precipitations',  type = float, default = 30)
@click.option('-m', '--manning', 'manning',  type = float, default = 0.033)
@click.option('-S', '--SFD', 'SFD', type = bool, default=False)
@click.option('-U', '--update_step', 'nupdate', type = int, default=10)
@click.option('-exp', '--experimental', 'experimental', type = bool, default=False, is_flag = True)
# @click.option('-h', '--help', 'help', type = bool, is_flag = True)
@click.argument('fname', type = str)
def graphflood_basic(fname,courant,dt,precipitations,manning,SFD,nupdate, experimental):
	print("EXPERIMENTAL IS ", experimental)

	P = precipitations /1e3/3600
	if(dt is not None and courant == 1e-3):
		print(RED + "WARNING, dt set to constant value AND courant, is that normal?\n\n" + RESET)

	print("+=+=+=+=+=+=+=+=+=+=+=+=+")
	print("+=+=+=GRAPHFLOOD+=+=+=+=+")
	print("+=+=+=+=+=+=+=+=+=+=+=+=+\n\n")

	print("Precipitation rates = ", P, " m/s (",precipitations," mm/h)")


	mod = scb.ModelHelper()
	mod.init_dem_model(fname, sea_level = 0., P = precipitations)

	# Number of nodes in the Y direction
	ny = mod.grid.ny
	nx = mod.grid.nx


	mod.courant = False if(dt is not None) else True
	mod.stationary = True 
	# mod.gf.flood.enable_Qwout_recording()

	# manning friction:
	mod.mannings = manning

	# Single flow solver?
	mod.SFD = SFD;


	mod.dt = dt if dt is not None else 1e-3
	mod.min_courant_dt = 1e-6
	mod.courant_number = courant

	ph = scb.PlotHelper(mod)
	ph.init_hw_plot(use_extent = False)

	update_fig = nupdate
	i = 0
	j = 0
	while True:
		i+=1
		mod.run() if(experimental == False) else mod.gf.flood.run_hydro_only()
		if(i % update_fig > 0):
			continue
		hw = mod.gf.flood.get_hw().reshape(mod.gf.grid.rshp)
		print("Running step", i)
		ph.update()




@click.command()
@click.argument('fname', type = str)
def _debug_1(fname):
	plt.ioff()
	dem = scb.raster2RGrid(fname)
	atlas = scb.Dplot.basemap(dem)
	atlas.fig.show()
	plt.pause(0.01)
	while(True):
		plt.pause(1)
		dem.add_random_noise(-10,10)
		atlas.update()
	input("press Enter to continue")



@click.command()
def haguid():
	from scabbard.haguid import launch_haGUId

	launch_haGUId()



@click.command()
@click.argument('fname', type = str)
def visu2Dnpy(fname):

	import sys
	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.backends.backend_qtagg import FigureCanvas
	from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

	# from PySide6 import QtWidgets, QtCore
	from matplotlib.backends.qt_compat import QtWidgets, QtCore
	import matplotlib
	from matplotlib.lines import Line2D

	import scabbard.phineas_helper as phelp

	# Ensure using the Qt6Agg backend with PySide6
	matplotlib.use('QtAgg')

	class MatplotlibWidget(QtWidgets.QWidget):
	
		def __init__(self, data):
			super(MatplotlibWidget, self).__init__()

			self.data = data

			self.cid = None
			self.cid2 = None

			self.initUI()
			self.setWindowTitle("Scabbard numpy 2D array explorer")
			 # Get the primary screen's geometry
			screen_geometry = QtWidgets.QApplication.primaryScreen().geometry()
			width, height = screen_geometry.width() * 0.8, screen_geometry.height() * 0.8

			# Calculate the position to center the window, then move it slightly down
			x = screen_geometry.width() * 0.1
			y = (screen_geometry.height() - height) / 2 + screen_geometry.height() * 0.05

			# Set the window geometry
			self.setGeometry(int(x), int(y), int(width), int(height))

		def initUI(self):
			# Main layout
			grid_layout = QtWidgets.QGridLayout(self)

			# Create two matplotlib figures and set their canvases
			self.figure1, self.ax1 = plt.subplots()
			self.canvas1 = FigureCanvas(self.figure1)
			self.figure2, self.ax2 = plt.subplots()
			self.canvas2 = FigureCanvas(self.figure2)

			# Add canvases and toolbars to the layout
			grid_layout.addWidget(self.canvas1, 1, 0)  # Row 0, Column 0
			grid_layout.addWidget(NavigationToolbar(self.canvas1, self), 0, 0)  # Row 1, Column 0
			grid_layout.addWidget(self.canvas2, 1, 1)  # Row 0, Column 1
			grid_layout.addWidget(NavigationToolbar(self.canvas2, self), 0, 1)  # Row 1, Column 1

			# Initial plots
			self.imshow1 = self.ax1.imshow(self.data, cmap='magma')
			self.colorbar1 = self.figure1.colorbar(self.imshow1, ax=self.ax1)
			self.plot1 = Line2D([0,self.data.shape[1]], [0, 0], color='r', linestyle='--')
			self.ax1.add_line(self.plot1)
			self.plot2 = self.ax2.plot(self.data[0,:], color = 'k', lw = 4)

			# Controls layout for the first plot
			controls_layout1 = QtWidgets.QHBoxLayout()
			grid_layout.addLayout(controls_layout1, 2, 0)  # Row 2, Column 0

			


			# Colormap selection for the first plot
			self.colormapComboBox1 = QtWidgets.QComboBox()
			self.colormapComboBox1.addItems(['magma', 'viridis', 'cividis', 'Blues', 'Reds', 'RdBu_r'])
			self.colormapComboBox1.currentIndexChanged.connect(lambda: self.update_plot('1'))
			controls_layout1.addWidget(self.colormapComboBox1)

			#RangeSelection
			self.rangeSlider1 = phelp.RangeSlider(self.data.min(), self.data.max(), ID = '1')
			# grid_layout.addWidget(self.rangeSlider1,3,0)
			self.rangeSlider1.rangeChanged.connect(self._update_crange)
			controls_layout1.addWidget(self.rangeSlider1)



			# Controls layout for the second plot
			controls_layout2 = QtWidgets.QHBoxLayout()
			grid_layout.addLayout(controls_layout2, 2, 1)  # Row 2, Column 1

			# Colormap selection for the second plot
			# self.colormapComboBox2 = QtWidgets.QComboBox()
			# self.colormapComboBox2.addItems(['magma', 'viridis', 'cividis', 'Blues', 'Reds', 'RdBu_r'])
			# self.colormapComboBox2.currentIndexChanged.connect(lambda: self.update_plot('2'))
			# controls_layout2.addWidget(self.colormapComboBox2)

			self.button = QtWidgets.QPushButton("Cross section",self)
			self.button.clicked.connect(self.switch_coordinate_picking)
			controls_layout2.addWidget(self.button)

			# Draw the canvases
			self.canvas1.draw()
			self.canvas2.draw()

		def update_plot(self, plt_ID):
			cmap1 = self.colormapComboBox1.currentText()
			# cmap2 = self.colormapComboBox2.currentText()

			if plt_ID == '1':
				self.imshow1.set_cmap(cmap1)


			# Redraw the canvases
			self.canvas1.draw()

		def _update_crange(self, vmin,vmax, plt_ID):
			if plt_ID == '1':
				self.imshow1.set_clim(vmin,vmax)


			# Redraw the canvases
			self.canvas1.draw()

		def switch_coordinate_picking(self):
			# Disconnect existing connections if any
			if self.cid:
				self.canvas1.mpl_disconnect(self.cid)
				self.cid = None
			else:
				# Connect the mouse click event to the on_motion function
				self.cid = self.canvas1.mpl_connect('motion_notify_event', self.on_motion)

		def plot_CS_from_coord(self, tx, ty):
			row = round(ty) #will not be true for maps, is just for test here
			self.plot2[0].set_ydata(self.data[row,:])
			self.ax2.relim()
			self.ax2.autoscale_view(True, True, True)
			self.plot1.set_ydata([row,row])
			self.canvas1.draw()
			self.canvas2.draw()
			# self.switch_coordinate_picking()


		def on_motion(self, event):
			# Call the function with the coordinates
			if event.xdata is not None and event.ydata is not None:
				self.cid2 = self.canvas1.mpl_connect('button_press_event', self.on_click)
				self.plot_CS_from_coord(event.xdata, event.ydata)

		def on_click(self, event):
			# Call the function with the coordinates
			if event.xdata is not None and event.ydata is not None:
				self.switch_coordinate_picking()
				if self.cid2:
					self.canvas1.mpl_disconnect(self.cid2)
					self.cid2 = None



	app = QtWidgets.QApplication(sys.argv)
	data = np.load(fname)
	if(len(data.shape) != 2):
		raise ValueError("Needs to be 2D numpy array")
	main = MatplotlibWidget(data)
	main.show()
	sys.exit(app.exec_())
