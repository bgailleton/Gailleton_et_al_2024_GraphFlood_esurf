
def launch_haGUId():

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
	from scabbard.gui_main_window import MainWindow
	from scabbard.gui_load_file import FileLoader
	from scabbard.gui_float_sliders import FloatSlider
	import scabbard.gui_maps_widgets as guimap
	import taichi as ti
	import scabbard as scb
	import numpy as np
	import matplotlib.pyplot as plt

	import threading

	# Create an Event object
	stop_event = threading.Event()

	
	fnamedem = f"/home/bgailleton/Desktop/code/FastFlood2.0/FastFlood2_Boris/graphflood/paper_scripts/data/green_river_{1}.tif"
	
	manning = 0.033
	env = scb.env_from_DEM(fnamedem)
	nx,ny = env.grid.nx, env.grid.ny
	dx,dy = env.grid.dx, env.grid.dy
	dt = 2e-3

	hmaxplot = 3

	ti.init(arch=ti.gpu)  # Initialize Taichi to use the CPU


	Z = ti.field(dtype=ti.f32, shape=(ny,nx))
	Z.from_numpy(env.grid.Z2D.astype(np.float32))
	hw = ti.field(dtype=ti.f32, shape=(ny,nx))

	QwA = ti.field(dtype=ti.f32, shape=(ny,nx))
	QwB = ti.field(dtype=ti.f32, shape=(ny,nx))
	QwC = ti.field(dtype=ti.f32, shape=(ny,nx))

	P = ti.field(dtype = ti.f32, shape = ())
	P[None] = np.float32(2e-5)







	@ti.kernel
	def init_field():
		for i, j in Z:
			hw[i,j] = 0

	@ti.kernel
	def stepinit():
		for i,j in Z:
			QwA[i,j] += P[None] * dx * dy
			QwB[i,j] = P[None] * dx * dy

	@ti.func
	def Zw(i,j) -> ti.f32:
		return Z[i,j] + hw[i,j]

	@ti.func
	def Sw(i,j, ir,jr)->ti.f32:
		return (Zw(i,j) - Zw(ir,jr))/dx

	@ti.func
	def Qw(i,j, tSw:ti.f32)->ti.f32:
		return dx/manning * ti.math.pow(hw[i,j],(5./3.)) * ti.math.pow(tSw,0.5)

	@ti.func
	def neighbour(i,j,k:int):
		ir,jr = -1,-1
		valid = True
		if(i == 0):
			if(j == 0 and k <= 1):
				valid = False
			elif(j == nx-1 and (k == 0 or k == 2)):
				valid = False
			elif(k==0):
				valid = False
		elif(j == 0 and k == 1):
			valid = False
		elif(j == nx-1 and k == 2):
			valid = False
		elif(i == ny-1):
			if(j == 0 and (k == 1 or k == 3)):
				valid = False
			elif(j == nx-1 and (k == 3 or k == 2)):
				valid = False
			elif(k==3):
				valid = False

		if(valid):
			if(k == 0):
				ir,jr = i-1, j
			if(k == 1):
				ir,jr = i, j-1
			if(k == 2):
				ir,jr = i, j+1
			if(k == 3):
				ir,jr = i+1, j
		return ir, jr



	@ti.kernel
	def compute_Qw():
		for i,j in Z:

			if(i == ny-1 or i ==0 or j==0 or j == nx-1 or Z[i,j] < 0.):
				continue

			Sws = ti.math.vec4(0.,0.,0.,0.)
			sumSw = 0.
			SSx = 0.
			SSy = 0.
			lockcheck = 0
			while(sumSw == 0.):
				lockcheck += 1
				for k in range(4):
					ir,jr = neighbour(i,j,k)
					if(ir == -1):
						continue

					tS = Sw(i,j,ir,jr)
					if(tS <= 0):
						continue

					if(k == 0 or k == 3):
						if(tS > SSy):
							SSy = tS
					else:
						if(tS > SSx):
							SSx = tS

					Sws[k] = tS
					sumSw += tS
				if(sumSw == 0.):
					hw[i,j] += 1e-4
				if(lockcheck > 10000):
					break

			gradSw = ti.math.sqrt(SSx*SSx + SSy*SSy)
			if(gradSw == 0):
				continue
			QwC[i,j] = dx/manning * ti.math.pow(hw[i,j], 5./3) *sumSw/ti.math.sqrt(gradSw)
			# print(QwC[i,j])
			for k in range(4):
				ir,jr = neighbour(i,j,k)
				if(ir == -1):
					continue
				ti.atomic_add(QwB[ir,jr], Sws[k]/sumSw * QwA[i,j])



	@ti.kernel
	def compute_hw():
		for i,j in Z:
			QwA[i,j] = QwB[i,j]
			if(i == ny-1 or i ==0 or j==0 or j == nx-1 or Z[i,j] < 0.):
				continue
			hw[i,j] = ti.math.max(0.,hw[i,j] + (QwA[i,j] - QwC[i,j]) * dt/(dx*dy) ) 



		# GUI setup

	init_field()

	thread = None

	def run_model(stop_event, mainwindow):
		i = 0
		while not stop_event.is_set():
			i+=1
			stepinit()
			compute_Qw()
			compute_hw()

			if(i % 1000 == 0):
				mainwindow.widgets['mainmap'].drapePlot.set_data(hw.to_numpy())
				mainwindow.widgets['mainmap'].canvas.draw()

	def valueChanged(newP, stop_event, mainwindow):
		# changing to
		# print("changing to ", newP * 1e-3/3600)
		global thread
		stop_event.set()
		thread.join()
		stop_event.clear()
		P[None] = np.float32(newP * 1e-3/3600)
		# Create and start the thread
		thread = threading.Thread(target=run_model, args=(stop_event, mainwindow,))
		thread.start()







	def init_model(main_window, guimap, stop_event, fname):
		global thread
		mainwindow.add_widget("mainmap", guimap.map_widget_from_fname_for_graphflood(fname), parentLayout = 'map')
		Z.from_numpy(mainwindow.widgets['mainmap'].env.grid.Z2D)

		# Create and start the thread
		thread = threading.Thread(target=run_model, args=(stop_event, mainwindow,))
		thread.start()




	# Ensure using the Qt6Agg backend with PySide6
	matplotlib.use('QtAgg')

	# Creating the application (required by Qt)
	app = QtWidgets.QApplication(sys.argv)

	# the mainwindow is the main Widget containing the others
	mainwindow = MainWindow(title = "HaGUId")

	main_style = """
		QWidget {
			background-color: #D3D3D3; /* Light Gray */
			border-style: solid;
			border-width: 2px;
			border-color: black;
		}

		QLabel {
			max-height: 20px;
			font-family: 'Courier New';
			font-size: 14pt;
			font-weight: normal;
			color: black;
			text-align: left;
			border: none;
		}

		QPushButton {
			border: 2px solid #8f8f91;
			border-radius: 6px;
			background-color: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
											  stop: 0 #f6f7fa, stop: 1 #dadbde);
			min-width: 80px;
		}

		QPushButton:pressed {
			background-color: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
											  stop: 0 #dadbde, stop: 1 #f6f7fa);
		}

		QPushButton:flat {
			border: none; /* no border for a flat push button */
		}

		QPushButton:default {
			border-color: navy; /* make the default button prominent */
		}

		QSlider {
			min-width: 50 px;
			border: none;
		}
		"""

	from scabbard.gui_stylesheets import dark_orange
	main_style = dark_orange

	##############################################################################################################################
	##############################################################################################################################
	# Starting by placing placeholders and their respective styles
	mainwindow.add_layout('menu', QtWidgets.QHBoxLayout(), 0,0, style = main_style)

	mainwindow.add_layout('toolbox', QtWidgets.QHBoxLayout(), 0,1, style = main_style)

	mainwindow.add_layout('map', QtWidgets.QHBoxLayout(), 1,0, style = main_style)

	mainwindow.add_layout('toolopt', QtWidgets.QHBoxLayout(), 1,1, style = main_style)

	# ttoolbox = QtWidgets.QLabel("toolopt")
	# mainwindow.layouts['toolopt'].addWidget(ttoolbox)
	ttoolopt = QtWidgets.QLabel("Current options: set precipitation rates")
	mainwindow.layouts['toolbox'].addWidget(ttoolopt)

	# Adjusting the stretch factors to control the fractions
	# These factors do not correspond directly of percentage of area but affect the distribution of space
	mainwindow.grid_layout.setColumnStretch(0, 5) # 30% width approximately
	mainwindow.grid_layout.setColumnStretch(1, 3) # 20% width approximately
	mainwindow.grid_layout.setRowStretch(0, 1) # 30% height approximately
	mainwindow.grid_layout.setRowStretch(1, 5) # 20% height approximately
	##############################################################################################################################
	##############################################################################################################################

	def add_scroll_bars(mainwindow):
		mainwindow.add_widget( 'P_slider_lab',QtWidgets.QLabel("Precipitation rates (1-300 mm/yrs)"),parentLayout = 'toolopt')
		mainwindow.add_widget( 'P_slider',FloatSlider(QtCore.Qt.Vertical,vmin = 1,vmax = 300, step = 1, multiplier = 50, label = "Precipitation rates"),parentLayout = 'toolopt')
		mainwindow.widgets['P_slider'].valMod.connect(lambda x: valueChanged(x,stop_event, mainwindow))


	mainwindow.add_widget("loader", FileLoader(),  parentLayout = 'menu')
	mainwindow.widgets['loader'].theChosenFile.connect(lambda x: (init_model(mainwindow, guimap, stop_event,x), add_scroll_bars(mainwindow) ,None)[-1])



	mainwindow.show()
	sys.exit(app.exec_())





# stop_event.set()

# # Wait for the thread to finish
# thread.join()






























# end of file