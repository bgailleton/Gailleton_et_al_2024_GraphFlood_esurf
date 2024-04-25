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
# Ensure using the Qt6Agg backend with PySide6
matplotlib.use('QtAgg')



class MainWindow(QtWidgets.QMainWindow):

	
	def __init__(self, title = "Main Window"):
		super(MainWindow, self).__init__()

		self.setWindowTitle(title)
		
		# Get the primary screen's geometry
		screen_geometry = QtWidgets.QApplication.primaryScreen().geometry()
		width, height = screen_geometry.width() * 0.9, screen_geometry.height() * 0.9

		# Calculate the position to center the window, then move it slightly down
		x = screen_geometry.width() * 0.05
		y = (screen_geometry.height() - height) / 2 + screen_geometry.height() * 0.02

		# Set the window geometry
		self.setGeometry(int(x), int(y), int(width), int(height))
		# grid layout is the easiest layout to organise everything by row and col
		self.centralWidget = QtWidgets.QWidget()
		self.grid_layout = QtWidgets.QGridLayout(self.centralWidget)
		# self.centralWidget.setLayout(self.grid_layout)
		self.setCentralWidget(self.centralWidget)
		self.widgets = {}
		self.layouts = {}
		self.layouts['main_grid'] = self.grid_layout


	def add_widget(self, ref, widget, row = None, col = None, rowmax = None, colmax = None, replace_ref = True, parentLayout = None):

		if ref in self.widgets.keys() and replace_ref:
			self.widgets[ref].deleteLater()
		# widget.setParentLayout(self)
		
		self.widgets[ref] = widget

		adada = self.layouts[parentLayout] if parentLayout else self.grid_layout


		if(rowmax):
			adada.addWidget(widget, row, col,rowmax, colmax)
		elif(row):
			adada.addWidget(widget, row, col)
		else:
			adada.addWidget(widget)


		adada.addStretch(1)

	def add_layout(self, ref, layout, row, col,rowmax = None, colmax = None, replace_ref = True, parentLayout = None, style = None):

		if ref in self.layouts.keys() and replace_ref:
			self.layouts[ref].deleteLater()
		# widget.setparentLayout(self)
		
		self.layouts[ref] = layout

		adada = self.layouts[parentLayout] if parentLayout else self.grid_layout

		twidget = QtWidgets.QWidget()
		twidget.setLayout(layout)
		if(style):
			twidget.setStyleSheet(style)


		if(rowmax):
			adada.addWidget(twidget, row, col,rowmax, colmax)
		else:
			adada.addWidget(twidget, row, col)






# end of file