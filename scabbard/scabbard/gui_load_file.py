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






class FileLoader(QtWidgets.QWidget):

	theChosenFile = QtCore.Signal((str))
	
	def __init__(self):
		super(FileLoader, self).__init__()
		self.layout = QtWidgets.QVBoxLayout(self)

		# Button to open the file dialog
		self.openButton = QtWidgets.QPushButton("Open File")
		self.openButton.clicked.connect(self.openFileDialog)
		self.layout.addWidget(self.openButton)
		self.setMaximumSize(100,100)

	def openFileDialog(self):
		# Open a file dialog and print the selected file path
		filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open File", "", "All Files (*)")
		self.theChosenFile.emit(filename)
		if filename:
			print(f"Selected file: {filename}")