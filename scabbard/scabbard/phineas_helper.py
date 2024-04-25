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



class FloatSlider(QtWidgets.QSlider):
	'''
	A widget that emulate a floating point slider
	'''

	def __init__(self, orientation, parent=None):
		super(FloatSlider, self).__init__(orientation, parent)

		self._multiplier = 100  # Define your own multiplier for float precision
		self._step = 0.5  # Define the step value you want for the slider

	def setFloatRange(self, min_val, max_val):
		self.setMinimum(int(min_val * self._multiplier))
		self.setMaximum(int(max_val * self._multiplier))

	def setFloatValue(self, value):
		intValue = int(value * self._multiplier)
		self.setValue(intValue)

	def floatValue(self):
		return self.value() / self._multiplier

	def setFloatStep(self, step):
		self._step = step
		self.setSingleStep(int(step * self._multiplier))


class RangeSlider(QtWidgets.QWidget):

	rangeChanged = QtCore.Signal((float, float, str))

	def __init__(self, minimum, maximum, parent=None, ID = ""):
		super(RangeSlider, self).__init__(parent)
		self.layout = QtWidgets.QVBoxLayout(self)

		self.minimum = minimum
		self.maximum = maximum



		self.ID = str(ID)

		self.minSlider = FloatSlider(QtCore.Qt.Horizontal, self)
		self.maxSlider = FloatSlider(QtCore.Qt.Horizontal, self)



		self.minSlider.setFloatRange(minimum,maximum)
		self.maxSlider.setFloatRange(minimum,maximum)

		self.minSlider.setFloatValue(minimum)
		self.maxSlider.setFloatValue(maximum)

		# self.minSlider.setMinimum(minimum)
		# self.minSlider.setMaximum(maximum)
		# self.maxSlider.setMinimum(minimum)
		# self.maxSlider.setMaximum(maximum)

		self.layout.addWidget(QtWidgets.QLabel("Min:"))
		self.layout.addWidget(self.minSlider)
		self.layout.addWidget(QtWidgets.QLabel("Max:"))
		self.layout.addWidget(self.maxSlider)

		self.minValue = minimum
		self.maxValue = maximum

		self.minSlider.valueChanged.connect(self.updateRange)
		self.maxSlider.valueChanged.connect(self.updateRange)

		

	def updateRange(self):
		minVal = self.minSlider.floatValue()
		maxVal = self.maxSlider.floatValue()

		if minVal > maxVal:
			minVal, maxVal = maxVal, minVal

		self.minValue = minVal
		self.maxValue = maxVal

		# Update the plot or imshow limits here
		# Make sure the sliders stay within each other's bounds
		# self.minSlider.setMaximum(self.maxValue)
		# self.maxSlider.setMinimum(self.minValue)
		self.minSlider.setFloatValue(minVal)
		self.maxSlider.setFloatValue(maxVal)

		self.rangeChanged.emit(self.minValue, self.maxValue, self.ID)