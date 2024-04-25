'''
This is a first attempt to make a map automator.
Probably will just test different way here before making a better one
B.G
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scabbard import Dfig, Dax, RGridDax



def basemap(Rgrid, **kwargs):

	fig, ax = plt.subplots()
	if("figsize" in kwargs.keys()):
		fig.set_size_inches(*figsize)

	tDax = RGridDax(Rgrid, ax)
	tdfig = Dfig(fig, tDax)
	return tdfig






























#end of files
