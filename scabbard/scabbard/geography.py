'''
Manages the geogrphical conversion and all
WILL EVOLVE QUITE A LOT
Not sure I know how I want ot handle that ATM
'''


class geog:
	'''
	Simple class wrapping the geographical information
	Will evolve
	'''

	def __init__(self, 
		xmin = 0., 
		ymin = 0.,
		xmax = 0.,
		ymax = 0.,
		crs = None

	):
		
		self.xmin = xmin
		self.ymin = ymin
		self.xmax = xmax
		self.ymax = ymax
		self.crs = crs


