"""
This class deals with loading raster informations
Authors: B.G.
"""
import numpy as np
import rasterio as rio
from rasterio.transform import from_bounds
import dagger as dag
from rasterio.transform import from_origin


def load_raster(fname):
	"""
	Load a raster array with different options. It uses rasterio that itself uses gdal.
	Arguments:
		fname (str): raster to load (path+file_name+format)
	Returns:
		A python dictionnary containing the following "key" -> val:
			"res" -> Resolution of the DEM
			"ncols" -> number of columns
			"nrows" -> number of rows
			"x_min" -> well x minimum
			"y_min" -> and y minimum
			"x_max" -> and x maximum
			"y_max" -> and x maximum
			"extent" -> extent combined in order to match matplotlib
			"array" -> numpy 2D array containing the data
			"crs" -> The crs string (geolocalisation)
			"nodata" -> list of nodata values
	Authors:
		B.G.
	Date:
		23/02/2019
	"""

	# Loading the raster with rasterio
	this_raster = rio.open(fname)

	# Initialising a dictionary containing the raster info for output
	out = {}
	# I am padding a no_data contour
	gt = this_raster.res
	out['dx'] = gt[0]
	out['dy'] = gt[1]
	out["nx"] = this_raster.width
	out["ny"] = this_raster.height
	out["x_min"] = this_raster.bounds[0]
	out["y_min"] = this_raster.bounds[1]
	out["x_max"] = this_raster.bounds[2]
	out["y_max"] = this_raster.bounds[3]
	corr = out['dx'] + out['dy']
	out["extent"] = [out["x_min"],out["x_max"]-corr,out["y_min"],out["y_max"]-corr]
	out["array"] = this_raster.read(1).astype(np.float64)
	try:
		out['crs'] = this_raster.crs['init']
	except (TypeError, KeyError) as e:
		out['crs'] = u'epsg:32601'
	out['nodata'] = this_raster.nodatavals

	
	
	# pixelSizeY =-gt[4]
	# (left=358485.0, bottom=4028985.0, right=590415.0, top=4265115.0)

	return out




def raster2graphcon(file_name):
	"""
	Ingest a raster with rasterio adn load a connector and graph from it
	"""
		
	# Loading DEM data with rasterio
	dem = load_raster(file_name)

	connector = dag.D8N(dem["nx"], dem["ny"], dem["dx"], dem["dy"], dem["x_min"], dem["y_min"])
	graph = dag.graph(connector)
	
	
	
	return connector, graph, dem

def raster2con(file_name):
	"""
	Ingest a raster with rasterio adn load a connector and graph from it
	"""
		
	# Loading DEM data with rasterio
	dem = load_raster(file_name)

	connector = dag.D8N(dem["nx"], dem["ny"], dem["dx"], dem["dy"], dem["x_min"], dem["y_min"])
	
	return connector, dem


# import rasterio
def save_raster(file_name, array, x_min, x_max, y_min, y_max, dx, dy, crs = None):
	# Define file name and file mode
	# file_name = "output.tif"
	file_mode = "w+"

	# Specify dimensions of the array
	height, width = array.shape
	count = 1  # Number of bands in the array
	dtype = array.dtype

	# Specify transformation parameters
	# x_min = 0  # X-coordinate of the top-left corner
	# y_max = 0  # Y-coordinate of the top-left corner
	x_res = dx  # Pixel size in the x-direction
	y_res = dy  # Pixel size in the y-direction
	transform = from_origin(x_min, y_max, dx, dy)

	if(crs is None):
		crs = "EPSG:35653"

	# Create output GeoTIFF file
	with rio.open(file_name, file_mode, driver='GTiff', height=height, width=width, dtype=dtype, count=count, transform=transform, nodata = -9999, crs = crs) as dst:
	    dst.write(array, 1)  # Write the array to the GeoTIFF file as the first band (band index is 1-based)
