from importlib import resources 
import numpy as np
import numba as nb
import scabbard as scb

import subprocess as sub

# Path to your JSON file
with resources.open_text('scabbard', '_blandplot.py') as config_file:
	blender_script = config_file.name



@nb.njit()
def _fill_arrays(data,displacement_array,color_array, intensity = 0.7):
    k=0
    for j in range(0,data.shape[1]+2):
        for i in range(0,data.shape[0]+2):
            i_data = i - 1
            j_data = j - 1

            i_minus = max(i_data-1,0)
            i_plus = min(i_data+1,data.shape[0]-1)
            j_minus = max(j_data-1,0)
            j_plus = min(j_data+1,data.shape[1]-1)
            
            i_data = min(max(i_data,0),data.shape[0]-1)
            j_data = min(max(j_data,0),data.shape[1]-1)
            
            if i >= 1 and j >= 1 and i <= data.shape[0] and j <= data.shape[1]:
                if np.isnan(data[i_data,j_data]):
                    displacement_array[k:k+4] = np.nan
                    color_array[k:k+4] = 1.0
                else:
                    displacement_array[k:k+4] = data[i_data,j_data]
                    # Get the RGBA color for the specified value
    #                rgba_color = np.array(sm.to_rgba(data[i_data,j_data]))
    #                rgba_color *= luminosity_scale
                    
                    color_array[k:k+4] = data[i_data,j_data] * intensity
                    # color_array[k:k+3] *= luminosity_scale
    #                color_array[k:k+4] = rgba_color
            else:
                displacement_array[k:k+4] = np.nan
                color_array[k:k+4] = 1.0
            k+=4
            
@nb.njit()
def _fill_arrays_cc(data, colors,displacement_array,color_array, intensity = 0.7):
    k=0
    for j in range(0,data.shape[1]+2):
        for i in range(0,data.shape[0]+2):
            i_data = i - 1
            j_data = j - 1

            i_minus = max(i_data-1,0)
            i_plus = min(i_data+1,data.shape[0]-1)
            j_minus = max(j_data-1,0)
            j_plus = min(j_data+1,data.shape[1]-1)
            
            i_data = min(max(i_data,0),data.shape[0]-1)
            j_data = min(max(j_data,0),data.shape[1]-1)
            
            if i >= 1 and j >= 1 and i <= data.shape[0] and j <= data.shape[1]:
                if np.isnan(data[i_data,j_data]):
                    displacement_array[k:k+4] = np.nan
                    color_array[k:k+4] = 1.0
                else:
                    displacement_array[k:k+4] = data[i_data,j_data]
                    # Get the RGBA color for the specified value
    #                rgba_color = np.array(sm.to_rgba(data[i_data,j_data]))
    #                rgba_color *= luminosity_scale
                    
                    color_array[k:k+4] = colors[i_data,j_data] 
                    # color_array[k:k+3] *= luminosity_scale
    #                color_array[k:k+4] = rgba_color
            else:
                displacement_array[k:k+4] = np.nan
                color_array[k:k+4] = 1.0
            k+=4
            



def plot_blender_from_array(_data, 
	fprefix = "__blend",
	save_prefix = "render",
	dx = 5.,
	gpu = True,
	perspective = True,
	ortho_scale = 2.,
	focal_length = 40,
	f_stop = 50000,
	shiftx = 0.,
	shifty = 0.,
	camera_tilt =  45.0,
	camera_rotation =  0.,
	sun_tilt =  25.0,
	sun_rotation =  315.0,
	sun_intensity =  0.2,
	exaggeration =  1.5,
	recast_minZ =  0.,
	recast_maxZ =  1.,
	number_of_subdivisions = 1000,
	res_x = 1920,
	res_y = 1080,
	samples = 50,
	intensity = 0.7,
	custom_colors = None

	):
	
	data = np.copy(_data)

	# create an image that holds the data
	displacement_array = np.zeros(((data.shape[0]+2)*(data.shape[1]+2)*4), dtype=np.float32)
	color_array = np.zeros(((data.shape[0]+2)*(data.shape[1]+2)*4), dtype=np.float32)

	# find minimum elevation
	data_min = np.min(data[data!=-9999])
	# find maximum elevation
	data_max = np.max(data[data!=-9999])
	# calculate total relief
	relief = data_max - data_min
	#scale values from 0 to 1
	data[data!=-9999] = (data[data!=-9999] - data_min) / relief
	# recasting the min max values
	data[data!=-9999] = (data[data!=-9999] - recast_minZ) / (recast_maxZ - recast_minZ)

	_fill_arrays(data,displacement_array,color_array, intensity) if custom_colors is None else _fill_arrays_cc(data, custom_colors, displacement_array, color_array, intensity)



	np.save(fprefix + ".npy", data.astype(np.float32))
	np.save(fprefix + "_displacement_array.npy", displacement_array)
	np.save(fprefix + "_color_array.npy", color_array)


	cmd = ""
	cmd += scb.config.query('blender')
	cmd += "  --background --python " + blender_script + " "
	cmd += "--gpu " if gpu else ""
	cmd += "--perspective " if perspective else ""
	cmd += f"--fprefix '{fprefix}' "
	cmd += f"--save_prefix '{save_prefix}' "
	cmd += f"--ortho_scale {ortho_scale} "
	cmd += f"--focal_length {focal_length} "
	cmd += f"--f_stop {f_stop} "
	cmd += f"--dx {dx} "
	cmd += f"--shiftx {shiftx} "
	cmd += f"--shifty {shifty} "
	cmd += f"--camera_tilt {camera_tilt} "
	cmd += f"--camera_rotation {camera_rotation} "
	cmd += f"--sun_tilt {sun_tilt} "
	cmd += f"--sun_rotation {sun_rotation} "
	cmd += f"--sun_intensity {sun_intensity} "
	cmd += f"--exaggeration {exaggeration} "
	cmd += f"--recast_minZ {recast_minZ} "
	cmd += f"--recast_maxZ {recast_maxZ} "
	cmd += f"--number_of_subdivisions {number_of_subdivisions} "
	cmd += f"--res_x {res_x} "
	cmd += f"--res_y {res_y} "
	cmd += f"--samples {samples} "
	cmd += f"--relief {relief} "



	sub.run(cmd, shell = True, check = False)





	# parser.add_argument('-G', '--gpu', action='store_true', help='GPU acceleration')
	# parser.add_argument('-P', '--perspective', action='store_true', help='Camera set to perspective or orthogonal')
	# parser.add_argument('-f', '--fprefix', type=str, help='file name')
	# parser.add_argument('-o', '--save_prefix', type=str, help='file name')
	# parser.add_argument('--ortho_scale', type=float, default = 2., help='file name')
	# parser.add_argument('--focal_length', type=float, default = 40., help='file name')
	# parser.add_argument('--f_stop', type=float, default = 500000., help='file name')
	# parser.add_argument('--dx', type=float, default = 30., help='file name')
	# parser.add_argument('--shiftx', type=float, default = 0., help='file name')
	# parser.add_argument('--shifty', type=float, default = 0., help='file name')
	# parser.add_argument('--camera_tilt', type = float, default = 45.0, help = '')
	# parser.add_argument('--camera_rotation', type = float, default = 0., help = '')
	# parser.add_argument('--sun_tilt', type = float, default = 25.0, help = '')
	# parser.add_argument('--sun_rotation', type = float, default = 315.0, help = '')
	# parser.add_argument('--sun_intensity', type = float, default = 0.2, help = '')
	# parser.add_argument('--exaggeration', type = float, default = 1.5, help = '')
	# parser.add_argument('--recast_minZ', type = float, default = 0., help = '')
	# parser.add_argument('--recast_maxZ', type = float, default = 1.4, help = '')
	# parser.add_argument('--number_of_subdivisions', type = int, default = 1000 , help = '')
	# parser.add_argument('--res_x', type = int, default = 1920 , help = '')
	# parser.add_argument('--res_y', type = int, default = 1080 , help = '')
	# parser.add_argument('--samples', type = int, default = 50 , help = '')
