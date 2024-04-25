################################################################################
# Import Libraries
import bpy
import numpy as np
import math
from mathutils import Matrix
import importlib
import sys
import os
import argparse


# Create the parser
parser = argparse.ArgumentParser(description='Process some arguments.')

# Add a boolean argument with both short (-b) and long (--flag) flags
# parser.add_argument('-b', "--blenderoptionignore1", action='store_true', help='GPU acceleration')
# parser.add_argument('-P', "--blenderoptionignore2", action='store_true', help='GPU acceleration')
parser.add_argument('--gpu', action='store_true', help='GPU acceleration')
parser.add_argument('--perspective', action='store_true', help='Camera set to perspective or orthogonal')
parser.add_argument('--fprefix', type=str, help='file name')
parser.add_argument('--save_prefix', type=str, help='file name')
parser.add_argument('--ortho_scale', type=float, default = 2., help='file name')
parser.add_argument('--focal_length', type=float, default = 40., help='file name')
parser.add_argument('--f_stop', type=float, default = 500000., help='file name')
parser.add_argument('--dx', type=float, default = 30., help='file name')
parser.add_argument('--shiftx', type=float, default = 0., help='file name')
parser.add_argument('--shifty', type=float, default = 0., help='file name')
parser.add_argument('--camera_tilt', type = float, default = 45.0, help = '')
parser.add_argument('--camera_rotation', type = float, default = 0., help = '')
parser.add_argument('--sun_tilt', type = float, default = 25.0, help = '')
parser.add_argument('--sun_rotation', type = float, default = 315.0, help = '')
parser.add_argument('--sun_intensity', type = float, default = 0.2, help = '')
parser.add_argument('--exaggeration', type = float, default = 1.5, help = '')
parser.add_argument('--recast_minZ', type = float, default = 0., help = '')
parser.add_argument('--recast_maxZ', type = float, default = 1.4, help = '')
parser.add_argument('--number_of_subdivisions', type = int, default = 1000 , help = '')
parser.add_argument('--res_x', type = int, default = 1920 , help = '')
parser.add_argument('--res_y', type = int, default = 1080 , help = '')
parser.add_argument('--samples', type = int, default = 50 , help = '')
parser.add_argument('--relief', type = float, default = 50. , help = '')







args, unknown = parser.parse_known_args()

print(args)
print('\n\n\n')
print(unknown)
print('\n\n\n')


################################################################################
# This sets the render enginer to CYCLES.
# In order for TopoBlender to render your topography correctly, you must use
# this render engine. This engine is optimized for GPUs, so if your computer
# lacks a GPU, TopoBlender may be slow.
GPU_boolean = args.gpu
delete_all = True

# define current scene
scn = bpy.context.scene
if(delete_all):
    # Select all objects in the scene
    bpy.ops.object.select_all(action='SELECT')

# Delete selected objects
bpy.ops.object.delete()

# check the render engine, and change to CYCLES
if not scn.render.engine == 'CYCLES':
    scn.render.engine = 'CYCLES'

#if cycles, change to gpu rendering if user selects GPU
if scn.render.engine == 'CYCLES':
    if GPU_boolean == 1:
        scn.cycles.device = 'GPU'


# Loading the file
data = np.load(args.fprefix + ".npy")
displacement_array = np.load(args.fprefix + "_displacement_array.npy")
color_array = np.load(args.fprefix + "_color_array.npy")

output = args.save_prefix + ".png"

ny,nx = data.shape
dx,dy = args.dx, args.dx

# Standard nodata value
nodata_value = -9999


##########camera parameters##########
#camera type
camera_type = 'perspective' if args.perspective else 'orthogonal' 

#orthogonal
ortho_scale = args.ortho_scale #when using orthogonal scale, increase to "zoom" out


focal_length = args.focal_length #mm when using perspective camera, increase to zoom in
f_stop = args.f_stop # affects depths of field, lower for a shallow dof, higher for wide dof
shift_x = args.shiftx # you may need to shift the camera to center the topo in the frame
shift_y = args.shifty # you may need to shift the camera to center the topo in the frame

#camera location
camera_tilt = args.camera_tilt # #degrees from horizontal
camera_rotation = args.camera_rotation # #camera location degrees CW from North
######################################


##########sun properties##############
sun_tilt = args.sun_tilt # #degrees from horizontal
sun_rotation = args.sun_rotation # #degrees CW from North
sun_intensity = args.sun_intensity # #sun intensity
######################################


#####landscape representation#########
number_of_subdivisions = args.number_of_subdivisions # of subvisions, more increases the detail
exaggeration = args.exaggeration # #vertical exaggeration
recast_minZ = args.recast_minZ #
recast_maxZ = args.recast_maxZ #
luminosity_scale = 0.7
######################################

#########render settings##############
res_x = args.res_x #x resolution of the render
res_y = args.res_y #y resolution of the render
samples = args.samples #number of samples that decides how "good" the render looks. more is better but takes longer
######################################

#convert nodata value to nan
data[data==nodata_value] = np.nan
# data y_length we use this to scale the topography
y_length = (2.0+float(data.shape[1])) * dx


# create an image that holds the data
# displacement_array = np.zeros(((data.shape[0]+2)*(data.shape[1]+2)*4), dtype=np.float32)
# color_array = np.zeros(((data.shape[0]+2)*(data.shape[1]+2)*4), dtype=np.float32)

# # Add a 1-pixel border around the entire topography and rescale the topgraphy
# # so minimum elevation is zero and the maximum elevation is 1.
# _fill_arrays(data,displacement_array,color_array, luminosity_scale)

# Create an image from the datafile        
displacement_image = bpy.data.images.new('displacement_data', data.shape[0]+2, data.shape[1]+2, alpha=False, float_buffer = True)
color_image = bpy.data.images.new('color_data', data.shape[0]+2, data.shape[1]+2, alpha=False, float_buffer = True)
# Fast way to set pixels (since 2.83)
displacement_image.pixels.foreach_set(displacement_array)
color_image.pixels.foreach_set(color_array)
# Pack the image into .blend so it gets saved with it
displacement_image.pack()
color_image.pack()
print("C") 
#make plane
plane_size = 1.0
topo_mesh = bpy.ops.mesh.primitive_plane_add(size=plane_size)
topo_obj = bpy.context.active_object

#Change the shape of the object to match data_aspect
topo_obj.scale = ((2.0+data.shape[0])/(2.0+data.shape[1]),1,1)

#add material
topo_mat = bpy.data.materials.new("topo_mat")
topo_mat.cycles.displacement_method = "DISPLACEMENT"
topo_mat.use_nodes = True

#calculate subdivisions
order_of_magnitude = math.floor(math.log10(number_of_subdivisions))
first_digit = int(np.round(number_of_subdivisions / (10.0 ** order_of_magnitude)))

topo_obj.data.materials.append(topo_mat)
bpy.ops.object.mode_set(mode="EDIT")
for i in range(0,order_of_magnitude):
    bpy.ops.mesh.subdivide(number_cuts=10)
bpy.ops.mesh.subdivide(number_cuts=first_digit)
bpy.ops.object.mode_set(mode="OBJECT")

#add image node - determines the displacement
displacement_image_node = topo_mat.node_tree.nodes.new("ShaderNodeTexImage")
#assign png to image to node
displacement_image_node.image = displacement_image
#change colorspace to b&w
displacement_image_node.image.colorspace_settings.name="Linear FilmLight E-Gamut"

#add image node - determines the color of the lanscape
color_image_node = topo_mat.node_tree.nodes.new("ShaderNodeTexImage")
#assign png to image to node
color_image_node.image = color_image
#change colorspace to b&w
color_image_node.image.colorspace_settings.name="Linear FilmLight E-Gamut"
    
#add displacement node - done
displacement_node = topo_mat.node_tree.nodes.new("ShaderNodeDisplacement")
displacement_node.inputs.get("Scale").default_value = exaggeration * plane_size * args.relief / y_length
displacement_node.inputs.get("Midlevel").default_value = 0.0

#add world sky - done
topo_world = bpy.data.worlds.new('topo_world')
topo_world.use_nodes = True
topo_world_node = topo_world.node_tree.nodes.new("ShaderNodeTexSky")
topo_world_node.sun_elevation = np.radians(sun_tilt)
topo_world_node.sun_rotation = np.radians(sun_rotation-90.0)
topo_world_node.sun_intensity = sun_intensity
topo_world.node_tree.links.new(topo_world_node.outputs['Color'], topo_world.node_tree.nodes['Background'].inputs[0])

#add camera - done
camera_distance = 2.0 #meters
cam = bpy.data.cameras.new('topo_cam')
cam_obj = bpy.data.objects.new('topo_cam',cam)
cam_obj.rotation_euler = (np.radians(90.0 - camera_tilt), np.radians(0), np.radians(270.0 - camera_rotation))
cam_obj.matrix_basis @= Matrix.Translation((0.0, 0.0, camera_distance))

if camera_type == 'orthogonal':
    bpy.data.cameras['topo_cam'].type = 'ORTHO'
    bpy.data.cameras['topo_cam'].ortho_scale = ortho_scale
elif camera_type == 'perspective':
    bpy.data.cameras['topo_cam'].type = 'PERSP'
    bpy.data.cameras['topo_cam'].lens = focal_length
    bpy.data.cameras['topo_cam'].shift_x = shift_x
    bpy.data.cameras['topo_cam'].shift_y = shift_y
    bpy.data.cameras['topo_cam'].dof.use_dof = True
    bpy.data.cameras['topo_cam'].dof.aperture_fstop = f_stop
    bpy.data.cameras['topo_cam'].dof.focus_distance = camera_distance
    
#connect nodes - done
topo_mat.node_tree.links.new(displacement_image_node.outputs["Color"], \
                             displacement_node.inputs["Height"])
topo_mat.node_tree.links.new(displacement_node.outputs["Displacement"], \
                             topo_mat.node_tree.nodes["Material Output"].inputs["Displacement"])
topo_mat.node_tree.links.new(color_image_node.outputs["Color"], \
                             topo_mat.node_tree.nodes["Principled BSDF"].inputs[0])

bpy.context.scene.collection.objects.link(cam_obj)                         
bpy.context.scene.camera = cam_obj    
bpy.context.scene.world = topo_world

#render settings
bpy.context.scene.render.engine = 'CYCLES'
bpy.context.scene.cycles.samples = samples
bpy.context.scene.render.resolution_x = int(res_x)
bpy.context.scene.render.resolution_y = int(res_y)
bpy.context.scene.render.filepath = output
bpy.ops.render.render(write_still=True)