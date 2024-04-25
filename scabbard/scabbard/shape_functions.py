'''
Collections of tools related to kernels
For example Kernel Density Estimator
Spreading a values following a kernel on a numpy array
...
'''

import numpy as np


def sinusoidal_wave(X, a, b, L):
    # Compute Y values for the sinusoidal wave
    Y = (b - a) / 2 * np.sin(2 * np.pi * X / L) + (a + b) / 2
    return Y

def generate_u_shaped_sloped_surface(nx, ny, dx, dy, slope=0.1, Umag = 0.8):
    """
    Generates a 2D numpy array describing a U-shaped surface that is concave (facing upwards)
    and has a slope making the northern part the highest point.
    
    Parameters:
    - nx, ny: Number of cells in the east-west (x) and north-south (y) directions.
    - dx, dy: Cell size in the east-west (x) and north-south (y) directions.
    - slope: Slope of the surface in the north-south direction, making the north higher.
    
    Returns:
    - 2D numpy array representing the surface elevation.
    """
    
    # Generate meshgrid for the positions
    x, y = np.meshgrid(np.arange(nx) * dx, np.arange(ny) * dy)
    
    # Center the x values around 0 for the U shape calculation
    x_centered = x - (nx * dx) / 2
    
    # Calculate the U shape using a negative quadratic function (for concave shape)
    u_shape = -x_centered**2
    
    # Normalize U shape to have values between 0 and 1
    u_shape_normalized = (u_shape - np.min(u_shape)) / (np.max(u_shape) - np.min(u_shape)) * Umag
    
    # Modify the slope effect to make the northern part the highest
    slope_effect = (ny * dy - y) * slope
    
    # Combine the U shape and the slope to generate the surface
    surface = u_shape_normalized*-1 + slope_effect
    
    # Normalize the surface to have a minimum of 0 (optional)
    surface -= np.min(surface)
    
    return surface

def gaussian_spread_on_1D(X = np.linspace(0, 100-1, 100), M = 50, x_c = 50, sigma = 10):

	# Calculate Gaussian
	gaussian = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((X - x_c) / sigma) ** 2)

	# Scale the Gaussian so its sum equals M
	gaussian *= M / np.sum(gaussian)

	# Verify the sum of the Gaussian is approximately M
	print("Sum of Gaussian values:", np.sum(gaussian))

	return gaussian
