import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage
import skimage.morphology
import argparse

# Python Imaging Library imports
import Image
import ImageDraw
import os

'''
This is a command line utility for generating the occupancy grid maps necessary for specifying out-of-bounds states.

'''

parser = argparse.ArgumentParser(description='Dilate maps for use in value iteration policy generation.')
parser.add_argument('--discretization', default=0, type=int, help='Number of theta slices to generate')
parser.add_argument('--map', default="./small_kitchen.png", help='Source map to use')
parser.add_argument('--out_path', default="./", help='File prefix to save dilated images with.')
parser.add_argument('--cell_size', default=0.01, type=float, help='Dimension of each pixel in the map (meters)')
parser.add_argument('--padding', default=0.0, type=float, help='Safety padding around the car footprint, only used if --discretization>0')
parser.add_argument('--dilate', default=0.01, type=float, help='Amount to dilate the map by, in meters. Only used if --discretization=0')
parser.add_argument('-s', '--save', help='Whether or not to save the output images', action='store_true')
parser.add_argument('-d', '--display', help='Whether or not to display the output images', action='store_true')

class Car(object):
	""" Wrapper class for utility functions related to the car dimensions """
	def __init__(self):
		# origin is the back right wheel of the car if the car is pointing along the positive x axis

		# 0.3 meters wide, 0.55 meters long
		self.dimensions = np.array([0.55, 0.3])

		# center of lidar, centered side to side and 0.14 meters from front
		self.lidar_position = np.array([0.41, 0.15])

		# point between the rear tires, used as the state in value iteration
		self.base_frame = np.array([0.11, 0.15])

	# theta is the heading of the car, center is the point on the 
	# car which should be at the center of the kernel, if none it
	# is set to the base_frame of the car. cell size is the width/height
	# of one pixel in the output kernel. Padding is how much space should be
	# added to the true dimensions of the car for safety margin (in meters)
	def footprint_kernel(self, theta=0, cell_size=0.01, padding=0.0, center=None):
		dims = self.dimensions / cell_size
		padding = padding / cell_size

		if center == None:
			center = self.base_frame / cell_size

		corners = np.array([[0,0], [dims[0]+padding,0], dims+padding, [0,dims[1]+padding], [0,0]])
		corners -= center
		corners -= padding / 2.0

		# rotation matrix
		R = np.array([[np.cos(-theta), -np.sin(-theta)],
                  	 [np.sin(-theta), np.cos(-theta)]])

		corners = np.dot(R, corners.T).T
		bounds = np.ceil(np.max(np.abs(corners), axis=0)).astype(int)
		kernel = np.zeros(bounds[::-1]*2)
		corners = corners + bounds
		
		# draw car footprint onto kernel
		img = Image.fromarray(kernel)
		draw = ImageDraw.Draw(img)
		draw.polygon([tuple(p) for p in corners], fill=1)
		kernel = np.asarray(img)
		return kernel

def rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.299, 0.587, 0.114])

if __name__ == '__main__':

	args = parser.parse_args()

	disc = args.discretization

	image = rgb2gray(scipy.ndimage.imread(args.map))

	if disc == 0:
		out_img = "small_kitchen_dilated.png"
		amount = args.dilate / args.cell_size

		og = np.zeros_like(image, dtype=float)
		og[image==0] = 255

		el = skimage.morphology.disk(amount)
		dilated = -skimage.morphology.dilation(og, selem=el)

		if args.save:
			scipy.misc.imsave(args.out_path, dilated)
		if args.display:
			plt.imshow(dilated, cmap="gray")
			plt.show()
	else:
		car = Car()
		thetas = np.linspace(0,2.0*np.pi, disc, endpoint=False)

		if args.save:
			summary = open(args.out_path + "_summary.yaml", "w")
			summary.write("discretization: " + str(disc) + "\n")
			summary.write("map: " + os.path.abspath(args.map) + "\n")
			summary.write("files: \n")
		i = 0
		for a in thetas:
			i += 1
			kernel = car.footprint_kernel(cell_size=args.cell_size, theta=a, padding=args.padding)

			og = np.zeros_like(image, dtype=float)
			og[image==0] = 255

			dilated = -skimage.morphology.dilation(og, selem=kernel)

			if args.save:
				file_path = args.out_path + "_" + str(i) + "_of_" + str(disc) + ".png"
				file_name = os.path.abspath(file_path)
				scipy.misc.imsave(file_path, dilated)
				summary.write("    -" + file_name + "\n")
			
			if args.display:
				plt.imshow(dilated, cmap="gray")
				plt.show()

		if args.save:
			summary.close()


		# information of interest:
		# - file list
		# - discretization