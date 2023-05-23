import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
from __future__ import print_function
import os
import h5py
import tensorflow as tf

# to check the impact of resolution in the 2d-map; 
# print the difference between the filledBins vector versus the map (map is restricted by resolution)
def print_points(filled_bins_array = None, map_array = None, mip_position_array = None):

    length = map_array.shape[0]
    distances_bins_list = []
    distances_map_list = []

    for i in range (1, length):
        
        filled_bins = filled_bins_array[i, :]
        map = map_array[i, :,:]
        mip_pos = mip_position_array[i, :]

        _mip_position = []
        #_mip_position.append(mip_position_array[])
        distances2 = []

        distances_bins = [norm(np.array(pos) - mip_pos) for pos in filled_bins]

        distances_map = []
        for y in range(map.shape[1]):
            for x in range(map.shape[2]):
                if map[i, y, x] == 1:
                    point = (x, y)
                    distance = np.linalg.norm(np.array(point) - mip_pos)
                    distances_map.append(distance)
        
        
        
        distances_bins_list.append(distances_bins)
        distances_map_list.append(distances_map)


    # Print the distances for each element in map_data_list
    print(f"Element {i+1} distances:")
    for j, (distances_bins, distances_map) in enumerate(zip(distances_bins_list, distances_map_list)):
        print(f"  Point {j+1}: Distance bins: {distances_bins}, Distance map: {distances_map}")
    print()




def plot_maps(filled_bins_array=None, map_array=None, mip_position_array=None, X_momentum=None, X_refractive_index=None, X_ckov=None, percentage_to_plot=5):
  """
  Args : filled_bins_array : array that holds the vectors of filled pads 
         map_array : 2d  map with a determined resolution (the points in the filled_bins_array element, just restricted by the resolution)
         mip_position_array : array of the MIP {x, y} positions

         TODO : add mass_category and actual mass?
  """
  

  #percentage_to_plot = 0.05 / 10

  # Calculate the starting index of the samples to plot
  num_samples = map_array.shape[0]
  start_index = -num_samples

  # Create a subplot with the number of rows based on the number of samples
  fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(8, 20))

  # Iterate over the samples and plot each map with information
  for i, ax in enumerate(axes):
      # Get the map and corresponding information
      map_data = map_array[start_index + i, :, :, 0]
      #mass_category = particle_vector[start_index + i].mass_category
      ckov = X_ckov[start_index + i, :]
      mip_position = mip_position_array[start_index + i,:]
      momentum = X_momentum[start_index + i, :]
      refractive_index = X_refractive_index[start_index + i,:]

      # Plot the map
      ax.imshow(map_data, cmap='gray')

      # Add a red dot at the MIP position
      ax.plot(mip_position[0], mip_position[1], 'ro')

      # Set the title with the information
      #ax.set_title(f"Mass: {mass_category}, CKOV: {ckov}, MIP Position: {mip_position}, Momentum: {momentum},  refractive_index: {refractive_index}")
      ax.set_title(f"CKOV: {ckov}, MIP Position: {mip_position}, Momentum: {momentum},  refractive_index: {refractive_index}")

      ax.axis('off')

  # Adjust the spacing between subplots
  plt.tight_layout()

  # Show the plot
  plt.show()
