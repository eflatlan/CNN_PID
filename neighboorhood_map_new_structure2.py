from __future__ import print_function
from statsmodels.stats.diagnostic import lilliefors
from scipy.stats import normaltest, anderson
import numpy as np
from tensorflow.keras.layers import Conv2D, BatchNormalization, LeakyReLU, MaxPooling2D, Dropout, Flatten
import tensorflow as tf
from tensorflow.keras.layers import Conv2D, BatchNormalization, LeakyReLU, MaxPooling2D, Dropout, Flatten
from tensorflow.keras.regularizers import l1, l2, l1_l2

from tensorflow.keras.layers import Conv3D, MaxPooling3D, BatchNormalization, LeakyReLU, Dropout, Flatten
from tensorflow.keras.regularizers import l1, l2, l1_l2


def filter_data(particle_mask, momentum_data, mass, refractive_index_data, highest_range):
    # Flatten the 2D arrays to 1D arrays
    momentum_data = momentum_data.flatten()
    refractive_index_data = refractive_index_data.flatten()

    min_momentum, max_momentum = highest_range
    print(f"Number of True values in particle_mask: {np.sum(particle_mask)}")

    mask_momentum = (momentum_data >= min_momentum) & (momentum_data <= max_momentum)

    mask_final = particle_mask & mask_momentum
    print(f"Number of True values in mask_final: {np.sum(mask_final)}")

    theta_filtered = calculate_theta(momentum_data[mask_final], mass, refractive_index_data[mask_final])

    print(f"Number of non-zero values in theta_filtered: {np.count_nonzero(theta_filtered)}")

    if len(theta_filtered) > 0:
        print(f"Min value of theta_filtered: {np.min(theta_filtered)}")
        print(f"Max value of theta_filtered: {np.max(theta_filtered)}")
        print(f"Differnce : {np.max(theta_filtered)- np.min(theta_filtered)}")
    else:
        print("theta_filtered is empty.")
    return theta_filtered

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
import os
import h5py
import tensorflow as tf

# to check the impact of resolution in the 2d-map;
# print the difference between the filledBins vector versus the map (map is restricted by resolution)
def print_points(filled_bins_array = None, map_array = None, mip_position_array = None, resolution = 10):

    length = map_array.shape[0]
    distances_bins_list = []
    distances_map_list = []

    print(f"filled_bins_array shape = {filled_bins_array.shape}")
    print(f"map_array shape = {map_array.shape}")
    print(f"mip_position_array shape = {mip_position_array.shape}")


    for i in range (1, length):

        filled_bins = np.array(filled_bins_array[i])
        map = np.array(map_array[i, :,:])
        mip_pos = np.array(mip_position_array[i, :])

        #print(f"filled_bins shape = {filled_bins.shape}")
        #print(f"map shape = {map.shape}")
        #print(f"mip_pos shape = {mip_pos.shape}")

        _mip_position = []
        #_mip_position.append(mip_position_array[])
        distances2 = []

        distances_bins = [norm(np.array(pos) - mip_pos) for pos in filled_bins]

        distances_map = []
        for y in range(map.shape[0]):
            for x in range(map.shape[1]):
                if map[y, x] == 1:
                    point = (x, y)
                    distance = np.linalg.norm(np.array(point) - mip_pos*resolution)
                    distances_map.append(distance)



        distances_bins_list.append(distances_bins)
        distances_map_list.append(distances_map)


    # Print the distances for each element in map_data_list
    print(f"Element {i+1} distances:")
    for j, (distances_bins, distances_map) in enumerate(zip(distances_bins_list, distances_map_list)):
        print(f"  Point {j+1}: Distance bins: {distances_bins}\n, Distance map: {distances_map}")
    print()




def plot_maps(filled_bins_array=None, map_array=None, mip_position_array=None, X_momentum=None, X_refractive_index=None, percentage_to_plot=5, resolution = 10):
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
      map_data = map_array[start_index + i, :, :]
      #mass_category = particle_vector[start_index + i].mass_category
      PDG = PDG[start_index + i]
      mip_position = mip_position_array[start_index + i]
      momentum = X_momentum[start_index + i]
      refractive_index = X_refractive_index[start_index + i]

      # Plot the map
      ax.imshow(map_data, cmap='gray')



      #try :
      # Add a red dot at the MIP position
      ax.plot(mip_position[0]*resolution, mip_position[1]*resolution, 'ro')
      #Except exception as e :
      #  print("caught non mip pos ")
      # Set the title with the information
      #ax.set_title(f"Mass: {mass_category}, CKOV: {ckov}, MIP Position: {mip_position}, Momentum: {momentum},  refractive_index: {refractive_index}")
      ax.set_title(f"PDG: {PDG}, MIP Position: {mip_position}, Momentum: {momentum},  refractive_index: {refractive_index}")

      ax.axis('off')

  # Adjust the spacing between subplots
  plt.tight_layout()

  # Show the plot
  plt.show()



import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

cmap = mcolors.LinearSegmentedColormap.from_list("", ["green", "red"])

def plot_individual_types2(idx, x_proton, x_pion, x_kaon, x_non):
    datasets = [x_proton, x_pion, x_kaon, x_non]
    titles = ["Proton", "Pion", "Kaon", "Non"]

    fig, ax = plt.subplots(2, 2, figsize=(14, 14))
    ax = ax.flatten()

    for i, (data, title) in enumerate(zip(datasets, titles)):
        x = data[idx, :, 0]
        y = data[idx, :, 1]
        charge = data[idx, :, 3]
        sc = ax[i].scatter(x, y, c=charge, cmap=cmap)
        plt.colorbar(sc, ax=ax[i], label='Charge', fraction=0.046, pad=0.04)
        ax[i].set_title(title)
        ax[i].set_xlim([0, 130])
        ax[i].set_ylim([0, 130])
        ax[i].set_aspect('equal', adjustable='box')

    plt.tight_layout()
    plt.show()

def plot_combined_types2(i, x_proton, x_pion, x_kaon, x_non):
    fig, ax = plt.subplots(figsize=(7, 7))

    # Plot X_non in red
    x_val = x_non[i, :, 0]
    y_val = x_non[i, :, 1]
    ax.scatter(x_val, y_val, c='red', label='Non', alpha=0.5)

    # Plot X_proton, X_pion, and X_kaon in green
    for x, label in zip([x_proton, x_pion, x_kaon], ['Proton', 'Pion', 'Kaon']):
        x_val = x[i, :, 0]
        y_val = x[i, :, 1]
        ax.scatter(x_val, y_val, c='green', label=label, alpha=0.5)

    ax.legend()
    ax.set_title("Combined Types")
    ax.set_xlim([0, 130])
    ax.set_ylim([0, 130])
    ax.set_aspect('equal', adjustable='box')

    plt.show()


import numpy as np
import numpy as np



# specie_probability : vector containing the specie for same prob as specie in candidate_positions
def extract_neighborhood_map_new(x, y, mip_positions, specie_probability, neighborhood_size, map_size):
    print(f"shape of hadron_candidates {hadron_candidates}")
    # get x, y


    
    #candidate_positions = hc.prob for hc in hadron_candidates
  
    #num_samples = candidate_positions.shape[0]


    print(f"mip_positions shape = {mip_positions.shape}")
    print(f"candidate_positions shape = {candidate_positions.shape}")

    cand_pos = np.asarray([x, y])

    #num_samples = candidate_positions.shape[0]
    num_candidates = candidate_positions.shape[1]

    # cand_pos = np.zeros((num_samples, num_candidates, 2))  # Create an array filled with zeros
    # cand_pos = candidate_positions[:, :, :2]

    print(f"cand_pos shape = {cand_pos.shape}")
    #mip_positions = mip_positions.reshape((num_samples, 1, 2))

    # Use np.tile to replicate along the second dimension (num_candidates)
    mip_positions_expanded = np.tile(mip_positions, (1, num_candidates, 1))
    # Use broadcasting to expand along the second dimension to get shape (8050, 915, 2)
    # Extend mip_positions to match the relevant dimensions of candidate_positions
    centered_positions = cand_pos - mip_positions_expanded


    # Now, extended_mip_positions has a shape of (8050, 1, 2, 4)

    # Perform the subtraction
    print(f"candidate_positions shape = {cand_pos.shape}")
    print(f"extended_mip_positions shape = {mip_positions_expanded.shape}")

    distances = cand_pos - mip_positions_expanded#[:, np.newaxis, :]

    print(f"candidate_positions shape = {candidate_positions.shape}")

    # Calculate distances between candidate positions and MIP positions
    #distances = candidate_positions - mip_positions[:, np.newaxis, :]

    # Calculate the norm of distances to get the Euclidean distance
    distances = np.linalg.norm(distances, axis=-1)

    # Create an empty map
    neighborhood_maps = np.zeros((num_samples, map_size, map_size))

    # Check if the candidate falls within the neighborhood
    mask = distances <= neighborhood_size

    
  
    # Convert centered positions to map indices and shift them to be around the center of the map
    map_indices = np.round(centered_positions[:, :, :2] + map_size // 2).astype(int)
    map_indices = np.clip(map_indices, 0, map_size - 1)

    # Update the map
    # ef :change was done here: take the mask of specie_probability ; take only within +- map_size of MIP for candidate of specie
    neighborhood_maps[np.arange(num_samples)[:, np.newaxis], map_indices[:, :, 1], map_indices[:, :, 0]] = specie_probability[mask]

    # Plotting
    for sample_idx in [1, 2]:
        plt.imshow(neighborhood_maps[sample_idx], cmap='gray')
        plt.colorbar()
        plt.title(f"Neighborhood Map for num_samples = {sample_idx}")
        plt.xlabel("X-axis")
        plt.ylabel("Y-axis")

        # Mark the MIP position (which should be at the center after centering)
        plt.scatter(map_size // 2, map_size // 2, c='red', marker='o')

        plt.show()

    return neighborhood_maps


""""# specie_probability : vector containing the specie for same prob as specie in candidate_positions
def extract_neighborhood_map_new_vectorized(hadron_candidates, mip_positions, specie_probability, neighborhood_size, map_size):

    # get x, y
    candidate_positions = [hc.x, hc.y] for hc in hadron_candidates
    #candidate_positions = hc.prob for hc in hadron_candidates
  
    num_samples = candidate_positions.shape[0]


    print(f"mip_positions shape = {mip_positions.shape}")
    print(f"candidate_positions shape = {candidate_positions.shape}")


    num_samples = candidate_positions.shape[0]
    num_candidates = candidate_positions.shape[1]

    cand_pos = np.zeros((num_samples, num_candidates, 2))  # Create an array filled with zeros
    cand_pos = candidate_positions[:, :, :2]

    # Assuming candidate_positions is an array of shape (8050, 915, 4)
    # Assuming mip_positions is an array of shape (8050, 1, 2, 1)
    print(f"cand_pos shape = {cand_pos.shape}")
    mip_positions = mip_positions.reshape((num_samples, 1, 2))

    # Use np.tile to replicate along the second dimension (num_candidates)
    mip_positions_expanded = np.tile(mip_positions, (1, num_candidates, 1))
    # Use broadcasting to expand along the second dimension to get shape (8050, 915, 2)
    # Extend mip_positions to match the relevant dimensions of candidate_positions
    centered_positions = cand_pos - mip_positions_expanded


    # Now, extended_mip_positions has a shape of (8050, 1, 2, 4)

    # Perform the subtraction
    print(f"candidate_positions shape = {candidate_positions.shape}")
    print(f"extended_mip_positions shape = {mip_positions_expanded.shape}")

    distances = cand_pos - mip_positions_expanded#[:, np.newaxis, :]

    print(f"candidate_positions shape = {candidate_positions.shape}")

    # Calculate distances between candidate positions and MIP positions
    #distances = candidate_positions - mip_positions[:, np.newaxis, :]

    # Calculate the norm of distances to get the Euclidean distance
    distances = np.linalg.norm(distances, axis=-1)

    # Create an empty map
    neighborhood_maps = np.zeros((num_samples, map_size, map_size))

    # Check if the candidate falls within the neighborhood
    mask = distances <= neighborhood_size

    
  
    # Convert centered positions to map indices and shift them to be around the center of the map
    map_indices = np.round(centered_positions[:, :, :2] + map_size // 2).astype(int)
    map_indices = np.clip(map_indices, 0, map_size - 1)

    # Update the map
    # ef :change was done here: take the mask of specie_probability ; take only within +- map_size of MIP for candidate of specie
    neighborhood_maps[np.arange(num_samples)[:, np.newaxis], map_indices[:, :, 1], map_indices[:, :, 0]] = specie_probability[mask]

    # Plotting
    for sample_idx in [1, 2]:
        plt.imshow(neighborhood_maps[sample_idx], cmap='gray')
        plt.colorbar()
        plt.title(f"Neighborhood Map for num_samples = {sample_idx}")
        plt.xlabel("X-axis")
        plt.ylabel("Y-axis")

        # Mark the MIP position (which should be at the center after centering)
        plt.scatter(map_size // 2, map_size // 2, c='red', marker='o')

        plt.show()

    return neighborhood_maps
"""


import tensorflow as tf
from tensorflow.keras import layers, models, Input

import tensorflow as tf
from tensorflow.keras import layers, models, Input


def create_cnn_model(input_shape=None, name="default_name"):
    # Create CNN layers
    cnn_input = Input(shape=input_shape, name=f"input{name}")
    x = layers.Conv2D(32, (3, 1), activation='relu', name=f"conv2D_1_{name}")(cnn_input)
    x = layers.MaxPooling2D((2, 1), name=f"MaxPooling2D_1_{name}")(x)
    x = layers.Conv2D(64, (5, 1), activation='relu', name=f"conv2D_2_{name}")(x)
    x = layers.MaxPooling2D((2, 1), name=f"MaxPooling2D_2_{name}")(x)
    x = layers.Conv2D(64, (7, 1), activation='relu', name=f"conv2D_3_{name}")(x)
    cnn_output = layers.Reshape((-1, 1), name=f"Reshape_{name}")(x)

    return models.Model(inputs=cnn_input, outputs=cnn_output)



