from __future__ import print_function
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




def plot_maps(filled_bins_array=None, map_array=None, mip_position_array=None, X_momentum=None, X_refractive_index=None, X_ckov=None, percentage_to_plot=5, resolution = 10):
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
      ckov = X_ckov[start_index + i]
      mip_position = mip_position_array[start_index + i]
      momentum = X_momentum[start_index + i]
      refractive_index = X_refractive_index[start_index + i]

      # Plot the map
      ax.imshow(map_data, cmap='gray')

      # Add a red dot at the MIP position
      ax.plot(mip_position[0]*resolution, mip_position[1]*resolution, 'ro')

      # Set the title with the information
      #ax.set_title(f"Mass: {mass_category}, CKOV: {ckov}, MIP Position: {mip_position}, Momentum: {momentum},  refractive_index: {refractive_index}")
      ax.set_title(f"CKOV: {ckov}, MIP Position: {mip_position}, Momentum: {momentum},  refractive_index: {refractive_index}")

      ax.axis('off')

  # Adjust the spacing between subplots
  plt.tight_layout()

  # Show the plot
  plt.show()


  
def calc_dist2mip(maps = None,  mip_positions = None, resolution = 10):
#        X_dist2mip = calc_dist2mip(maps = X_map, mip_positions = X_mip_position, resolution = resolution)

    length = maps.shape[0]
    distances_map_list = []



    #i = 0
    for i in range (length):

        map = np.array(maps[i, :,:])
        mip_pos = np.array(mip_positions[i, :]).copy()

        #print(f"filled_bins shape = {filled_bins.shape}")
        #print(f"map shape = {map.shape}")
        #print(f"mip_pos shape = {mip_pos.shape}")

        _mip_position = []


        # start vectorization
        # Assuming your map is a numpy array
        indices = np.where(map == 1)

        # Now indices[0] contains the y-indices and indices[1] contains the x-indices
        points = np.stack(indices, axis=-1)  # Shape is [num_points, 2]
        mip_pos = np.array([mip_pos[1], mip_pos[0]]) * resolution
        distances = np.linalg.norm(points - mip_pos, axis=-1)
        distances = distances[distances < 20*resolution]  # Filter out large distances
                                                          # later do this instead by imposing the masshypothesis
        # end vectorization

        #distances_map = []
        #for y in range(map.shape[0]):
        #    for x in range(map.shape[1]):
        #        if map[y, x] == 1:
        #            point = (x, y)
        #            distance = np.linalg.norm(np.array(point) - mip_pos*resolution)
        #            if distance < 80*resolution: # dont add if distance is unreasonably large
        #              distances_map.append(distance)
      
      
        distances_map_list.append(distances)
       # print(f"calc_dist2mip : distances : {distances}")
       # print(f"calc_dist2mip : mip_pos * resolution : {mip_pos}")
       # print(f"calc_dist2mip : points : {points}")

    # TODO: NB this should be removed, if i dont add it distances_map_list is one elem shorter than the other dataframes
      
    #distances_map_list.append(temp)


    print(f"maps shape = {maps.shape}")
    print(f"mip_position_array shape = {mip_positions.shape}")
    print(f"distances_map_list shape = {np.array(distances_map_list, dtype=object).shape}")

    return distances_map_list



import matplotlib.pyplot as plt
import random

def plot_random_element(X_train_map):
    index = random.randint(0, len(X_train_map) - 1)  # Pick a random index
    element = X_train_map[index, :, :, 0]  # Retrieve the element

    plt.figure(figsize=(8, 6))
    plt.imshow(element, cmap='viridis', origin='lower')
    plt.title(f"Random Element from X_train_map (Index {index})")
    plt.colorbar(label='Intensity')
    plt.xlabel('X Axis')
    plt.ylabel('Y Axis')
    plt.show()


    

def plot5(X_test_map, particle_vector):

  # Plotting random maps with information

  # Select 5 random indices from the test data
  random_indices = np.random.choice(range(X_test_map.shape[0]), size=5, replace=False)

  # Create a subplot with 5 rows and 1 column
  fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(8, 20))

  # Iterate over the random indices and plot each map with information
  for i, index in enumerate(random_indices):
      # Get the map and corresponding information
      map_data = X_test_map[index, :, :, 0]
      mass_category = particle_vector[index].mass_category
      ckov = particle_vector[index].ckov
      mip_position = particle_vector[index].mip_position
      momentum = particle_vector[index].momentum
      
      # Plot the map
      axes[i].imshow(map_data, cmap='gray')
      
      # Add a red dot at the MIP position
      axes[i].plot(mip_position[0], mip_position[1], 'ro')
      
      # Set the title with the information    a
      #x.set_title(f"Mass: {mass_category}, CKOV: {ckov}, MIP Position: {mip_position:.4f}, Momentum: {momentum:.4f}")
      mip_pos = f"{mip_position:.4f}"
      axes[i].set_title(f"Mass: {mass_category}, CKOV: {ckov}, MIP Position: {mip_pos}, Momentum: {momentum}")
      axes[i].axis('off')

  # Adjust the spacing between subplots
  plt.tight_layout()

  # Show the plot
  plt.show()




  def create_lr_scheduler(num_epochs = 10):

  start_lr = 0.1
  end_lr = 5e-6
  exp_decay = -np.log(end_lr/start_lr) / num_epochs # Calculate decay rate based on start and end learning rate

  lr_scheduler = tf.keras.callbacks.LearningRateScheduler(lambda epoch: start_lr * np.exp(-exp_decay * epoch)) 
  return lr_scheduler





def plot_lr(num_epochs = 10, history = None):
  div = num_epochs/4
  lrs = 1e-4 * (10 ** (np.arange(num_epochs)/div))
  plt.figure(figsize=(10, 7))
  plt.semilogx(lrs, history.history["loss"]) # we want the x-axis (learning rate) to be log scale
  plt.xlabel("Learning Rate")
  plt.ylabel("Loss")


  plt.title("Learning rate vs. loss")



def create_lr_scheduler(num_epochs = 10):

  start_lr = 0.1
  end_lr = 5e-6
  exp_decay = -np.log(end_lr/start_lr) / num_epochs # Calculate decay rate based on start and end learning rate

  lr_scheduler = tf.keras.callbacks.LearningRateScheduler(lambda epoch: start_lr * np.exp(-exp_decay * epoch)) 
  return lr_scheduler





def plot_lr(num_epochs = 10, history = None):
  div = num_epochs/4
  lrs = 1e-4 * (10 ** (np.arange(num_epochs)/div))
  plt.figure(figsize=(10, 7))
  plt.semilogx(lrs, history.history["loss"]) # we want the x-axis (learning rate) to be log scale
  plt.xlabel("Learning Rate")
  plt.ylabel("Loss")


  plt.title("Learning rate vs. loss");



#def plot_worst_(model, y_test, X_test_map, X_test_momentum, X_test_refractive_index, X_test_ckov, X_test_mip_position, y_pred):
def plot_maps(filled_bins_array=None, map_array=None, mip_position_array=None, X_momentum=None, X_refractive_index=None, X_ckov=None, percentage_to_plot=5, resolution = 10):
  #  print("Shape of y_pred: ", y_pred.shape)
  # 1. Predict labels on validation data
  #plot_worst(model, y_test, X_test["X_test_map"], X_test["X_test_momentum"], X_test["X_test_refractive_index"], X_test["X_test_ckov"], X_test["X_test_mip_position"], y_pred_test)

  # 2. Calculate the difference between predicted and actual labels
  losses = tf.keras.losses.categorical_crossentropy(y_test, y_pred).numpy()

  # Sort the indices of the losses from highest to lowest
  sorted_indices = np.argsort(losses)[::-1]

  # Get the indices of the worst performing 10%
  worst_10_percent_indices = sorted_indices[:int(0.1*len(sorted_indices))]

  # Create figure and axes
  num_plots = len(worst_10_percent_indices)
  #fig, axes = plt.subplots(num_plots, 1, figsize=(8, 20))
  fig, axes = plt.subplots(num_plots,figsize=(8, 20))

  # Define mass categories
  mass_categories = ["pion", "kaon", "proton"]

  # 3. Create plots for these cases, including their feature information and predicted vs actual labels
  for i, index in enumerate(worst_10_percent_indices):
      # Get the map and corresponding information
      map_data = X_test_map[index, :, :]
      actual_mass_category = mass_categories[np.argmax(y_test[index])]

      print(f"y_test[index] = {y_test[index]}")

      predicted_mass_category = mass_categories[np.argmax(y_pred[index])]
      ckov = X_test_ckov[index]
      mip_position = X_test_mip_position[index]
      momentum = X_test_momentum[index]
      refractive_index = X_test_refractive_index[index]
      
      mass_actual = momentum * np.sqrt(refractive_index**2 * np.cos(ckov)*np.cos(ckov) - 1)
      
      # Check if the value is NaN (invalid Cherenkov angle)
      if np.isnan(mass_actual):
          mass_actual = "Invalid"

      # Plot the map
      axes[i].imshow(map_data, cmap='gray')

      # Add a red dot at the MIP position
      axes[i].plot(mip_position[0]*resolution, mip_position[1]*resolution, 'ro')

      # Set the title with the information
      axes[i].set_title(f"Actual Mass")#: {actual_mass_category}, Predicted Mass: {predicted_mass_category},\nMass: {mass_actual}, Mass_prob = {y_pred[index]} \nCKOV: {ckov}, MIP Position: {mip_position}, \nMomentum: {momentum}, Refractive Index: {refractive_index}")
      #
      axes[i].set_title(f"Actual Mass: {actual_mass_category}, Predicted Mass: {predicted_mass_category},\nMass: {mass_actual}, Mass_prob = {y_pred[index]} \nCKOV: {ckov}, MIP Position: {mip_position}, \nMomentum: {momentum}, Refractive Index: {refractive_index}")

      #axes[i].set_title(f"Actual Mass: {actual_mass_category}, Predicted Mass: {predicted_mass_category}, Mass: {mass_actual}\nCKOV: {ckov}, MIP Position: {mip_position}, Momentum: {momentum}, Refractive Index: {refractive_index}")
      axes[i].axis('off')

      print("\n")
      print(f"  Actual Mass: {actual_mass_category}, Predicted Mass: {predicted_mass_category},\n Mass: {mass_actual}, Mass_prob = {y_pred[index]} \n CKOV: {ckov}, MIP Position: {mip_position}, \n  Momentum: {momentum}, Refractive Index: {refractive_index}")
  # Adjust the spacing between subplots
  plt.tight_layout()

  # Show the plot
  plt.show()
def create_circular_mask(center, size, mean_radius, std):
    H, W = size
    Y, X = np.ogrid[:H, :W]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = (dist_from_center >= mean_radius - 2*std) & (dist_from_center <= mean_radius + 2*std)
    return mask


def extract_segment_around_mip(mip_positions = None, window_sizes = None, maps = None, std = 7):
    """
    Args : mip_positions : array of the MIP-positions, (num_samples, {x, y})
           window_sizes : radius of the segments (num_samples, 3) 3 = number of particle classes
           maps : the photon hit-maps (num_samples, 144*resolution, 160*resolution)
           std : standard deviatons to be applied to the ring-radius 
           
    Returns : the extracted regions for each of the different radiuses (m_pion, m_kaon, m_proton), masked with the map
    """
    
    windows = []

    for mip_position, window_size, map in zip(mip_positions, window_sizes, maps):
        radius = np.mean(window_size)
        #std = np.std(window_size)

        mask = create_circular_mask(mip_position, map.shape, radius, std)

        window = np.where(mask, map, 0)

        windows.append(window)

    return np.array(windows)




# create a map, the resolution is the "inverse" 
def create_map(filledBins=None, resolution=4):
    map_shape = (int(144 * resolution), int(160 * resolution))
    map_data = np.zeros(map_shape, dtype=np.int32)
   
    if filledBins is not None:
        filledBins_np = np.array(filledBins)
        indices = (filledBins_np * resolution).astype(int)
        map_data[indices[:, 1], indices[:, 0]] = 1

    return map_data