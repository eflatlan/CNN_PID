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

def build_species_layers2(input_map, filters, filter_sizes, stride_arr, alpha, dropout_rate, l1_reg=0.0, l2_reg=0.0):

    if l1_reg > 0.0 and l2_reg > 0.0:
        reg = l1_l2(l1=l1_reg, l2=l2_reg)
    elif l1_reg > 0.0:
        reg = l1(l1_reg)
    elif l2_reg > 0.0:
        reg = l2(l2_reg)
    else:
        reg = None

    layers = input_map

    for i in range(len(filters)):
        layers = Conv3D(filters[i], 
                        filter_sizes[i],  
                        strides=stride_arr[i], 
                        padding='same', 
                        kernel_regularizer=reg)(layers)

        
        layers = BatchNormalization()(layers)
        layers = LeakyReLU(alpha=alpha)(layers)
        layers = MaxPooling3D((2, 2, 2))(layers)  # Updated to MaxPooling3D
        layers = Dropout(dropout_rate)(layers)

    flat_map = Flatten()(layers)
    return flat_map


def build_species_layers(input_map, filters, filter_sizes, stride_arr, alpha, dropout_rate, l1_reg=0.0, l2_reg=0.0):

    if l1_reg > 0.0 and l2_reg > 0.0:
        reg = l1_l2(l1=l1_reg, l2=l2_reg)
    elif l1_reg > 0.0:
        reg = l1(l1_reg)
    elif l2_reg > 0.0:
        reg = l2(l2_reg)
    else:
        reg = None
    layers = input_map
    for i in range(len(filters)):
        layers = Conv2D(filters[i], 
                        (filter_sizes[i], filter_sizes[i]), 
                        strides=(stride_arr[i], stride_arr[i]), 
                        padding='same', 
                        kernel_regularizer=reg)(layers)
        layers = BatchNormalization()(layers)
        layers = LeakyReLU(alpha=alpha)(layers)
        layers = MaxPooling2D((2, 2))(layers)
        layers = Dropout(dropout_rate)(layers)
    flat_map = Flatten()(layers)
    return flat_map
def eval_data(data):
    print(f"Data shape: {data.shape}")

    # Check normality using normaltest from scipy
    stat, p = normaltest(data)
    print(f"Normaltest P-value: {p}")
    if p > 0.05:
        print('normaltest > Probably Gaussian')
    else:
        print('normaltest > Probably not Gaussian')

    # Check normality using lilliefors from statsmodels
    stat, p = lilliefors(data)
    print(f"Lilliefors P-value: {p}")
    if p > 0.05:
        print('lilliefors > Probably Gaussian')
    else:
        print('lilliefors > Probably not Gaussian')

    # Check normality using anderson from scipy
    result = anderson(data)
    print(f"Anderson Statistic: {result.statistic}")
    for i in range(len(result.critical_values)):
        sl, cv = result.significance_level[i], result.critical_values[i]
        if result.statistic < cv:
            print(f'anderson > Probably Gaussian at the {sl}% level')
        else:
            print(f'anderson > Probably not Gaussian at the {sl}% level')



# Function to calculate theta
def calculate_theta(momentum, mass, refractive_index):
    return np.arccos(np.sqrt(momentum**2 + mass**2) / (momentum * refractive_index))



def find_highest_frequency_momentum_range(momentum_data, min_momentum=2.5, max_momentum=3.0, window=0.05):
    highest_frequency = 0
    highest_range = (0, 0)

    for start in np.arange(min_momentum, max_momentum - window + 0.01, window):
        end = start + window
        count = np.sum((momentum_data >= start) & (momentum_data < end))

        if count > highest_frequency:
            highest_frequency = count
            highest_range = (start, end)

    return highest_range

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

# # Function to plot histogram and Gaussian fit
# def plot_hist_and_fit(ax, data, color, label):
#     n, bins, patches = ax.hist(data, bins=40, alpha=0.6, label=label, color=color, density=True)
#     mu, std = norm.fit(data)
#     x = np.linspace(min(data), max(data), 100)
#     p = norm.pdf(x, mu, std)
#     ax.plot(x, p, color=color, linestyle='dashed', linewidth=2)
#     #ax.text(0.6, 0.8, f'Std: {std:.2f}', transform=ax.transAxes)
#     for i in [-3, 3]:
#         ax.axvline(mu + i * std, color=color, linestyle='dashed', linewidth=1)

#     std_mrad = std * 1e3

#     # Plot vertical lines and add text for 2 std and 3 std
#     for num_std in [2, 3]:
#         ax.axvline(mu - num_std * std, color='purple', linestyle='--')
#         ax.axvline(mu + num_std * std, color='purple', linestyle='--')
#         ax.text(mu + num_std * std  , 1.085 * max(n), f'{num_std} std', color='purple')

#     # Display standard deviation in mRad
#     ax.text(0.5, 0.875, f'\u03C3 = {std_mrad:.2f} [mRad]', transform=ax.transAxes, color='black')
#     ax.set_title(label)

# # Constants
# mass_pion = 0.1396  # GeV/c^2
# mass_kaon = 0.4937  # GeV/c^2
# mass_proton = 0.9383  # GeV/c^2


# highest_range = find_highest_frequency_momentum_range(X_train_momentum)
# print("Highest frequency momentum range:", highest_range)

# highest_range = (2.895, 2.9)
# # Filter data using this momentum range
# theta_train_pion_filtered = filter_data(y_train[:, 0] == 1, X_train_momentum, mass_pion, X_train_refractive_index, highest_range)
# theta_train_kaon_filtered = filter_data(y_train[:, 1] == 1, X_train_momentum, mass_kaon, X_train_refractive_index, highest_range)
# theta_train_proton_filtered = filter_data(y_train[:, 2] == 1, X_train_momentum, mass_proton, X_train_refractive_index, highest_range)

# # Calculate theta for training and test sets
# # theta_train_pion_filtered = filter_data(y_train[:, 0] == 1, X_train_momentum, mass_pion, X_train_refractive_index)
# # theta_train_kaon_filtered = filter_data(y_train[:, 1] == 1, X_train_momentum, mass_kaon, X_train_refractive_index)
# # theta_train_proton_filtered = filter_data(y_train[:, 2] == 1, X_train_momentum, mass_proton, X_train_refractive_index)

# theta_test_pion_filtered = filter_data(y_test[:, 0] == 1, X_test_momentum, mass_pion, X_test_refractive_index, highest_range)
# theta_test_kaon_filtered = filter_data(y_test[:, 1] == 1, X_test_momentum, mass_kaon, X_test_refractive_index, highest_range)
# theta_test_proton_filtered = filter_data(y_test[:, 2] == 1, X_test_momentum, mass_proton, X_test_refractive_index, highest_range)

# # Create 2x3 subplots for histograms
# fig, axs = plt.subplots(2, 3, figsize=(18, 12))

# # Plot histograms with Gaussian fits for the training set
# plot_hist_and_fit(axs[0, 0], theta_train_proton_filtered, 'b', 'Proton (Training)')
# plot_hist_and_fit(axs[0, 1], theta_train_kaon_filtered, 'g', 'Kaon (Training)')
# plot_hist_and_fit(axs[0, 2], theta_train_pion_filtered, 'r', 'Pion (Training)')

# # Plot histograms with Gaussian fits for the test set
# plot_hist_and_fit(axs[1, 0], theta_test_proton_filtered, 'b', 'Proton (Test)')
# plot_hist_and_fit(axs[1, 1], theta_test_kaon_filtered, 'g', 'Kaon (Test)')
# plot_hist_and_fit(axs[1, 2], theta_test_pion_filtered, 'r', 'Pion (Test)')

# # Labels and titles
# for i in range(2):
#     for j in range(3):
#         axs[i, j].set_xlabel('Theta (radians)')
#         axs[i, j].set_ylabel('Frequency')
#         axs[i, j].legend()

# plt.tight_layout()
# plt.show()








# # # Masses in GeV/c^2
# # mass_pion = 0.1396
# # mass_kaon = 0.4937
# # mass_proton = 0.9383

# # Calculate theta for training set and test set
# def calculate_theta(momentum, mass, refractive_index):
#     return np.arccos(np.sqrt(momentum**2 + mass**2) / (momentum * refractive_index))

# # Training set
# theta_train_pion = calculate_theta(X_train_momentum[y_train[:, 0] == 1], mass_pion, X_train_refractive_index[y_train[:, 0] == 1])
# theta_train_kaon = calculate_theta(X_train_momentum[y_train[:, 1] == 1], mass_kaon, X_train_refractive_index[y_train[:, 1] == 1])
# theta_train_proton = calculate_theta(X_train_momentum[y_train[:, 2] == 1], mass_proton, X_train_refractive_index[y_train[:, 2] == 1])

# # Test set
# theta_test_pion = calculate_theta(X_test_momentum[y_test[:, 0] == 1], mass_pion, X_test_refractive_index[y_test[:, 0] == 1])
# theta_test_kaon = calculate_theta(X_test_momentum[y_test[:, 1] == 1], mass_kaon, X_test_refractive_index[y_test[:, 1] == 1])
# theta_test_proton = calculate_theta(X_test_momentum[y_test[:, 2] == 1], mass_proton, X_test_refractive_index[y_test[:, 2] == 1])

# fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# # Training set
# axs[0].scatter(X_train_momentum[y_train[:, 0] == 1], theta_train_pion, label='Pion', alpha=0.6)
# axs[0].scatter(X_train_momentum[y_train[:, 1] == 1], theta_train_kaon, label='Kaon', alpha=0.6)
# axs[0].scatter(X_train_momentum[y_train[:, 2] == 1], theta_train_proton, label='Proton', alpha=0.6)
# axs[0].set_title('Training Set')
# axs[0].set_xlabel('Momentum (GeV/c)')
# axs[0].set_ylabel('Cherenkov Ange (radians)')
# axs[0].legend()

# # Test set
# axs[1].scatter(X_test_momentum[y_test[:, 0] == 1], theta_test_pion, label='Pion', alpha=0.6)
# axs[1].scatter(X_test_momentum[y_test[:, 1] == 1], theta_test_kaon, label='Kaon', alpha=0.6)
# axs[1].scatter(X_test_momentum[y_test[:, 2] == 1], theta_test_proton, label='Proton', alpha=0.6)
# axs[1].set_title('Test Set')
# axs[1].set_xlabel('Momentum (GeV/c)')
# axs[1].set_ylabel('Cherenkov Ange (radians)')
# axs[1].legend()

# plt.show()


##### Ckov angle histogram for a certain momentum range

# def filter_data(particle_mask, momentum_data, mass, refractive_index_data):
#     # Flatten the 2D arrays to 1D arrays
#     momentum_data = momentum_data.flatten()
#     refractive_index_data = refractive_index_data.flatten()

#     print(f"Particle mask shape: {particle_mask.shape}")  # Should be (10720,)
#     print(f"Momentum data shape: {momentum_data.shape}")  # Should be (10720,)

#     mask_momentum = (momentum_data >= 2.8) & (momentum_data <= 3)
#     print(f"Mask momentum shape: {mask_momentum.shape}")  # Should be (10720,)

#     mask_final = particle_mask & mask_momentum
#     print(f"Mask final shape: {mask_final.shape}")  # Should be (10720,)

#     theta_filtered = calculate_theta(momentum_data[mask_final], mass, refractive_index_data[mask_final])
#     return theta_filtered



# # Filtering data for each particle type for training set
# theta_train_pion_filtered = filter_data(y_train[:, 0] == 1, X_train_momentum, mass_pion, X_train_refractive_index)
# theta_train_kaon_filtered = filter_data(y_train[:, 1] == 1, X_train_momentum, mass_kaon, X_train_refractive_index)
# theta_train_proton_filtered = filter_data(y_train[:, 2] == 1, X_train_momentum, mass_proton, X_train_refractive_index)

# # Filtering data for each particle type for test set
# theta_test_pion_filtered = filter_data(y_test[:, 0] == 1, X_test_momentum, mass_pion, X_test_refractive_index)
# theta_test_kaon_filtered = filter_data(y_test[:, 1] == 1, X_test_momentum, mass_kaon, X_test_refractive_index)
# theta_test_proton_filtered = filter_data(y_test[:, 2] == 1, X_test_momentum, mass_proton, X_test_refractive_index)

# Create a 3x2 grid for the histograms

# Histograms for the training set
# from scipy.stats import norm

# # Function to plot histogram and Gaussian fit
# def plot_hist_and_fit(ax, data, color, title):
#     # Plot histogram
#     n, bins, patches = ax.hist(data, bins=200, alpha=0.6, color=color, label=title)

#     # Fit a Gaussian
#     mu, std = norm.fit(data)
#     xmin, xmax = plt.xlim()
#     x = np.linspace(xmin, xmax, 100)
#     p = norm.pdf(x, mu, std)
#     ax.plot(x, p * max(n), 'k', linewidth=2, color='orange')

#     # Plot vertical lines and add text for 2 std and 3 std
#     for num_std in [2, 3]:
#         ax.axvline(mu - num_std * std, color='purple', linestyle='--')
#         ax.axvline(mu + num_std * std, color='purple', linestyle='--')
#         ax.text(mu + num_std * std, 0.8 * max(n), f'{num_std} std', color='purple')

#     # Display standard deviation
#     ax.text(0.6, 0.8, f'Std: {std:.2f}', transform=ax.transAxes, color='black')
#     ax.set_title(title)


# # Create 3x2 subplots
# fig1, axs1 = plt.subplots(2, 3, figsize=(18, 12))
# fig1.suptitle('In momentum range 2.8 to 3 GeV/c')

# # Create 2x1 subplots for combined histograms
# fig2, axs2 = plt.subplots(2, 1, figsize=(12, 12))

# def plot_histogram_with_fit(ax, data, label, color):
#     # Plot histogram
#     n, bins, patches = ax.hist(data, bins=200, alpha=0.6, label=label, color=color)

#     # Fit a Gaussian
#     mu, std = norm.fit(data)
#     xmin, xmax = plt.xlim()
#     x = np.linspace(xmin, xmax, 100)
#     p = norm.pdf(x, mu, std)
#     ax.plot(x, p * max(n), 'k', linewidth=2, color='orange')

#     # Plot vertical lines and add text for 2 std and 3 std
#     for num_std in [2, 3]:
#         ax.axvline(mu - num_std * std, color='purple', linestyle='--')
#         ax.axvline(mu + num_std * std, color='purple', linestyle='--')
#         ax.text(mu + num_std * std, 0.8 * max(n), f'{num_std} std', color='purple')

#     # Display standard deviation
#     ax.text(0.6, 0.8, f'Std: {std:.4f}', transform=ax.transAxes, color='black')

# # Individual particles for training set
# # plot the fit and dta in the same
# # for these it should be 2x1; in the top all the hist for train, on bottom all for the test
# plot_histogram_with_fit(axs1[0, 0], theta_train_proton_filtered, 'Proton', 'b')
# plot_histogram_with_fit(axs1[0, 1], theta_train_kaon_filtered, 'Kaon', 'g')
# plot_histogram_with_fit(axs1[0, 2], theta_train_pion_filtered, 'Pion', 'r')

# # Individual particles for test set
# plot_histogram_with_fit(axs1[1, 0], theta_test_proton_filtered, 'Proton', 'b')
# plot_histogram_with_fit(axs1[1, 1], theta_test_kaon_filtered, 'Kaon', 'g')
# plot_histogram_with_fit(axs1[1, 2], theta_test_pion_filtered, 'Pion', 'r')
# # plot the fit and dta in the same

# # Combined histogram for training and test set
# plot_histogram_with_fit(axs2[0], theta_train_proton_filtered, 'Proton', 'b')
# plot_histogram_with_fit(axs2[0], theta_train_kaon_filtered, 'Kaon', 'g')
# plot_histogram_with_fit(axs2[0], theta_train_pion_filtered, 'Pion', 'r')

# plot_histogram_with_fit(axs2[1], theta_test_proton_filtered, 'Proton', 'b')
# plot_histogram_with_fit(axs2[1], theta_test_kaon_filtered, 'Kaon', 'g')
# plot_histogram_with_fit(axs2[1], theta_test_pion_filtered, 'Pion', 'r')

# # Show the plots
# plt.show()
# fig, axs = plt.subplots(2, 3, figsize=(18, 12))



# # add the liens for 2 and 3 std dev also here : and indicate the values by prtining on figuer :
# # Histograms for the training set with Gaussian fits
# plot_hist_and_fit(axs[0, 0], theta_train_proton_filtered, 'b', 'Proton (Training Set)')
# plot_hist_and_fit(axs[0, 1], theta_train_kaon_filtered, 'g', 'Kaon (Training Set)')
# plot_hist_and_fit(axs[0, 2], theta_train_pion_filtered, 'r', 'Pion (Training Set)')

# # Histograms for the test set with Gaussian fits
# plot_hist_and_fit(axs[1, 0], theta_test_proton_filtered, 'b', 'Proton (Test Set)')
# plot_hist_and_fit(axs[1, 1], theta_test_kaon_filtered, 'g', 'Kaon (Test Set)')
# plot_hist_and_fit(axs[1, 2], theta_test_pion_filtered, 'r', 'Pion (Test Set)') # add text to show +- 2 adnd 3 std-dev

# # Labels and title
# for i in range(2):
#     for j in range(3):
#         axs[i, j].set_xlabel('Theta (radians)')
#         axs[i, j].set_ylabel('Frequency')
#         axs[i, j].legend()

# plt.tight_layout()
# plt.show()



# #def plot_worst_(model, y_test, X_test_map, X_test_momentum, X_test_refractive_index, X_test_ckov, X_test_mip_position, y_pred):
# def plot_maps(filled_bins_array=None, map_array=None, mip_position_array=None, X_momentum=None, X_refractive_index=None, X_ckov=None, percentage_to_plot=5, resolution = 10):
#   #  print("Shape of y_pred: ", y_pred.shape)
#   # 1. Predict labels on validation data
#   #plot_worst(model, y_test, X_test["X_test_map"], X_test["X_test_momentum"], X_test["X_test_refractive_index"], X_test["X_test_ckov"], X_test["X_test_mip_position"], y_pred_test)

#   # 2. Calculate the difference between predicted and actual labels
#   losses = tf.keras.losses.categorical_crossentropy(y_test, y_pred).numpy()

#   # Sort the indices of the losses from highest to lowest
#   sorted_indices = np.argsort(losses)[::-1]

#   # Get the indices of the worst performing 10%
#   worst_10_percent_indices = sorted_indices[:int(0.1*len(sorted_indices))]

#   # Create figure and axes
#   num_plots = len(worst_10_percent_indices)
#   #fig, axes = plt.subplots(num_plots, 1, figsize=(8, 20))
#   fig, axes = plt.subplots(num_plots,figsize=(8, 20))

#   # Define mass categories
#   mass_categories = ["pion", "kaon", "proton"]

#   # 3. Create plots for these cases, including their feature information and predicted vs actual labels
#   for i, index in enumerate(worst_10_percent_indices):
#       # Get the map and corresponding information
#       map_data = map_array[index, :, :]
#       actual_mass_category = mass_categories[np.argmax(y_test[index])]

#       print(f"y_test[index] = {y_test[index]}")

#       predicted_mass_category = mass_categories[np.argmax(y_pred[index])]
#       ckov = X_test_ckov[index]
#       mip_position = X_test_mip_position[index]
#       momentum = X_test_momentum[index]
#       refractive_index = X_test_refractive_index[index]

#       mass_actual = momentum * np.sqrt(refractive_index**2 * np.cos(ckov)*np.cos(ckov) - 1)

#       # Check if the value is NaN (invalid Cherenkov angle)
#       if np.isnan(mass_actual):
#           mass_actual = "Invalid"

#       # Plot the map
#       axes[i].imshow(map_data, cmap='gray')

#       # Add a red dot at the MIP position
#       axes[i].plot(mip_position[0]*resolution, mip_position[1]*resolution, 'ro')

#       # Set the title with the information
#       axes[i].set_title(f"Actual Mass")#: {actual_mass_category}, Predicted Mass: {predicted_mass_category},\nMass: {mass_actual}, Mass_prob = {y_pred[index]} \nCKOV: {ckov}, MIP Position: {mip_position}, \nMomentum: {momentum}, Refractive Index: {refractive_index}")
#       #
#       axes[i].set_title(f"Actual Mass: {actual_mass_category}, Predicted Mass: {predicted_mass_category},\nMass: {mass_actual}, Mass_prob = {y_pred[index]}, MIP Position: {mip_position}, \nMomentum: {momentum}, Refractive Index: {refractive_index}")

#       #axes[i].set_title(f"Actual Mass: {actual_mass_category}, Predicted Mass: {predicted_mass_category}, Mass: {mass_actual}\nCKOV: {ckov}, MIP Position: {mip_position}, Momentum: {momentum}, Refractive Index: {refractive_index}")
#       axes[i].axis('off')

#       print("\n")
#       print(f"  Actual Mass: {actual_mass_category}, Predicted Mass: {predicted_mass_category},\n Mass: {mass_actual}, Mass_prob = {y_pred[index]} , MIP Position: {mip_position}, \n  Momentum: {momentum}, Refractive Index: {refractive_index}")
#   # Adjust the spacing between subplots
#   plt.tight_layout()

#   # Show the plot
#   plt.show()



def create_lr_scheduler(num_epochs=10, warmup_epochs=50):


    # endre warmup_epochs

    start_lr = 0 # Starting learning rate is 0 for the warm-up
    peak_lr = 0.001/(5*1.2) # was /5 # The learning rate we reach at the end of the warm-up
    end_lr = 5e-6   # The final learning rate at the end of training

    # Calculate decay rate based on peak and end learning rate
    # Notice we adjust the total number of epochs by subtracting the warm-up period
    exp_decay = -np.log(end_lr / peak_lr) / (num_epochs - warmup_epochs)

    # Warm-up followed by exponential decay
    def schedule(epoch):
        if epoch < warmup_epochs:
            return start_lr + ((peak_lr - start_lr) / (warmup_epochs - 1)) * epoch
        else:
            return peak_lr * np.exp(-exp_decay * (epoch - warmup_epochs))

    return tf.keras.callbacks.LearningRateScheduler(schedule)




def plot_lr(num_epochs = 10, history = None):
  div = num_epochs/4
  lrs = 1e-4 * (10 ** (np.arange(num_epochs)/div))
  plt.figure(figsize=(10, 7))
  plt.semilogx(lrs, history.history["loss"]) # we want the x-axis (learning rate) to be log scale
  plt.xlabel("Learning Rate")
  plt.ylabel("Loss")


  plt.title("Learning rate vs. loss");



#@title Default title text

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

def extract_neighborhood_map(candidate_positions, mip_positions, neighborhood_size, map_size):
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
    neighborhood_maps[np.arange(num_samples)[:, np.newaxis], map_indices[:, :, 1], map_indices[:, :, 0]] = mask

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




def plot_worst(model, y_test, X_test_map, X_test_momentum, X_test_refractive_index, X_test_ckov, X_test_mip_position, y_pred):
  # 1. Predict labels on validation data

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
      axes[i].plot(mip_position[0], mip_position[1], 'ro')

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
