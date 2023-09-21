from __future__ import print_function

from sklearn.metrics import precision_recall_curve
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
from itertools import cycle

import sys

print(sys.getrecursionlimit()) # Prints 1000

print_vals = False
from numpy.linalg import norm
from tensorflow.keras.backend import expand_dims
from tensorflow.keras.preprocessing.sequence import pad_sequences
from sklearn.metrics import precision_recall_curve, confusion_matrix

from scipy.signal import find_peaks

import os
import h5py
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Activation, Input, Conv2D, Lambda, Flatten, Dense, concatenate, BatchNormalization, MaxPooling2D, Dropout, LeakyReLU, Masking, Embedding

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelBinarizer
import matplotlib.pyplot as plt
from tensorflow.keras import regularizers

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
import matplotlib.pyplot as plt


from tensorflow.keras.callbacks import EarlyStopping

# create a callback
early_stopping = EarlyStopping(
    monitor='val_loss', # you can monitor 'val_loss' or 'val_accuracy'
    patience=100, # stop training if the monitored quantity does not improve for 50 epochs
    restore_best_weights=True, # restore model weights from the epoch with the best value
)

class Constants:
    PION_MASS = 0.1396
    KAON_MASS = 0.4937
    PROTON_MASS = 0.938

np.set_printoptions(precision=4)

@staticmethod
def calculate_mass(momentum, refractiveIndex, ckov):
    """ args : momentum, refractiveIndex, ckov
        returns : mass
    """
    mass = momentum * np.sqrt((refractiveIndex * np.cos(ckov))**2 - 1)
    return mass


class ParticleDataUtils:

    def __init__(self, filenames =  [], percentage_to_read = 100):
        self.filename = filenames
        self.percentage_to_read = percentage_to_read
        self.particle_vector = []
        self.load_data(filenames=filenames)

        #
        # jeg flyttet
        self.particle_info = self.process_data(self.particle_vector, self.percentage_to_read)
        self.num_particles = len(self.particle_info)


        # new scalers to be created :
        self.phi_scaler, self.phi_stats = self.create_scalar_scaler("phiP")
        self.theta_scaler, self.theta_stats = self.create_scalar_scaler("thetaP")
        self.refractive_index_scaler, self.refractive_index_stats = self.create_scalar_scaler("refractiveIndex")
        self.momentum_scaler, self.momentum_stats = self.create_scalar_scaler("momentum")
        self.mCluSize_scaler, self.mCluSize_stats = self.create_scalar_scaler("mCluSize")
        self.mCluCharge_scaler, self.mCluCharge_stats = self.create_scalar_scaler("mCluCharge")
        print("Created scaler for scalars")

        # 2D
        self.mip_scaler, self.mip_stats = self.create_2D_scaler("mip_position")
        self.rad_scaler, self.rad_stats = self.create_2D_scaler("rad_position")

        # # vector of 2D
        self.proton_scalers, self.proton_stats = self.create_vector_scaler("proton_candidates")
        self.kaon_scalers, self.kaon_stats = self.create_vector_scaler("kaon_candidates")
        self.pion_scalers, self.pion_stats = self.create_vector_scaler("pion_candidates")



        #self.ckov_scaler, self.ckov_stats = self.create_scaler("ckov")
        #self.distances_scaler, self.distances_stats = self.create_scaler("distances")def classify_candidates(candidates_data):


    class Candidate2:
        def __init__(self, x_values, y_values, chi2_values, q_values, xe_values, ye_values, candStatus_values):
            self.x_values = x_values
            self.y_values = y_values
            self.chi2_values = chi2_values
            self.q_values = q_values
            self.xe_values = xe_values
            self.ye_values = ye_values
            self.candStatus_values = candStatus_values


    class ParticleInfo: # p
        def __init__(self,  momentum, refractiveIndex, xRad, yRad, xMIP, yMIP, thetaP, phiP, mCluCharge, mCluSize, non_candidates, pion_candidates, kaon_candidates, proton_candidates, mTrackPdg):
            self.momentum = momentum # this dhould be with
            self.refractiveIndex = refractiveIndex # with
            self.xRad = xRad # with
            self.yRad = yRad # with

            self.xMIP = xMIP # with
            self.yMIP = yMIP # with

            self.mCluCharge = mCluCharge # with
            self.mCluSize = mCluSize # with

            self.thetaP = thetaP# with
            self.phiP = phiP# with
            self.non_candidates = non_candidates # with the field candStatus is a int that is 0..7, please make it categorical
            self.rad_position = [xRad, yRad]
            self.mip_position = [xMIP, yMIP]


            self.pion_candidates = pion_candidates # pion_candidates = [1 if (int(candStatus) & 4) == 4 else 0 for candStatus in candsCombined]
            self.kaon_candidates = kaon_candidates # = kaon_candidates [1 if (int(candStatus) & 2) == 2 else 0 for candStatus in candsCombined]
            self.proton_candidates = proton_candidates # = proton_candidates[1 if (int(candStatus) & 1) == 1 else 0 for candStatus in candsCombined]


            self.mTrackPdg = mTrackPdg # with

            abs_mTrackPdg = abs(self.mTrackPdg)  # Take the absolute value

            # Set particleType based on absolute PDG code
            if abs_mTrackPdg == 211:
                self.particleType = 'pion'
            elif abs_mTrackPdg == 321:
                self.particleType = 'kaon'
            elif abs_mTrackPdg == 2212:
                self.particleType = 'proton'
            else:
                self.particleType = 'other'
                #print(f"pdg type was other : {abs_mTrackPdg}")


        @staticmethod
        def infer_mass_category_from_ckov(momentum, refractiveIndex, ckov):
            mass = momentum * np.sqrt((refractiveIndex * np.cos(ckov))**2 - 1)

            mass_category = "unknown"
            if abs(mass - Constants.PION_MASS) < 1e-4:
                mass_category = "pion"
            elif abs(mass - Constants.KAON_MASS) < 1e-4:
                mass_category = "kaon"
            elif abs(mass - Constants.PROTON_MASS) < 1e-4:
                mass_category = "proton"
            if print_vals:
              print(f"\ninfer_mass_category_from_ckov :  momentum = {momentum}|  mass_calc = {mass} |  mass_category={mass_category} | refractiveIndex = {refractiveIndex} | ckov = {ckov}")
            return mass_category

        @staticmethod
        def infer_mass_category(mass):
            if abs(mass - Constants.PION_MASS) < 1e-6:
                return "pion"
            elif abs(mass - Constants.KAON_MASS) < 1e-6:
                return "kaon"
            elif abs(mass - Constants.PROTON_MASS) < 1e-6:
                return "proton"
            else:
                return "unknown"

        def __str__(self):
            if print_vals:
              return (f"ParticleInfo(momentum={self.momentum} | mass={self.mass} |  mass_category={self.mass_category} | "
                      f"refractiveIndex={self.refractiveIndex} | ckov={self.ckov} | rad_position={len(self.rad_position)}, "
                      f"mip_position={self.mip_position})")

    #  ''' def calculate_distances_to_mip(self):
    #       """Calculate Euclidean distances from all filled bins to MIP position"""
    #       filledBins_np = np.array(self.filledBins)
    #       mip_position_np = np.array(self.mip_position)

    #       distances = np.linalg.norm(filledBins_np - mip_position_np, axis=1)
    #       return distances'''


    def load_data(self, filenames):
      drive_path = '/content/drive/MyDrive/Colab Notebooks/CERN_ML/CNN_PID/'

      max_length_nested = 0
      num_particles = 0
      file_num = 0
      for filename in filenames:
        file_path = os.path.join(drive_path, filename)
        print(f"Reading file {file_num} : {filename}")
        with h5py.File(file_path, 'r') as file:
            for i, group_name in enumerate(file):
                group = file[group_name]
                num_particles = num_particles + 1
                current_length = len(group['ye_values'][...])
                if current_length > max_length_nested:
                    max_length_nested = current_length
                    print(f" i {i} max_length_nested {max_length_nested}")


      # Lists to store scalar and array-like attributes
      momentum_list = np.zeros((num_particles, 1))
      refractiveIndex_list = np.zeros((num_particles, 1))
      xRad_list = np.zeros((num_particles, 1))
      yRad_list = np.zeros((num_particles, 1))
      xMIP_list = np.zeros((num_particles, 1))
      yMIP_list = np.zeros((num_particles, 1))
      thetaP_list = np.zeros((num_particles, 1))
      phiP_list = np.zeros((num_particles, 1))
      mCluCharge_list = np.zeros((num_particles, 1))
      mCluSize_list = np.zeros((num_particles, 1))
      mTrackPdg_list = np.zeros((num_particles, 1))


      x_values_data_list = np.zeros((num_particles, max_length_nested))
      y_values_data_list =  np.zeros((num_particles, max_length_nested))
      chi2_values_data_list =  np.zeros((num_particles, max_length_nested))
      q_values_data_list =  np.zeros((num_particles, max_length_nested))
      xe_values_data_list =  np.zeros((num_particles, max_length_nested))
      ye_values_data_list =  np.zeros((num_particles, max_length_nested))
      candStatus_values_data_list =  np.zeros((num_particles, max_length_nested))
      size_clu_lst =  np.zeros((num_particles, max_length_nested))


      index_particle = 0
      for filename in filenames:
        file_path = os.path.join(drive_path, filename)
        print(f"Reading file {file_num} : {filename}")
        with h5py.File(file_path, 'r') as file:
            for i, group_name in enumerate(file):
                group = file[group_name]

                # Store scalar attributes into lists
                momentum_list[index_particle] = group.attrs['Momentum']
                refractiveIndex_list[index_particle] = group.attrs['RefractiveIndex']
                xRad_list[index_particle] = group.attrs['xRad']
                yRad_list[index_particle] = group.attrs['yRad']
                xMIP_list[index_particle] = group.attrs['xMip']
                yMIP_list[index_particle] = group.attrs['yMip']
                thetaP_list[index_particle] = group.attrs['ThetaP']
                phiP_list[index_particle] = group.attrs['PhiP']
                mCluCharge_list[index_particle] = group.attrs['CluCharge']
                mCluSize_list[index_particle] = group.attrs['CluSize']
                mTrackPdg_list[index_particle] = group.attrs['TrackPdg']

                actual_length = len(group['x_values'][...])
                x_values_data_list[index_particle, :actual_length] = group['x_values'][...]

                actual_length = len(group['y_values'][...])
                y_values_data_list[index_particle, :actual_length] = group['y_values'][...]

                actual_length = len(group['chi2_values'][...])
                chi2_values_data_list[index_particle, :actual_length] = group['chi2_values'][...]

                actual_length = len(group['q_values'][...])
                q_values_data_list[index_particle, :actual_length] = group['q_values'][...]

                actual_length = len(group['xe_values'][...])
                xe_values_data_list[index_particle, :actual_length] = group['xe_values'][...]

                actual_length = len(group['ye_values'][...])
                ye_values_data_list[index_particle, :actual_length] = group['ye_values'][...]

                actual_length = len(group['candStatus_values'][...])
                candStatus_values_data_list[index_particle, :actual_length] = group['candStatus_values'][...]


                actual_length = len(group['mSize_values'][...])
                size_clu_lst[index_particle, :actual_length] = group['mSize_values'][...]

                index_particle += 1

                # cehck the length of current entry inn   ye_values_data_list and update max_length_nested if logner
                current_length = len(group['ye_values'][...])
                if current_length > max_length_nested:
                    max_length_nested = current_length
                    print(f" i {i} max_length_nested {max_length_nested}")

      # Perform your vectorized operations here
      # ...

      # Create ParticleInfo objects

      # Assuming that classify_candidates_with_pad_sequences can work on lists
      # of data arrays, and it pads them to a uniform size
      pion_candidates, kaon_candidates, proton_candidates, non_candidates, cand_combined = classify_candidates_with_pad_sequences(
          x_values_data_list, y_values_data_list, size_clu_lst,
          q_values_data_list, candStatus_values_data_list,
          max_length_nested,xMIP_list,yMIP_list

      )
      print(f"pion_candidates shape {pion_candidates.shape}")

      print(f"Dtype : {pion_candidates.dtype}")  # Output will be something like: int64

      particle_vector = [None] * len(momentum_list)

      MIP_list = np.hstack([xMIP_list, yMIP_list])

      # Reshape the array to (N, 1, 2)
      MIP_list_reshaped = MIP_list[:, np.newaxis, :]

      # Extract only the x and y coordinates from pion_candidates
      pion_candidates_xy = pion_candidates[:, :, :2]

      # Compute squared differences
      diff = np.sum((pion_candidates_xy - MIP_list_reshaped)**2, axis=2)
      r_max = 35
      r_min = 0.5
      # Count the number of points within radius r (5 to 40)
      non_padded_mask = np.any(pion_candidates[:, :, :2] != 0, axis=2)
      within_r_mask = (diff >= r_min**2) & (diff <= r_max**2) & non_padded_mask

      count_within_r = np.sum(within_r_mask, axis=1)


      for i in range(len(momentum_list)):
          particle_info = ParticleDataUtils.ParticleInfo(
              momentum_list[i], refractiveIndex_list[i], xRad_list[i], yRad_list[i],
              xMIP_list[i], yMIP_list[i], thetaP_list[i], phiP_list[i],
              mCluCharge_list[i], mCluSize_list[i], non_candidates[i], pion_candidates[i],
              kaon_candidates[i], proton_candidates[i], mTrackPdg_list[i]
          )
          abs_pdg = abs(mTrackPdg_list[i])

          if i == 0:
            print(f"pion_candidates[i] shape {pion_candidates[i].shape}")
            print(f"Dtype : {pion_candidates[i].dtype}")  # Output will be something like: int64


          #if mTrackPdg_list[i] in ['pion', 'kaon', 'proton']:
          if abs_pdg in [211, 321, 2212]:
            cnt_min = 3
            cnt_max = 50
            particle_vector[i] = particle_info
            self.particle_vector.append(particle_info)
            # if count_within_r[i] > cnt_min and count_within_r[i] < cnt_max:

            #   if mCluCharge_list[i] == 200:
            #     a =1 #print(f"mCluCharge_list[i] {mCluCharge_list[i]}")
            #   else :
            #     #print(f"mCluCharge_list[i] {mCluCharge_list[i]}")
            #     particle_vector[i] = particle_info
            #     self.particle_vector.append(particle_info)
            # else:
            #   print(f"Number withing ", count_within_r[i])


          #else :
          #  print(f"PDG Type: {mTrackPdg_list[i]}")

      print(f"Slenght particle_vector {len(self.particle_vector)}")



    def process_data(self, particle_vector, percentage):

        # Calculate the number of particles based on the percentage
        num_particles = int(len(self.particle_vector) * (percentage / 100.0))

        # Slice the particle_vector to the desired percentage
        particle_vector = self.particle_vector[:num_particles]
        return particle_vector

    def create_scalar_scaler(self, feature):
        try:
            if feature in ["momentum", "refractiveIndex", "phiP", "thetaP", "mCluSize", "mCluCharge"]:  # added mCluSize and mCluCharge
                values = np.array([getattr(info, feature) for info in self.particle_info]).reshape(-1, 1)
                if values.size == 0:
                    raise ValueError(f"Empty values array for feature: {feature}")
                scaler = StandardScaler()
                scaled_values = scaler.fit_transform(values)
                values = scaled_values

                stats = {
                    "mean": scaler.mean_[0],
                    "std": scaler.scale_[0]
                }
                return scaler, stats
            else:
                raise ValueError(f"Invalid feature: {feature}")

        except Exception as e:
            print(f"An error occurred in create_scalar_scaler: {e}")
            raise



    def create_2D_scaler(self, feature):
        try:
            if feature in ["mip_position", "rad_position"]:
                values = np.array([getattr(info, feature) for info in self.particle_info])
                if values.size == 0:
                    raise ValueError(f"Empty values array for feature: {feature}")

                scaler_x = StandardScaler()
                scaler_y = StandardScaler()



                print(f"create_2D_scalervalues.shape : {values.shape}")

                n_samples, dim, r = values.shape
                print(f"create_2D_scalervalues.shape : {values[:, 0].shape}")

                feature_values_x = values[:, 0].reshape(n_samples, 1)
                scaled_values_x = scaler_x.fit_transform(feature_values_x).reshape(n_samples, 1)
                values[:, 0] = scaled_values_x

                feature_values_y = values[:, 1].reshape(n_samples, 1)
                scaled_values_y = scaler_y.fit_transform(feature_values_y).reshape(n_samples, 1)
                values[:, 1] = scaled_values_y


                # Store both the scalers
                scalers = {'x': scaler_x, 'y': scaler_y}
                stats = {
                    'x': {
                        "mean": scaler_x.mean_[0],
                        "std": scaler_x.scale_[0]
                    },
                    'y': {
                        "mean": scaler_y.mean_[0],
                        "std": scaler_y.scale_[0]
                    }
                }
                return scalers, stats
        except Exception as e:
            print(f"An error occurred in create_2D_scaler: {e}")
            raise

    def create_vector_scaler(self, feature):
        try:
            if feature in ["pion_candidates", "kaon_candidates", "proton_candidates"]:
                values = np.array([getattr(info, feature) for info in self.particle_info])
                if values.size == 0:
                    raise ValueError(f"Empty values array for feature: {feature}")

                n_samples, n_clusters, n_features = values.shape
                scalers = []

                for i in range(n_features):
                    scaler = StandardScaler()
                    feature_values = values[:, :, i].reshape(-1, 1)
                    scaled_values = scaler.fit_transform(feature_values).reshape(n_samples, n_clusters)
                    values[:, :, i] = scaled_values
                    scalers.append(scaler)
                stats = []
                for scaler in scalers:
                    stat = {
                        "mean": scaler.mean_[0],
                        "std": scaler.scale_[0]
                }
                    stats.append(stat)

                return scalers, stats

        except Exception as e:
            print(f"An error occurred in create_vector_scaler: {e}")
            raise



# TODO : denne skal mulgiens fjernes helt ?
# create a map, the resolution is the "inverse"
def create_map(filledBins=None, resolution=4):
    # Add an offset to your map shape calculation to handle edge cases
    offset = 0
    map_shape = (int(144 * resolution + offset), int(160 * resolution + offset))
    map_data = np.zeros(map_shape, dtype=np.int32)
    if filledBins is not None:
        filledBins_np = np.array(filledBins)
        indices = (filledBins_np * resolution).astype(int)
        #print(f"create_map : indices shape : {np.array(indices).shape}")

        map_data[indices[:, 1], indices[:, 0]] = 1

        ind = np.wherme(map_data == 1)
        #print(f"create_map : ind shape : {np.array(ind).shape}")
    return map_data
