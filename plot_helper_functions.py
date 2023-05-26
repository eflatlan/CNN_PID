from __future__ import print_function
import numpy as np



from sklearn.metrics import precision_recall_curve, confusion_matrix
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
from itertools import cycle
import matplotlib.pyplot as plt
from numpy.linalg import norm
import os
import h5py
import tensorflow as tf
def plot_training_history(history=None, vector_of_weights=None, vector_of_weights2=None, dropout=None, y_pred_train=None, y_pred_test=None, y_train_true=None, y_test_true=None, relu_alpha = None):
    try:       
        fig, axs = plt.subplots(3, 2, figsize=(16, 18))
        # Plot Training and Validation Loss
        axs[0, 0].plot(history.history["loss"], label="Train Loss")
        axs[0, 0].plot(history.history["val_loss"], label="Validation Loss")
        axs[0, 0].set_xlabel("Epochs")
        axs[0, 0].set_ylabel("Loss")
        axs[0, 0].legend()

        # Plot Training and Validation Accuracy
        axs[0, 1].plot(history.history["accuracy"], label="Train Accuracy")
        axs[0, 1].plot(history.history["val_accuracy"], label="Validation Accuracy")
        axs[0, 1].set_xlabel("Epochs")
        axs[0, 1].set_ylabel("Accuracy")
        axs[0, 1].legend()

        # Add text for vector_of_weights, vector_of_weights2, and dropout
        axs[0, 0].text(0.05, 0.05, f'Weights: {len(vector_of_weights)}', transform=axs[0, 0].transAxes, fontsize=12)
        axs[0, 1].text(0.05, 0.05, f'Weights2: {len(vector_of_weights2)}', transform=axs[0, 1].transAxes, fontsize=12)
        axs[0, 1].text(0.05, 0.95, f'Dropout: {dropout} LR alpha = {relu_alpha}', transform=axs[0, 1].transAxes, fontsize=12)

        # Binarize the output
        y_train_bin = label_binarize(y_train_true, classes=[0, 1, 2])
        y_test_bin = label_binarize(y_test_true, classes=[0, 1, 2])

        # For each class
        for i in range(3):
            precision_train, recall_train, _ = precision_recall_curve(y_train_bin[:, i], y_pred_train[:, i])
            axs[1, 0].plot(recall_train, precision_train, lw=2, label='class {}'.format(i))

        axs[1, 0].set_xlabel("recall")
        axs[1, 0].set_ylabel("precision")
        axs[1, 0].legend(loc="best")
        axs[1, 0].set_title("precision vs. recall curve Train")

        for i in range(3):
            precision, recall, _ = precision_recall_curve(y_test_bin[:, i], y_pred_test[:, i])
            axs[1, 1].plot(recall, precision, lw=2, label='class {}'.format(i))

        axs[1, 1].set_xlabel("recall")
        axs[1, 1].set_ylabel("precision")
        axs[1, 1].legend(loc="best")
        axs[1, 1].set_title("precision vs. recall curve Test")

        # Add text for vector_of_weights, vector_of_weights2, and dropout
        axs[1, 1].text(0.05, 0.05, f'Weights: {len(vector_of_weights)}', transform=axs[1, 1].transAxes, fontsize=12)
        axs[1, 1].text(0.05, 0.95, f'Weights2: {len(vector_of_weights2)}\nDropout: {dropout}\nLR alpha = {relu_alpha}', transform=axs[1, 1].transAxes, fontsize=12)

        # Compute confusion matrices
        cm_train = confusion_matrix(y_train_true.argmax(axis=1), y_pred_train.argmax(axis=1))
        cm_test = confusion_matrix(y_test_true.argmax(axis=1), y_pred_test.argmax(axis=1))

        axs[2, 0].imshow(cm_train, cmap='Blues', interpolation='nearest')
        axs[2, 0].set_xticks(np.arange(2))
        axs[2, 0].set_yticks(np.arange(2))
        axs[2, 0].set_xticklabels(['0', '1'])
        axs[2, 0].set_yticklabels(['0', '1'])
        axs[2, 0].set_xlabel('Predicted label')
        axs[2, 0].set_ylabel('True label')
        axs[2, 0].set_title('Confusion Matrix (Train)')

        axs[2, 1].imshow(cm_test, cmap='Blues', interpolation='nearest')
        axs[2, 1].set_xticks(np.arange(2))
        axs[2, 1].set_yticks(np.arange(2))
        axs[2, 1].set_xticklabels(['0', '1'])
        axs[2, 1].set_yticklabels(['0', '1'])
        axs[2, 1].set_xlabel('Predicted label')
        axs[2, 1].set_ylabel('True label')
        axs[2, 1].set_title('Confusion Matrix (Test)')

        # Add text for vector_of_weights, vector_of_weights2, and dropout
        axs[2, 0].text(0.05, 0.95, f'Weights: {len(vector_of_weights)}\nWeights2: {len(vector_of_weights2)}\nDropout: {dropout}\nLR alpha = {relu_alpha}', transform=axs[2, 0].transAxes, fontsize=12)
        axs[2, 1].text(0.05, 0.95, f'Weights: {len(vector_of_weights)}\nWeights2: {len(vector_of_weights2)}\nDropout: {dropout}\nLR alpha = {relu_alpha}', transform=axs[2, 1].transAxes, fontsize=12)


        plt.tight_layout()
        plt.show()
        print(f'Weights: {vector_of_weights}\nWeights2: {vector_of_weights2}\nDropout: {dropout}\nLR alpha = {relu_alpha}')
    except Exception as e:
        print(f"plot_training_history failed with {e}")

def plot_dist2mip_histograms(X_test_dist2mip = None, resolution = None):
        
    if X_test_dist2mip is None: 
        print(f"plot_dist2mip_histograms got empty X_test_dist2mip")

    if resolution is None: 
        print(f"plot_dist2mip_histograms got empty resolution")

    else:
        try:
            # Define the bin edges for each histogram.
            bins1 = np.arange(0, self.resolution*10.1, 0.01)
            bins2 = np.arange(4, self.resolution*8.1, 0.01)
            bins3 = np.arange(0, self.resolution*2.1, 0.01)
    

            print(f'There are {num_zeros} zeros in X_train_dist2mip.')
            print(f"X_train_dist2mip shape = {np.array(X_train_dist2mip).shape}")
            fig, axs = plt.subplots(1, 3, figsize=(15, 5))

            # Plot the first histogram.
            axs[0].hist(X_test_dist2mip, bins=bins1, edgecolor='black')
            axs[0].set_title('Histogram 1 of X_test_dist2mip')
            axs[0].set_xlabel('Values')
            axs[0].set_ylabel('Frequency')

            # Plot the second histogram.
            axs[1].hist(X_test_dist2mip, bins=bins2, edgecolor='black')
            axs[1].set_title('Histogram 2 of X_test_dist2mip')
            axs[1].set_xlabel('Values')
            axs[1].set_ylabel('Frequency')
    
            # Plot the third histogram.
            axs[2].hist(X_test_dist2mip, bins=np.arange(0, .01, 0.001), edgecolor='black')
            axs[2].set_title('Histogram 3 of X_test_dist2mip')
            axs[2].set_xlabel('Values')
            axs[2].set_ylabel('Frequency')

            # Display the plot.
            plt.tight_layout()
            plt.show()


            # # Define the bin edges for each histogram.
            # bins11 = np.arange(self.resolution*0.5, self.resolution*10.1, 0.1)
            # bins22 = np.arange(self.resolution*0.5, self.resolution*2.1, 0.1)
            # bins33 = np.arange(self.resolution*0.5, self.resolution*5.1, 0.1)

            # # Create histograms without plotting, to identify peaks.
            # hist1, edges1 = np.histogram(X_test_dist2mip, bins=bins11)
            # hist2, edges2 = np.histogram(X_test_dist2mip, bins=bins22)
            # hist3, edges3 = np.histogram(X_test_dist2mip, bins=bins33)

            # # Find peaks (you might want to adjust the parameters).
            # peaks1, _ = find_peaks(hist1, height=10) # adjust the height as per your requirement
            # peaks2, _ = find_peaks(hist2, height=10)
            # peaks3, _ = find_peaks(hist3, height=10)

            # # Create a figure with 3 subplots (one row, three columns).
            # fig, axs = plt.subplots(1, 3, figsize=(15, 5))

            # # Plot the first histogram with peaks.
            # axs[0].hist(X_test_dist2mip, bins=bins11, edgecolor='black')
            # axs[0].plot(edges1[peaks1], hist1[peaks1], "ro") # Peaks marked in red
            # axs[0].set_title('Histogram 1 with Peaks')
            # axs[0].set_xlabel('Values')
            # axs[0].set_ylabel('Frequency')

            # # Plot the second histogram with peaks.
            # axs[1].hist(X_test_dist2mip, bins=bins22, edgecolor='black')
            # axs[1].plot(edges2[peaks2], hist2[peaks2], "ro") # Peaks marked in red
            # axs[1].set_title('Histogram 2 with Peaks')
            # axs[1].set_xlabel('Values')
            # axs[1].set_ylabel('Frequency')

            # # Plot the third histogram with peaks.
            # axs[2].hist(X_test_dist2mip, bins=bins33, edgecolor='black')
            # axs[2].plot(edges3[peaks3], hist3[peaks3], "ro") # Peaks marked in red
            # axs[2].set_title('Histogram 3 with Peaks')
            # axs[2].set_xlabel('Values')
            # axs[2].set_ylabel('Frequency')

            # # Display the plot.
            # plt.tight_layout()
            # plt.show()

            # for i in range(5):
            #   fig, axs = plt.subplots(1, 3, figsize=(15, 5))
            #   temp = X_test_dist2mip[:, i]
            #   max_value = np.max(temp)
            #   print(f"max_value = {max_value}")
            #   print(f"temp shape = {np.array(temp).shape}")
            
            #   # Plot the first histogram.
            #   axs[0].hist(temp, bins=bins1, edgecolor='black')
            #   axs[0].set_title(f'{i} X_test_dist2mip')
            #   axs[0].set_xlabel('Values')
            #   axs[0].set_ylabel('Frequency')

            #   # Plot the second histogram.
            #   axs[1].hist(temp, bins=bins2, edgecolor='black')
            #   axs[1].set_title(f'{i} X_test_dist2mip')
            #   axs[1].set_xlabel('Values')
            #   axs[1].set_ylabel('Frequency')

            #   # Plot the third histogram.
            #   axs[2].hist(temp, bins=bins3, edgecolor='black')
            #   axs[2].set_title(f'{i} X_test_dist2mip')
            #   axs[2].set_xlabel('Values')
            #   axs[2].set_ylabel('Frequency')

            #   # Display the plot.
            #   plt.tight_layout()
            #   plt.show()

            # for i in range(5):
            #   fig, axs = plt.subplots(1, 3, figsize=(15, 5))
            #   temp = X_test_dist2mip[i, :]
            #   max_value = np.max(temp)
            #   print(f"max_value = {max_value}")
            #   print(f"temp shape = {np.array(temp).shape}")

            #   # Plot the first histogram.
            #   axs[0].hist(temp, bins=bins1, edgecolor='black')
            #   axs[0].set_title(f'{i} a X_test_dist2mip')
            #   axs[0].set_xlabel('a Values')
            #   axs[0].set_ylabel('Frequency')

            #   # Plot the second histogram.
            #   axs[1].hist(temp, bins=bins2, edgecolor='black')
            #   axs[1].set_title(f'{i} a X_test_dist2mip')
            #   axs[1].set_xlabel('a Values')
            #   axs[1].set_ylabel('a Frequency')

            #   # Plot the third histogram.
            #   axs[2].hist(temp, bins=bins3*10, edgecolor='black')
            #   axs[2].set_title(f'{i} a X_test_dist2mip')
            #   axs[2].set_xlabel('a Values')
            #   axs[2].set_ylabel('a Frequency')

            #   # Display the plot.
            #   plt.tight_layout()
            #   plt.show()
        except Exception as e:
            print(f"plot_helper_functions.py function plot_dist2mip_histograms failed with error : {e}")
