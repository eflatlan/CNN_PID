from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import label_binarize
from sklearn.metrics import precision_recall_curve, confusion_matrix

def plot_confusion_matrix(ax, cm):
    """Utility function to plot the confusion matrix."""
    ax.imshow(cm, cmap='Blues', interpolation='nearest')
    ax.set_xticks(np.arange(3))
    ax.set_yticks(np.arange(3))
    ax.set_xticklabels(['0', '1', '2'])
    ax.set_yticklabels(['0', '1', '2'])
    for x in range(3):
        for y in range(3):
            percent = cm[x, y] / np.sum(cm[x, :]) * 100  # Percentage formula
            ax.text(y, x, f"{cm[x, y]} ({percent:.1f}%)", ha='center', va='center', color='red')

def plot_training_history(history, y_pred_train, y_pred_test, y_train_true, y_test_true):

    fig, axs = plt.subplots(4, 4, figsize=(25, 20))
    
    # Plots for Overall data
    labels = ['Overall Loss', 'Overall Accuracy', 'Overall P-R Curve', 'Overall Confusion Matrix']

    # Loss & Accuracy
    axs[0, 0].plot(history.history['loss'], label="Train")
    axs[0, 0].plot(history.history['val_loss'], label="Validation")
    axs[0, 0].legend()
    axs[0, 0].set_title(labels[0])

    axs[0, 1].plot(history.history['accuracy'], label="Train")
    axs[0, 1].plot(history.history['val_accuracy'], label="Validation")
    axs[0, 1].legend()
    axs[0, 1].set_title(labels[1])
    
    # P-R Curve for overall data
    y_train_bin = label_binarize(y_train_true, classes=[0, 1, 2])
    y_test_bin = label_binarize(y_test_true, classes=[0, 1, 2])
    for i in range(3):
        precision_train, recall_train, _ = precision_recall_curve(y_train_bin[:, i], y_pred_train[:, i])
        axs[0, 2].plot(recall_train, precision_train, lw=2, label=f"Class {i}")
        axs[0, 2].set_title(labels[2])
        axs[0, 2].legend()
    
    # Confusion Matrix
    cm = confusion_matrix(y_train_true.argmax(axis=1), y_pred_train.argmax(axis=1))
    axs[0, 3].set_title(labels[3])
    plot_confusion_matrix(axs[0, 3], cm)

    # Plots for each class separately
    for idx in range(3):
        axs[1, idx].plot(history.history['loss'], label="Train")
        axs[1, idx].plot(history.history['val_loss'], label="Validation")
        axs[1, idx].legend()
        axs[1, idx].set_title(f"Species {idx+1} - Loss")

        axs[2, idx].plot(history.history['accuracy'], label="Train")
        axs[2, idx].plot(history.history['val_accuracy'], label="Validation")
        axs[2, idx].legend()
        axs[2, idx].set_title(f"Species {idx+1} - Accuracy")

        # P-R curve for each species
        precision_train, recall_train, _ = precision_recall_curve(y_train_bin[:, idx], y_pred_train[:, idx])
        axs[3, idx].plot(recall_train, precision_train, lw=2)
        axs[3, idx].set_title(f"Species {idx+1} - Precision vs Recall")

    plt.tight_layout()
    plt.show()
        
