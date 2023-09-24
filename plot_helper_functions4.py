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

def plot_training_history(history, y_pred_train, y_pred_test, y_train_true, y_test_true,
                          vector_of_weights=None, vector_of_weights2=None, dropout=None, relu_alpha=None):

    fig, axs = plt.subplots(10, 1, figsize=(12, 40))
    
    # Parameters for annotation
    params = []
    if vector_of_weights:
        params.append(f"VoW: {vector_of_weights}")
    if vector_of_weights2:
        params.append(f"VoW2: {vector_of_weights2}")
    if dropout:
        params.append(f"Dropout: {dropout}")
    if relu_alpha:
        params.append(f"ReLU alpha: {relu_alpha}")
    annotation = '\n'.join(params)

    # Add annotation to each axis
    for ax in axs:
        ax.annotate(annotation, (0.05, 0.95), textcoords='axes fraction',
                    bbox=dict(boxstyle="square", fc="white", ec="black"),
                    va="top", ha="left")

    # Overall Loss & Accuracy
    axs[0].plot(history.history['loss'], label="Train")
    axs[0].plot(history.history['val_loss'], label="Validation")
    axs[0].legend()
    axs[0].set_title("Overall Loss")

    axs[1].plot(history.history['accuracy'], label="Train")
    axs[1].plot(history.history['val_accuracy'], label="Validation")
    axs[1].legend()
    axs[1].set_title("Overall Accuracy")

    # P-R Curves (Train & Test)
    y_train_bin = label_binarize(y_train_true, classes=[0, 1, 2])
    y_test_bin = label_binarize(y_test_true, classes=[0, 1, 2])
    for i in range(3):
        precision_train, recall_train, _ = precision_recall_curve(y_train_bin[:, i], y_pred_train[:, i])
        precision_test, recall_test, _ = precision_recall_curve(y_test_bin[:, i], y_pred_test[:, i])
        axs[2].plot(recall_train, precision_train, lw=2, label=f"Class {i} Train")
        axs[3].plot(recall_test, precision_test, lw=2, label=f"Class {i} Test")
    axs[2].legend()
    axs[2].set_title("Precision-Recall Curve (Train)")
    axs[3].legend()
    axs[3].set_title("Precision-Recall Curve (Test)")

    # Confusion Matrices (Train & Test)
    cm_train = confusion_matrix(y_train_true.argmax(axis=1), y_pred_train.argmax(axis=1))
    axs[4].set_title("Training Confusion Matrix")
    plot_confusion_matrix(axs[4], cm_train)

    cm_test = confusion_matrix(y_test_true.argmax(axis=1), y_pred_test.argmax(axis=1))
    axs[5].set_title("Testing Confusion Matrix")
    plot_confusion_matrix(axs[5], cm_test)

    # Loss & Accuracy plots for each species
    for idx in range(3):
        axs[6 + idx].plot(history.history['loss'], label="Train")
        axs[6 + idx].plot(history.history['val_loss'], label="Validation")
        axs[6 + idx].legend()
        axs[6 + idx].set_title(f"Species {idx+1} - Loss")

        axs[7 + idx].plot(history.history['accuracy'], label="Train")
        axs[7 + idx].plot(history.history['val_accuracy'], label="Validation")
        axs[7 + idx].legend()
        axs[7 + idx].set_title(f"Species {idx+1} - Accuracy")

    plt.tight_layout()
    plt.show()
