from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import label_binarize
from sklearn.metrics import precision_recall_curve, confusion_matrix
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import label_binarize
from sklearn.metrics import precision_recall_curve, confusion_matrix

def plot_confusion_matrix(ax, cm, title):
    ax.imshow(cm, cmap='Blues', interpolation='nearest')
    ax.set_xticks(np.arange(3))
    ax.set_yticks(np.arange(3))
    ax.set_xticklabels(['0', '1', '2'])
    ax.set_yticklabels(['0', '1', '2'])
    ax.set_title(title)
    for x in range(3):
        for y in range(3):
            percent = cm[x, y] / np.sum(cm[x, :]) * 100  # Percentage formula
            ax.text(y, x, f"{cm[x, y]} ({percent:.1f}%)", ha='center', va='center', color='red')

def plot_training_history(history, y_pred_train, y_pred_test, y_train_true, y_test_true):

    fig, axs = plt.subplots(5, 2, figsize=(20, 40))

    # Overall Loss & Accuracy
    axs[0, 0].plot(history.history['loss'], label="Train")
    axs[0, 0].plot(history.history['val_loss'], label="Validation")
    axs[0, 0].legend()
    axs[0, 0].set_title("Overall Loss")

    axs[0, 1].plot(history.history['accuracy'], label="Train")
    axs[0, 1].plot(history.history['val_accuracy'], label="Validation")
    axs[0, 1].legend()
    axs[0, 1].set_title("Overall Accuracy")

    # P-R Curves
    y_train_bin = label_binarize(y_train_true, classes=[0, 1, 2])
    y_test_bin = label_binarize(y_test_true, classes=[0, 1, 2])
    for i in range(3):
        precision_train, recall_train, _ = precision_recall_curve(y_train_bin[:, i], y_pred_train[:, i])
        precision_test, recall_test, _ = precision_recall_curve(y_test_bin[:, i], y_pred_test[:, i])
        axs[1, 0].plot(recall_train, precision_train, lw=2, label=f"Class {i} Train")
        axs[1, 1].plot(recall_test, precision_test, lw=2, label=f"Class {i} Test")
    axs[1, 0].legend()
    axs[1, 0].set_title("Precision-Recall Curve (Train)")
    axs[1, 1].legend()
    axs[1, 1].set_title("Precision-Recall Curve (Test)")

    # Confusion Matrices
    cm_train = confusion_matrix(y_train_true.argmax(axis=1), y_pred_train.argmax(axis=1))
    plot_confusion_matrix(axs[2, 0], cm_train, "Training Confusion Matrix")

    cm_test = confusion_matrix(y_test_true.argmax(axis=1), y_pred_test.argmax(axis=1))
    plot_confusion_matrix(axs[2, 1], cm_test, "Testing Confusion Matrix")

    # Loss & Accuracy plots for each species
    for idx in range(3):
        axs[3 + idx, 0].plot(history.history['loss'], label="Train")
        axs[3 + idx, 0].plot(history.history['val_loss'], label="Validation")
        axs[3 + idx, 0].legend()
        axs[3 + idx, 0].set_title(f"Species {idx+1} - Loss")

        axs[3 + idx, 1].plot(history.history['accuracy'], label="Train")
        axs[3 + idx, 1].plot(history.history['val_accuracy'], label="Validation")
        axs[3 + idx, 1].legend()
        axs[3 + idx, 1].set_title(f"Species {idx+1} - Accuracy")

    plt.tight_layout()
    plt.show()

