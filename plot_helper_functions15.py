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


def plot_hist(X_train=None, X_test=None, description=None):
    def freedman_diaconis_bins(data):
        data = np.asarray(data)
        data_range = np.nanmax(data) - np.nanmin(data)
        iqr = np.subtract(*np.percentile(data, [75, 25]))
        bin_width = 2 * iqr * (len(data) ** -0.33)
        return int(data_range / bin_width) if bin_width > 0 else 1

    def check_nan_inf(arr):
        return np.isnan(arr).any(), np.isinf(arr).any()

    for label, variable in {**X_train, **X_test}.items():
        has_nan, has_inf = check_nan_inf(variable)
        if has_nan or has_inf:
            warnings.warn(f"{label} contains NaN or Inf values. This may lead to issues.")

    fig, axs = plt.subplots(2, 5, figsize=(25, 10))
    fig.suptitle(f"Training and Testing Data: {description}", fontsize=20)
    fig1, axs1 = plt.subplots(2, 3, figsize=(18, 10))
    fig1.suptitle(f"2D Maps and Projections: {description}", fontsize=20)
    fig2, axs2 = plt.subplots(2, 4, figsize=(16, 10))
    fig2.suptitle(f"Pion Candidate Distributions: {description}", fontsize=20)

    def plot_routine(variables, row_idx):
        axs_idx = 0
        for label, variable in variables.items():
            bins = freedman_diaconis_bins(variable)
            range_val = (np.nanmin(variable), np.nanmax(variable))
            
            if label == 'MIP Position' or label == 'Rad Position':
                variable = np.asarray(variable).reshape(-1, 2)
                mask = (variable[:, 0] != 0) & (variable[:, 1] != 0)
                variable = variable[mask]
                axs1[row_idx, 0].scatter(variable[:, 0], variable[:, 1], marker='o', s=10)
                axs1[row_idx, 0].set_title(f"{'Train' if row_idx == 0 else 'Test'} 2D Map: {label}")
                axs1[row_idx, 1].hist(variable[:, 0], bins=bins, range=range_val, edgecolor='black')
                axs1[row_idx, 1].set_title(f"{label} X")
                axs1[row_idx, 2].hist(variable[:, 1], bins=bins, range=range_val, edgecolor='black')
                axs1[row_idx, 2].set_title(f"{label} Y")
            elif label == 'Pion Candidates':
                variable = np.asarray(variable)
                mask = np.any(variable != 0, axis=-1)
                filtered_variable = variable[mask]
                variable = np.asarray(variable)
                ranges = [(0,125),(0,125),(0,11),(0,200)]
                bins = [10, 0, 1, 20]
                for i in range(4):
                    axs2[row_idx, i].hist(variable[:, :, i][variable[:, :, i] > 0], bins=bins[i], range=ranges[i], edgecolor='black')
                    titles = ['Pion X', 'Pion Y', 'Pion Size', 'Pion Charge']
                    axs2[row_idx, i].set_title(f"{'Train' if row_idx == 0 else 'Test'} {titles[i]}")
                    
            elif label == 'mCluSize':
                axs[row_idx, axs_idx].hist(variable, bins=13, range=(0, 12), edgecolor='black')
                axs[row_idx, axs_idx].set_title(f"{'Train' if row_idx == 0 else 'Test'} {label}")

            elif label == 'Momentum':
                axs[row_idx, axs_idx].hist(variable, bins=50, range=(0, 5), edgecolor='black')
                axs[row_idx, axs_idx].set_title(f"{'Train' if row_idx == 0 else 'Test'} {label}")
            else:
                axs[row_idx, axs_idx].hist(variable, bins=bins, range=range_val, edgecolor='black')
                axs[row_idx, axs_idx].set_title(f"{'Train' if row_idx == 0 else 'Test'} {label}")
                axs_idx += 1

    X_train_variables = {
        'Refractive Index': X_train["X_train_refractive_index"],
        'Momentum': X_train["X_train_momentum"],
        'Phi': X_train["X_train_phi"],
        'Theta': X_train["X_train_theta"],
        'mCluSize': X_train["X_train_mCluSize"],
        'MIP Position': X_train["X_train_mip_position"],
        'Rad Position': X_train["X_train_rad_position"],
        'Pion Candidates': X_train["X_train_pion_candidates"],
    }

    X_test_variables = {
        'Refractive Index': X_test["X_test_refractive_index"],
        'Momentum': X_test["X_test_momentum"],
        'Phi': X_test["X_test_phi"],
        'Theta': X_test["X_test_theta"],
        'mCluSize': X_test["X_test_mCluSize"],
        'MIP Position': X_test["X_test_mip_position"],
        'Rad Position': X_test["X_test_rad_position"],
        'Pion Candidates': X_test["X_test_pion_candidates"],
    }

    plot_routine(X_train_variables, 0)
    plot_routine(X_test_variables, 1)

    plt.show()
