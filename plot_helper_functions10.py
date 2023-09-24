python
Copy code
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import label_binarize
from sklearn.metrics import precision_recall_curve, confusion_matrix
import warnings

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
                for i in range(4):
                    axs2[row_idx, i].hist(variable[:, :, i][variable[:, :, i] > 0], bins=bins, range=range_val, edgecolor='black')
                    titles = ['Pion X', 'Pion Y', 'Pion Size', 'Pion Charge']
                    axs2[row_idx, i].set_title(f"{'Train' if row_idx == 0 else 'Test'} {titles[i]}")
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

def plot_training_history(history, y_pred_train, y_pred_test, y_train_true, y_test_true,
                          vector_of_weights=None, vector_of_weights2=None, dropout=None, relu_alpha=None):

    fig, axs = plt.subplots(5, 2, figsize=(15, 40))
    
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
    for ax in axs.ravel():
        ax.annotate(annotation, (0.05, 0.95), textcoords='axes fraction',
                    bbox=dict(boxstyle="square", fc="white", ec="black"),
                    va="top", ha="left")

    # Plot training & validation accuracy values
    axs[0, 0].plot(history.history['accuracy'])
    axs[0, 0].plot(history.history['val_accuracy'])
    axs[0, 0].set_title('Model accuracy')
    axs[0, 0].set_ylabel('Accuracy')
    axs[0, 0].set_xlabel('Epoch')
    axs[0, 0].legend(['Train', 'Test'], loc='upper left')

    # Plot training & validation loss values
    axs[0, 1].plot(history.history['loss'])
    axs[0, 1].plot(history.history['val_loss'])
    axs[0, 1].set_title('Model loss')
    axs[0, 1].set_ylabel('Loss')
    axs[0, 1].set_xlabel('Epoch')
    axs[0, 1].legend(['Train', 'Test'], loc='upper left')

    # Assuming binary classification
    precision_train, recall_train, _ = precision_recall_curve(y_train_true, y_pred_train)
    precision_test, recall_test, _ = precision_recall_curve(y_test_true, y_pred_test)
    axs[1, 0].plot(recall_train, precision_train, lw=2, label='Train')
    axs[1, 0].plot(recall_test, precision_test, lw=2, label='Test')
    axs[1, 0].set_xlabel('Recall')
    axs[1, 0].set_ylabel('Precision')
    axs[1, 0].legend()
    axs[1, 0].set_title('Precision-Recall Curve')

    # Confusion matrix for training data
    cm_train = confusion_matrix(y_train_true, np.round(y_pred_train))
    axs[2, 0].imshow(cm_train, interpolation='nearest', cmap=plt.cm.Blues)
    axs[2, 0].set_title('Train Confusion Matrix')
    axs[2, 0].set_ylabel('True label')
    axs[2, 0].set_xlabel('Predicted label')
    
    # Confusion matrix for test data
    cm_test = confusion_matrix(y_test_true, np.round(y_pred_test))
    axs[2, 1].imshow(cm_test, interpolation='nearest', cmap=plt.cm.Blues)
    axs[2, 1].set_title('Test Confusion Matrix')
    axs[2, 1].set_ylabel('True label')
    axs[2, 1].set_xlabel('Predicted label')

    # Further plots can be added based on your requirements

    plt.tight_layout()
    plt.show()
