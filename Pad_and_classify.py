from __future__ import print_function


from tensorflow.keras.preprocessing.sequence import pad_sequences
import numpy as np

from tensorflow.keras.preprocessing.sequence import pad_sequences
import numpy as np

from tensorflow.keras.preprocessing.sequence import pad_sequences
import numpy as np

import numpy as np

def pad_and_stack(sequences, max_length=None):
    # Your existing code
    try:
        # Try padding, if max_length is not None, pad or truncate to that length
        padded_sequences = pad_sequences(sequences, maxlen=max_length, padding='post', dtype='float32')  # Changed dtype to 'floar32'
    except ValueError:
        # Fallback: manually pad with zeros
        max_len = max_length if max_length is not None else max(len(seq) for seq in sequences)
        padded_sequences = np.array([np.pad(seq, (0, max_len - len(seq)), 'constant', constant_values=0) for seq in sequences])

    return padded_sequences





def classify_candidates_with_pad_sequences(x_values_data, y_values_data, q_values_data, mCluSize_lista, candStatus_values_data, max_length_nested, xmip_list, ymip_list):


    # Pad the sequences
    print(f"x_values_data Type of sequences: {type(x_values_data)}, shape of sequences: {np.shape(x_values_data)}")
    print(f"Type of sequences: {type(y_values_data)}, shape of sequences: {np.shape(y_values_data)}")
    print(f"Type of sequences: {type(mCluSize_lista)}, shape of sequences: {np.shape(mCluSize_lista)}")
    print(f"Type of sequences: {type(q_values_data)}, shape of sequences: {np.shape(q_values_data)}")
    print(f"Type of sequences: {type(candStatus_values_data)}, shape of sequences: {np.shape(candStatus_values_data)}")

    print(type(x_values_data), np.shape(x_values_data))
    if np.ndim(x_values_data) == 1:
        x_values_data = np.expand_dims(x_values_data, axis=-1)



    # xmip, ymip):
     # xmip, ymip):

     
    #xmip_list = pad_and_stack(xmip_list, max_length=max_length_nested)
    #ymip_list = pad_and_stack(ymip_list, max_length=max_length_nested)

    x_padded = pad_and_stack(x_values_data, max_length=max_length_nested)
    y_padded = pad_and_stack(y_values_data, max_length=max_length_nested)
    #chi2_padded = pad_and_stack(chi2_values_data, max_length=max_length_nested)
    q_padded = pad_and_stack(q_values_data, max_length=max_length_nested)
    size_padded = pad_and_stack(mCluSize_lista, max_length=max_length_nested)

    #xe_padded = pad_and_stack(xe_values_data, max_length=max_length_nested)
    #ye_padded = pad_and_stack(ye_values_data, max_length=max_length_nested)
    candStatus_padded = pad_and_stack(candStatus_values_data, max_length=max_length_nested)
    candStatus_padded = candStatus_padded.astype(int)
    # Stack the data into a single array
    #padded_data = np.stack([x_padded, y_padded, chi2_padded, q_padded, xe_padded, ye_padded, candStatus_padded], axis=-1)
    padded_data = np.stack([x_padded, y_padded, q_padded, size_padded], axis=-1)

    # Create masks for different particle types based on 'candStatus'
    pion_mask = (candStatus_padded & 4).astype(bool)
    kaon_mask = (candStatus_padded & 2).astype(bool)
    proton_mask = (candStatus_padded & 1).astype(bool)

    # Create masks that are True for all elements
    pion_mask = np.ones_like(candStatus_padded, dtype=bool)
    kaon_mask = np.ones_like(candStatus_padded, dtype=bool)
    proton_mask = np.ones_like(candStatus_padded, dtype=bool)
    non_mask = np.ones_like(candStatus_padded, dtype=bool)

    # Initialize arrays for particle candidates
    pion_candidates = np.zeros_like(padded_data)
    kaon_candidates = np.zeros_like(padded_data)
    proton_candidates = np.zeros_like(padded_data)
    non_candidates = np.zeros_like(padded_data)

    print(f"shape pion_candidates : {np.shape(pion_candidates)}")
    print(f"shape pion_mask: {np.shape(pion_mask)}")
    print(f"shape padded_data : {np.shape(padded_data)}")
    print(f"shape candStatus_padded : {np.shape(candStatus_padded)}")



    # Initialize arrays for particle candidates
    shape_data = np.shape(padded_data)
    pion_candidates = np.zeros(shape_data)
    kaon_candidates = np.zeros(shape_data)
    proton_candidates = np.zeros(shape_data)
    non_candidates = np.zeros(shape_data)



    x_np = np.array(xmip_list)
    y_np = np.array(ymip_list)
    diffx = np.array(x_padded - x_np)
    diffy = np.array(y_padded - y_np)

    dist = np.sqrt(diffx*diffx + diffy*diffy)
    dist_to_mip_mask = dist > 5

    # Print the first 5 elements of the 174-length dimension for the first entry in the 7-length dimension
    print("First 5 of x_np for the first entry:", x_np[0, :5])
    print("First 5 of padded_data[:,:,0] for the first entry:", padded_data[0, :5, 0])
    print("First 5 of diffx for the first entry:", diffx[0, :5])

    print("First 5 of y_np for the first entry:", y_np[0, :5])
    print("First 5 of padded_data[:,:,1] for the first entry:", padded_data[0, :5, 1])
    print("First 5 of diffy for the first entry:", diffy[0, :5])

    print("First 5 of dist for the first entry:", dist[0, :5])
    print("First 5 of dist_to_mip_mask for the first entry:", dist_to_mip_mask[0, :5])




    print(f"shape dist_to_mip_mask: {np.shape(dist_to_mip_mask)}")
    print(f"shape dist : {np.shape(dist)}")
    print(f"shape diffx : {np.shape(diffx)}")

    pion_mask = (candStatus_padded & 4).astype(bool) 
    kaon_mask = (candStatus_padded & 2).astype(bool)
    proton_mask = (candStatus_padded & 1).astype(bool)
    non_mask = (candStatus_padded == 0).astype(bool)  # Create a mask for non-candidates

    # Assuming pion_candidates, kaon_candidates are initialized and have the same shape as padded_data
    pion_candidates[pion_mask & dist_to_mip_mask] = padded_data[pion_mask & dist_to_mip_mask]
    kaon_candidates[kaon_mask & dist_to_mip_mask] = padded_data[kaon_mask & dist_to_mip_mask]
    proton_candidates[proton_mask & dist_to_mip_mask] = padded_data[proton_mask & dist_to_mip_mask]


    non_candidates[non_mask] = padded_data[~(pion_mask | kaon_mask | proton_mask)]
    non_candidates[non_mask] = padded_data[non_mask]

    # Return the candidates and the original status data
    return pion_candidates, kaon_candidates, proton_candidates, non_candidates, candStatus_padded
    #return padded_data, padded_data, padded_data, padded_data, padded_data
