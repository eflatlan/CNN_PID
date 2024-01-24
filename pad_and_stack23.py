




def pad_and_stack2(sequences, max_length=None):
	try:
		# Try padding, if max_length is not None, pad or truncate to that length
		padded_sequences = pad_sequences(sequences, maxlen=max_length, padding='post', dtype='float32')  # Changed dtype to 'floar32'
		#print("padded ok")
	except ValueError:
		# Fallback: manually pad with zeros
		max_len = max_length if max_length is not None else max(len(seq) for seq in sequences)
		padded_sequences = np.array([np.pad(seq, (0, max_len - len(seq)), 'constant', constant_values=0) for seq in sequences])
		#print("revert to other pad")

	return padded_sequences


def pad_and_stack3(sequences, max_length=None):
		# Check if sequences is a single numpy array and not a list of sequences.
		# If it's a single array, wrap it in a list.
		if isinstance(sequences, np.ndarray) and sequences.ndim == 1:
				sequences = [sequences]  # Wrap the single array in a list to create a list of sequences.

		try:
				# Try padding, if max_length is not None, pad or truncate to that length.
				padded_sequences = pad_sequences(sequences, maxlen=max_length, padding='post', dtype='float32')  # dtype was corrected to 'float32'.
				#print("padded ok")
		except ValueError as e:
				print(f"pad_and_stack3 : ValueError: {e}")  # Print the error for debugging.
				# Fallback: manually pad with zeros if there's an error.
				max_len = max_length if max_length is not None else max(len(seq) for seq in sequences)
				padded_sequences = np.array([np.pad(seq, (0, max_len - len(seq)), 'constant', constant_values=0) for seq in sequences])
				#print("reverted to manual padding")

		if padded_sequences.shape[1] == max_length:
				padded_sequences = padded_sequences.transpose()


		return padded_sequences

import numpy as np
import numpy as np
