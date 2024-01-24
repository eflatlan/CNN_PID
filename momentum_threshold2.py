import os
from tensorflow.keras.callbacks import EarlyStopping

# create a callback
early_stopping = EarlyStopping(
    monitor='val_loss', # you can monitor 'val_loss' or 'val_accuracy'
    patience=100, # stop training if the monitored quantity does not improve for 50 epochs
    restore_best_weights=True, # restore model weights from the epoch with the best value
)

from tensorflow.keras.preprocessing.sequence import pad_sequences
import numpy as np
from tensorflow.keras.preprocessing.sequence import pad_sequences
import numpy as np

import numpy as np

import numpy as np
MASS_PION = 0.1396
MASS_KAON = 0.4937
MASS_PROTON = 0.938

# Squared masses
MASS_PION_SQ = MASS_PION * MASS_PION
MASS_KAON_SQ = MASS_KAON * MASS_KAON
MASS_PROTON_SQ = MASS_PROTON * MASS_PROTON
REF_INDEX_FREON = 1.29  # Given refraction index
REF_INDEX_FREON_SQ = REF_INDEX_FREON * REF_INDEX_FREON

def threshold_momentum(pdg_code, p):
	"""
	Calculate the threshold momentum based on the given PDG code.

	:param pdg_code: PDG code of the particle.
	p : momentum of track
	(tbd : refindex)

	:return : boolean value for Cherenkov radiation.
	"""

	# Determine mass based on PDG code
	if abs(pdg_code) == 211:
		mass = MASS_PION
	elif abs(pdg_code) == 321:
		mass = MASS_KAON
	elif abs(pdg_code) == 2212:
		mass = MASS_PROTON
	else:
		raise ValueError(f"Unsupported PDG code: {pdg_code}")

	p_lim = mass/(np.sqrt(REF_INDEX_FREON_SQ-1))
	#print(f" p_lim {p_lim} p {p}")
	return p_lim < p



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


def count_non_zero_charges_vectorized(datasets):
	return [np.count_nonzero(data[:, :, 3], axis=1) for data in datasets]

import numpy as np
import numpy as np
