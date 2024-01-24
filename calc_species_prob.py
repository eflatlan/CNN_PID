import numpy as np
def calc_species_prob(theta_c_hyps, sigma_ring_padded, theta_cer_padded, mip_mask):


		verbose = False

		if verbose:
			print(f"shape of mip_mask : {mip_mask.shape}")

		"""
		Calculate the species probability per photon per track, and also compute the
		likelihood per track and across all tracks.

		Args:
		theta_c_hyps (np.array): Cherenkov angles per species, shape (n_tracks_max, 5).
		sigma_ring_padded (np.array): Standard deviations, shape (max_length, 1).
		theta_cer_padded (np.array): Cherenkov angles, shape (max_length, 1).
		mip_mask : mask of size geq 3 AND charge geq 200



		Returns:
		tuple: A tuple containing:
				- p_specie_per_track: The probability of each species per photon per track,
															shape (n_tracks_max, max_length, 5).
				- L_track: The likelihood of each photon per track, shape (n_tracks_max, max_length).
				- L_all_tracks: The likelihood of each photon across all tracks, shape (max_length,).
		"""

		#momentums = np.tile(momentums, (1, 5))

		# Expand the dimensions of theta_c_hyps to broadcast across the max_length dimension
		MAX_TRACKS_IN_EVENT = 4
		#theta_c_hyps_expanded = np.repeat(theta_c_hyps[:, np.newaxis, :], repeats=theta_cer_padded.shape[0], axis=1)  # (9, 550, 5)


		theta_c_hyps_newdim = theta_c_hyps[:, np.newaxis, :]

		# Calculate the padding needed for the first dimension
		padding_needed = max(MAX_TRACKS_IN_EVENT - theta_c_hyps_newdim.shape[0], 0)
		#print("theta_c_hyps_newdim:", theta_c_hyps_newdim.shape)  # Should be (MAX_TRACKS_IN_EVENT, 1, 5) if padding was needed		verbose = True

		# Pad the first dimension only
		pad_width = [(padding_needed, 0), (0, 0), (0, 0)]
		theta_c_hyps_expanded = np.pad(theta_c_hyps_newdim, pad_width=pad_width, mode='constant')

		# Check the new shape
		if verbose:
			print("New shape:", theta_c_hyps_expanded.shape)  # Should be (MAX_TRACKS_IN_EVENT, 1, 5) if padding was needed		verbose = True
			print(f"shape of theta_c_hyps_expanded : {theta_c_hyps_expanded.shape}")


		mip_mask = mip_mask.T
		mip_mask =  np.repeat(mip_mask[:, :, np.newaxis], repeats=theta_c_hyps.shape[1], axis=2) # shape will be (1, max_length, 1)


		if verbose:
			print(f"shape of mip_mask : {mip_mask.shape}")


		sigma_ring_padded = sigma_ring_padded.T

		theta_cer_padded = theta_cer_padded.T

		if verbose:
			print(f"shape of theta_cer_padded : {theta_cer_padded.shape}")


		## 		masses_shaped =  np.tile(masses, (9, 1))

		# Ensure sigma_ring_padded and theta_cer_padded are expanded to have a last dimension of 1 for broadcasting
		sigma_ring_padded_expanded =  np.repeat(sigma_ring_padded[:, :, np.newaxis], repeats=theta_c_hyps.shape[1], axis=2) # shape will be (1, max_length, 1)
		theta_cer_padded_expanded =  np.repeat(theta_cer_padded[:, :, np.newaxis], repeats=theta_c_hyps.shape[1], axis=2) # shape will be (1, max_length, 1)

		if verbose:
			print(f"shape of theta_cer_padded_expanded : {theta_cer_padded_expanded.shape}")
			print(f"shape of sigma_ring_padded_expanded : {sigma_ring_padded_expanded.shape}")

		# Calculate the z-score
		num_tracks = sigma_ring_padded_expanded.shape[0]

		mip_mask =  np.repeat(mip_mask[:, :, :], repeats=num_tracks, axis=0) # shape will be (1, max_length, 1)

		if verbose:
			print(f"shape of mip_mask : {mip_mask.shape}")


		theta_cer_padded_expanded =  np.repeat(theta_cer_padded[:, :, np.newaxis], repeats=theta_c_hyps.shape[1], axis=2) # shape will be (1, max_length, 1)


		#		# The shape of z will be (n_tracks, max_length, 5) after broadcasting

		#z_track_number[<number>, :] = # take sum of all other tracks
		# sum of [n_tracks, :, 5] # so sum over the 5 dimensions and the number of tracks !
		# also count the numbers
		#print(f"theta_c_hyps_expanded: {theta_c_hyps_expanded[:5,:10,:5]}")

		sigma_mask = (sigma_ring_padded_expanded > 0) & (sigma_ring_padded_expanded < 10)
		theta_mask = theta_cer_padded_expanded >= 0


		sigma_ring_padded_expanded[~sigma_mask] = 0
		theta_cer_padded_expanded[~theta_mask] = 0

		# Now z_sum_over_tracks_and_species will have the shape (max_length,)
		#print(f"theta_c_hyps_expanded: {theta_c_hyps_expanded[:5,:,:5]}")


		sigma_ring_padded_expanded = np.clip(sigma_ring_padded_expanded, 0, 0.02) #
		z = (theta_c_hyps_expanded[:num_tracks,:,:] - theta_cer_padded_expanded) / sigma_ring_padded_expanded  # Broadcasting




		# only accept +- 2 std-devs
		z_mask = (z < 2) & (z > -2)
		z[mip_mask == True] = 0

		z[z_mask == False] = 0



		valid_counts_per_photon = np.count_nonzero(z, axis=(0, 2))

		n_tracks = z.shape[0]
		for n in range(n_tracks):
			z_for_track = z[n, ...]
			#print(f"t_num{n} of {n_tracks} :  nonzero nonzero_count_vector {np.count_nonzero(z_for_track, axis = 1)}")
			if verbose : 
				print(f"t_num{n} of {n_tracks} : nonzero nonzero_count_vector {np.count_nonzero(z_for_track, axis = 0)}")


		if verbose:
			print(f"shape of theta_c_hyps_expanded: {theta_c_hyps_expanded.shape}")
			print(f"shape of theta_cer_padded_expanded: {theta_cer_padded_expanded.shape}")
			print(f"shape of sigma_ring_padded_expanded: {sigma_ring_padded_expanded.shape}")

			print(f"shape of z-score array: {z.shape}")
			print(f"SHAPE  Valid counts per photon: {valid_counts_per_photon.shape}")

		#print(f"Valid counts per photon:  {valid_counts_per_photon[:5]}")

		#HM = sum_i  1/(1/z_i) # harmomic mean
		# probability :
		# w_i = e^{-z_i} /(Sum_(i=1...n)e^{-z_j}) # weight

		from scipy.stats import norm
		# p-value from error-function of z-score

		#print(f"mip_mask array: {mip_mask[:5,:10,:5]}")

		#print(f"z-score array: {z[:5,:10,:5]}")

		# find related x-score med andre tracks
		p_value = norm.pdf(z)



		#print(f"nonzero p_value: {np.count_nonzero(p_value)}")

		nonzero_count_vector = np.count_nonzero(p_value, axis=(0, 1))

		# Reshaping to 5x1 if necessary
		nonzero_count_vector_5x1 = nonzero_count_vector.reshape(-1, 1)
		#print(f"nonzero p_value: {nonzero_count_vector_5x1.flatten()}")

		p_value[z_mask == False] = 0



		nonzero_count_vector = np.count_nonzero(p_value, axis=(0, 1))

		# Reshaping to 5x1 if necessary
		nonzero_count_vector_5x1 = nonzero_count_vector.reshape(-1, 1)


		# sett p verdier til MIP til 0
		p_value[mip_mask == True] = 0

		nonzero_count_vector = np.count_nonzero(p_value, axis=(0, 1))

		# Reshaping to 5x1 if necessary
		nonzero_count_vector_5x1 = nonzero_count_vector.reshape(-1, 1)






		if verbose:
			print(f"nonzero p_value[z_mask]: {nonzero_count_vector_5x1.flatten()} | nonzero p_value[z_mask]: {np.count_nonzero(p_value)}")

			print(f"nonzero p_value[mip_mask]: {nonzero_count_vector_5x1.flatten()} | nonzero p_value[mip_mask]: {np.count_nonzero(p_value)}")

			print(f"shape of z-score array: {z.shape}")
			# print(f"SHAPE  Valid counts per photon: {valid_counts_per_photon.shape}")

		#print(f"p_value-score array: {p_value[:5,:,:5]}")
		# Normalize these likelihoods to get the probability of each species within a track (p_specie)
		# (num_tracks, num_photons, num_species) -> (num_tracks, num_photons, num_species)
		p_specie_per_track_norm = p_value / np.sum(p_value, axis=2, keepdims=True)
		#print(f"p_specie_per_track-score array: {p_specie_per_track[:5,:10,:5]}")
		p_specie_per_track = p_value

		# L_track : likelihood : sum of species-likelihood per track
		L_track = np.sum(p_value, axis=2)  # Shape: (num_tracks, num_photons)
		#print(f"L_track array: {L_track[:10,:5]}")

		# Sum L_track across all tracks for each photon:
		# This represents the total likelihood of each photon being part of any track
		L_all_tracks = np.sum(L_track, axis=0)  # Shape: (num_photons,)

		return p_specie_per_track_norm, p_specie_per_track, L_track, L_all_tracks, z_mask

		# def	calc_species_prob(theta_c_hyps, sigma_ring_padded, theta_cer_padded)
		# 		"""
		# 				# max_length, num_tracks_in_event, 5 : probability per photon per track

		# 		"""
		# 		print(f"shape theta_c_hyps {theta_c_hyps.shape} |  sigma_ring_padded {sigma_ring_padded.shape} } |  theta_cer_padded {theta_cer_padded.shape} ")
		# 		z = (theta_c_hyps - theta_cer_padded) / sigma_ring_padded
