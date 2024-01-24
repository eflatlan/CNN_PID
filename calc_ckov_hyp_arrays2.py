
import numpy as np
# Particle masses in GeV/c^2
MASS_ELECTRON = 0.000511
MASS_MUON = 0.1057
MASS_PION = 0.1396
MASS_KAON = 0.4937
MASS_PROTON = 0.9383
masses = np.array([MASS_ELECTRON, MASS_MUON, MASS_PION, MASS_KAON, MASS_PROTON])
masses = masses.reshape(1, 5)

NON_VALID_PION_COUNT, NON_VALID_KAON_COUNT, NON_VALID_PROTON_COUNT = 0, 0, 0
GLOBAL_VALID_PION_COUNT, GLOBAL_VALID_KAON_COUNT, GLOBAL_VALID_PROTON_COUNT = 0, 0, 0


def calc_ckov_hyp_arrays(event_data_dict,):
		"""
		# Calculate the ckov hyps
		# shape : num_tracks_in_event, 5 || : 5 species
		theta_c_hyps = calc_ckov_hyp_arrays(momentums, ref_indexes)
		"""
		# fiels in event_data_dict
		#		'Momentum': np.zeros((max_tracks_in_event, 1)),
		#		'RefractiveIndex': np.zeros((max_tracks_in_event, 1)),
		# cos_theta_c = np.sqrt(p^2 + m^2)/(p*n)
		verbose = False

		if verbose :
			print(f"calc_ckov_hyp_arrays")


		try:
			CluCharge  = event_data_dict['CluCharge']#  # Flattening to 1D if necessary
		except Exception as e:
			print(f"calc_ckov_hyp_arrays CluCharge failed w {e}")
		if verbose:
			print(f" CluCharge  {CluCharge.flatten() } ")

		pdg  = event_data_dict['TrackPdg']#.flatten()  # Flattening to 1D if necessary

		if verbose:
			print(f" pdg  {pdg.flatten()  } ")

		refractive_indices = event_data_dict['RefractiveIndex']#.flatten()  # Flattening to 1D if necessary

		momentums = event_data_dict['Momentum']#.flatten()  # Flattening to 1D if necessary
		verbose = False
		if verbose:
			print(f" momentums {momentums.flatten()} ")
			print(f" refractive_indices {refractive_indices.flatten()} ")


		if verbose:
			print(f" shape momentums {momentums.shape} | refractive_indices {refractive_indices.shape}")
			print(f" shape momentums {momentums.shape} | refractive_indices {refractive_indices.shape}")



		momentums = np.tile(momentums, (1, 5))
		refractive_indices = np.tile(refractive_indices, (1, 5))

		if verbose:
			print(f" shape momentums {momentums.shape} | refractive_indices {refractive_indices.shape}")

		# to see ckov radiation per specie
		# if p > p_lim for a given specie : set mask to 0

		#print(f" refractive_indices {refractive_indices} ")

		p_lim = masses/np.sqrt(refractive_indices**2 - 1)



		# denne må endres!
		masses_shaped =  np.tile(masses, (MAX_TRACKS_IN_EVENT, 1))


		if verbose:
			print(f"momentums shape {momentums.shape}")

			print(f"masses_shaped shape {masses_shaped.shape}")

		#print(f"shape p_lim {p_lim.shape} | shape masses_shaped {masses_shaped.shape}")


		p_mask = momentums > p_lim
		#print(f"shape p_lim {p_lim.shape} | shape p_mask {p_mask.shape} | shape momentums {momentums.shape}")
		# the mask to check for radiation is now if p_mask is 1, if its 1 its radiation


		p_squared = momentums ** 2
		p_squared_plus_m_squared = p_squared + masses_shaped ** 2  # Broadcasting to shape (len(momentums), 5)
		#print(f"shape p_squared {p_squared.shape} | shape p_squared_plus_m_squared {p_squared_plus_m_squared.shape} ")

		# Calculate cos(theta_c) for all hypotheses
		cos_theta_c_hyps = np.sqrt(p_squared_plus_m_squared) / (momentums * refractive_indices)
		#print(f"shape p_squared {p_squared.shape} | shape cos_theta_c_hyps {cos_theta_c_hyps.shape} ")

		# Calculate theta_c from cos(theta_c)
		# Note: Where cos(theta_c) > 1 due to numerical issues, it implies no Cherenkov radiation is possible
		#
		# So we can clip the values to 1 to avoid NaNs from the arccos function

		theta_c_hyps = np.zeros_like(cos_theta_c_hyps)
		theta_c_hyps[p_mask] = np.arccos(cos_theta_c_hyps[p_mask])
		# print(f" momentums {momentums} ")
		# print(f" masses_shaped {masses_shaped} ")

		# print(f" theta_c_hyps {theta_c_hyps} ")
		# print(f" cos_theta_c_hyps {cos_theta_c_hyps} ")
		# print(f" p_mask {p_mask} ")

		return theta_c_hyps
