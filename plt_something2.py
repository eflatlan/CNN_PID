import numpy as np
import matplotlib.pyplot as plt


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


# def calc_ckov_hyp_array_group(group):
# 		"""
# 		# Calculate the ckov hyps
# 		# shape : num_tracks_in_event, 5 || : 5 species
# 		theta_c_hyps = calc_ckov_hyp_arrays(momentums, ref_indexes)
# 		"""
# 		# fiels in event_data_dict
# 		#		'Momentum': np.zeros((max_tracks_in_event, 1)),
# 		#		'RefractiveIndex': np.zeros((max_tracks_in_event, 1)),
# 		# cos_theta_c = np.sqrt(p^2 + m^2)/(p*n)


# 		abs_pdg = group.attrs["TrackPdg"]
# 		momentum = group.attrs["Momentum"]
# 		PhiP = group.attrs["PhiP"]
# 		ThetaP = group.attrs["ThetaP"]
# 		xMip = group.attrs["xMip"]
# 		yMip = group.attrs["yMip"]


# 		if verbose:
# 			print(f" pdg  {pdg.flatten()  } ")


# 		verbose = True
# 		if verbose:
# 			print(f" momentums {momentum.flatten()} ")
# 			print(f" refractive_indices {refractive_indices.flatten()} ")


# 		if verbose:
# 			print(f" shape momentums {momentum.shape} | refractive_indices {refractive_indices.shape}")
# 			print(f" shape momentums {momentum.shape} | refractive_indices {refractive_indices.shape}")



# 		momentum = np.tile(momentum, (1, 5))
# 		refractive_indices = np.tile(refractive_indices, (1, 5))

# 		if verbose:
# 			print(f" shape momentums {momentums.shape} | refractive_indices {refractive_indices.shape}")

# 		# to see ckov radiation per specie
# 		# if p > p_lim for a given specie : set mask to 0

# 		#print(f" refractive_indices {refractive_indices} ")

# 		p_lim = masses/np.sqrt(refractive_indices**2 - 1)



# 		# denne må endres!
# 		print(f"momentums shape {momentums.shape}")
# 		masses_shaped =  np.tile(masses, (MAX_TRACKS_IN_EVENT, 1))
# 		print(f"masses_shaped shape {masses_shaped.shape}")

# 		#print(f"shape p_lim {p_lim.shape} | shape masses_shaped {masses_shaped.shape}")


# 		p_mask = momentums > p_lim
# 		#print(f"shape p_lim {p_lim.shape} | shape p_mask {p_mask.shape} | shape momentums {momentums.shape}")
# 		# the mask to check for radiation is now if p_mask is 1, if its 1 its radiation


# 		p_squared = momentums ** 2
# 		p_squared_plus_m_squared = p_squared + masses_shaped ** 2  # Broadcasting to shape (len(momentums), 5)
# 		#print(f"shape p_squared {p_squared.shape} | shape p_squared_plus_m_squared {p_squared_plus_m_squared.shape} ")

# 		# Calculate cos(theta_c) for all hypotheses
# 		cos_theta_c_hyps = np.sqrt(p_squared_plus_m_squared) / (momentums * refractive_indices)
# 		#print(f"shape p_squared {p_squared.shape} | shape cos_theta_c_hyps {cos_theta_c_hyps.shape} ")

# 		# Calculate theta_c from cos(theta_c)
# 		# Note: Where cos(theta_c) > 1 due to numerical issues, it implies no Cherenkov radiation is possible
# 		#
# 		# So we can clip the values to 1 to avoid NaNs from the arccos function

# 		theta_c_hyps = np.zeros_like(cos_theta_c_hyps)
# 		theta_c_hyps[p_mask] = np.arccos(cos_theta_c_hyps[p_mask])
# 		# print(f" momentums {momentums} ")
# 		# print(f" masses_shaped {masses_shaped} ")

# 		# print(f" theta_c_hyps {theta_c_hyps} ")
# 		# print(f" cos_theta_c_hyps {cos_theta_c_hyps} ")
# 		# print(f" p_mask {p_mask} ")

# 		return theta_c_hyps


MASS_ELECTRON = 0.000511
MASS_MUON = 0.1057
MASS_PION = 0.1396
MASS_KAON = 0.4937
MASS_PROTON = 0.9383
masses = np.array([MASS_ELECTRON, MASS_MUON, MASS_PION, MASS_KAON, MASS_PROTON])
masses = masses.reshape(1, 5)

NON_VALID_PION_COUNT, NON_VALID_KAON_COUNT, NON_VALID_PROTON_COUNT = 0, 0, 0
GLOBAL_VALID_PION_COUNT, GLOBAL_VALID_KAON_COUNT, GLOBAL_VALID_PROTON_COUNT = 0, 0, 0

verbose = True

def calc_ckov_hyp_array_group(group):
		"""
		# Calculate the ckov hyps
		# shape : num_tracks_in_event, 5 || : 5 species
		theta_c_hyps = calc_ckov_hyp_arrays(momentums, ref_indexes)
		"""
		# fiels in event_data_dict
		#		'Momentum': np.zeros((max_tracks_in_event, 1)),
		#		'RefractiveIndex': np.zeros((max_tracks_in_event, 1)),
		# cos_theta_c = np.sqrt(p^2 + m^2)/(p*n)


		pdg = group.attrs["TrackPdg"]
		momentum = group.attrs["Momentum"]
		PhiP = group.attrs["PhiP"]
		ThetaP = group.attrs["ThetaP"]
		xMip = group.attrs["xMip"]
		yMip = group.attrs["yMip"]

		refractive_indices = group.attrs["RefractiveIndex"]

		if True:
			print(f" pdg  {pdg.flatten()  } ")


		verbose = True
		if True:
			print(f" momentums {momentum.flatten()} ")
			print(f" refractive_indices {refractive_indices.flatten()} ")


		if True:
			print(f" shape momentums {momentum.shape} | refractive_indices {refractive_indices.shape}")
			print(f" shape momentums {momentum.shape} | refractive_indices {refractive_indices.shape}")


		momentum = np.tile(momentum, (1, 5))
		refractive_indices = np.tile(refractive_indices, (1, 5))

		if True:
			print(f" shape momentums {refractive_indices.shape} | refractive_indices {refractive_indices.shape}")

		# to see ckov radiation per specie
		# if p > p_lim for a given specie : set mask to 0

		#print(f" refractive_indices {refractive_indices} ")

		p_lim = masses/np.sqrt(refractive_indices**2 - 1)



		# denne må endres!
		print(f"momentums shape {momentum.shape}")
		masses_shaped =  np.tile(masses, (1, 1))
		print(f"masses_shaped shape {masses_shaped.shape}")

		#print(f"shape p_lim {p_lim.shape} | shape masses_shaped {masses_shaped.shape}")


		p_mask = momentum > p_lim
		#print(f"shape p_lim {p_lim.shape} | shape p_mask {p_mask.shape} | shape momentums {momentums.shape}")
		# the mask to check for radiation is now if p_mask is 1, if its 1 its radiation


		p_squared = momentum ** 2
		p_squared_plus_m_squared = p_squared + masses_shaped ** 2  # Broadcasting to shape (len(momentums), 5)
		#print(f"shape p_squared {p_squared.shape} | shape p_squared_plus_m_squared {p_squared_plus_m_squared.shape} ")

		# Calculate cos(theta_c) for all hypotheses
		cos_theta_c_hyps = np.sqrt(p_squared_plus_m_squared) / (momentum * refractive_indices)
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

def plt_something(group):
    verbose = True
    print(f"called plot_something")
    pdg = group.attrs["TrackPdg"]
    if np.abs(pdg) in [211, 321, 2212]:
        momentum = group.attrs["Momentum"]
        PhiP = group.attrs["PhiP"]
        ThetaP = group.attrs["ThetaP"]
        xMip = group.attrs["xMip"]
        yMip = group.attrs["yMip"]
        ckov_hyps_cp = np.copy(calc_ckov_hyp_array_group(group)[:, 2:]) # theta_c_hyps should be shape ((max_tracks_in_event, 3))
        print(f"ckov_hyps_cp shape {ckov_hyps_cp.shape}")
        #keys = ["x_values", "y_values", "sigmaRingValues", "thetaCerValues"]
        species_names = ['Pions', 'Kaons', 'Protons']


        fig, axs = plt.subplots(4, 3, figsize=(12, 12))
        # plt.subplots_adjust(wspace=0.2)  # Adjust the horizontal spacing between subplots


        # plt2  = plt.figure(figsize=(18, 6))  # Adjust the figure size as needed for horizontal subplots
        # plt2.subplots_adjust(wspace=0.2)  # Adjust the horizontal spacing between subplots


        x = np.copy(np.asarray(group["x_values"][...]))
        y = np.copy(np.asarray(group["y_values"][...]))
        sigma = np.copy(np.asarray(group["sigmaRingValues"][...]))
        theta_cer = np.copy(np.asarray(group["thetaCerValues"][...]))
        phi = np.copy(np.asarray(group["phiCerValues"][...]))


        for i, specie in enumerate(species_names):
            print(f"\n\n\n ckov_hyps_cp {ckov_hyps_cp[0, :]}")
		# # 		padded_data = np.stack([x_padded, y_padded, q_padded, size_padded, phi_cer_padded, theta_cer_padded, sigma_ring_padded, pion_prob_per_specie, kaon_prob_per_specie, proton_prob_per_specie, L_track, L_all_tracks], axis=-1)
            if species_names[i] == "Pions":
                ckov = ckov_hyps_cp[0 , 0]

            elif species_names[i] == "Kaons":
                ckov = ckov_hyps_cp[0, 1]


            elif species_names[i] == "Protons":
                ckov = ckov_hyps_cp[0, 2]

            print(f"specie {species_names[i]} | ckov {ckov}")



            sigma[sigma < 0] = 0

            sigma[sigma > 0.025] = 0.025
            from scipy.stats import norm
            theta_cer[theta_cer < 0] = 0

            n_std =  np.asarray((theta_cer - ckov)/sigma)
            z_mask = np.abs(n_std) < 2

            # print(f"z_mask shape {z_mask.shape}")
            p = norm.pdf(n_std)
            print(f"mask new cond {p[z_mask]}")
		    # print(f"shapes p {p.shape} | x {x.shape}")
		    # print(f"x : {x[z_mask]} | y : {y[z_mask]} | p : {p[z_mask]} | theta_cer : {theta_cer[z_mask]}")

		    # # Plot with mask

            try:
                axs[0, i].scatter(x[z_mask], y[z_mask], c=p[z_mask], cmap='viridis', vmin=0)
                axs[0, i].scatter(xMip, yMip, color='red', marker='o', s=10)  # 'o' marker is a red dot, s is marker size
                axs[0, i].set_title(f'Probability per photon for {specie} (With Mask)')
                axs[0, i].set_xlabel('X')
                axs[0, i].set_ylabel('Y')
                axs[0, i].set_xlim(0, 135)
                axs[0, i].set_ylim(0, 135)
            except Exception as e :
                print(f"failed gettin figure with {e}")



            try:
                #Plot without mask
                axs[1, i].scatter(x, y, c=p, cmap='viridis', vmin=0)
                axs[1, i].set_title(f'Probability per photon for {specie} (Without Mask)')
                axs[1, i].set_xlabel('X')
                axs[1, i].set_ylabel('Y')
                axs[1, i].set_xlim(0, 135)
                axs[1, i].set_ylim(0, 135)
                axs[1, i].scatter(xMip, yMip, color='red', marker='o', s=10)  # 'o' marker is a red dot, s is marker size
            except Exception as e :
                print(f"failed gettin figure with {e}")



            import matplotlib.colors as colors

		    # # Plot with mask
            #min_t = theta_cer[z_mask].min()
            #max_t = theta_cer[z_mask].max()


            try:
                sc = axs[2, i].scatter(x[z_mask], y[z_mask], c=theta_cer[z_mask], cmap='viridis', norm=colors.LogNorm(vmin=0.2, vmax=1))
                axs[2, i].scatter(xMip, yMip, color='red', marker='o', s=10)  # 'o' marker is a red dot, s is marker size
                axs[2, i].set_title(f'Theta-Cherenkov per photon for {specie} (With Mask)')
                axs[2, i].set_xlabel('X')
                axs[2, i].set_ylabel('Y')
                axs[2, i].set_xlim(0, 135)
                axs[2, i].set_ylim(0, 135)
                plt.colorbar(sc, ax=axs[2, i])
            except Exception as e :
                print(f"failed gettin figure with {e}")


            try:
                #Plot without mask
                sc=axs[3, i].scatter(x, y, c=theta_cer, cmap='viridis', norm=colors.LogNorm(vmin=0.2, vmax=1))
                axs[3, i].set_title(f'Theta-Cherenkov per {specie} (Without Mask)')
                axs[3, i].set_xlabel('X')
                axs[3, i].set_ylabel('Y')
                axs[3, i].set_xlim(0, 135)
                axs[3, i].set_ylim(0, 135)
                axs[3, i].scatter(xMip, yMip, color='red', marker='o', s=10)  # 'o' marker is a red dot, s is marker size
                plt.colorbar(sc, ax=axs[3, i])
            except Exception as e :
                print(f"failed gettin figure with {e}")

            #plt.show()

            print(f"z_mask shape {z_mask.shape}")
            #p = norm.pdf(n_std)
            #print(f"mask new cond {p[z_mask]}")
            #print(f"shapes p {p.shape} | x {x.shape}")
            print(f"x : {x[z_mask]} | y : {y[z_mask]} | p : {p[z_mask]} | theta_cer : {theta_cer[z_mask]}")

        plt.show()


    plt.show()

