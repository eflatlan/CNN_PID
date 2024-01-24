
import numpy as np
mass_Pion = 0.1396
mass_Kaon = 0.4937
mass_Proton = 0.938
mass_Pion_sq = mass_Pion * mass_Pion
mass_Kaon_sq = mass_Kaon * mass_Kaon
mass_Proton_sq = mass_Proton * mass_Proton


def calculate_cherenkov_tracks(p, n):
    p_sq = np.power(p, 2)

    # Calculations for min_ckov using n_min
    cos_ckov_denom = p * n
    cos_ckov_Pion = np.sqrt(p_sq + mass_Pion_sq) / cos_ckov_denom
    cos_ckov_Kaon = np.sqrt(p_sq + mass_Kaon_sq) / cos_ckov_denom
    cos_ckov_Proton = np.sqrt(p_sq + mass_Proton_sq) / cos_ckov_denom

    pion_ckov = np.arccos(cos_ckov_Pion)#.reshape(-1, 1)
    kaon_ckov = np.arccos(cos_ckov_Kaon)#.reshape(-1, 1)
    proton_ckov = np.arccos(cos_ckov_Proton)#.reshape(-1,1)
