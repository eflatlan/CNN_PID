


def plot_ckov_tracks(lower_momentum = 0.5, upper_momentum = 5):
    p = np.linspace(.5, 5, 1000)[:, None]  # Reshape p to be a column vector for broadcasting
    n = 1.3
    MASS_ELECTRON = 0.000511
    MASS_MUON = 0.1057
    MASS_PION = 0.1396
    MASS_KAON = 0.4937
    MASS_PROTON = 0.9383
    masses = np.array([MASS_ELECTRON, MASS_MUON, MASS_PION, MASS_KAON, MASS_PROTON])

    p_lim = masses/np.sqrt(n**2 - 1)
    p_squared = p ** 2
    # The broadcasting here is correct, and will result in a (1000, 5) shape array
    p_squared_plus_m_squared = p_squared + masses ** 2

    p_mask = p > p_lim  # Broadcasting comparison
    cos_theta_c_hyps = np.sqrt((p_squared_plus_m_squared)) / (p * n)
    #theta_c_hyps = np.arccos(cos_theta_c_hyps)

    print(f"{cos_theta_c_hyps.shape}")
    theta_c_hyps = np.zeros_like(cos_theta_c_hyps)
    theta_c_hyps[p_mask] = np.arccos(cos_theta_c_hyps[p_mask])
    print(f"{theta_c_hyps.shape}")

    import matplotlib.pyplot as plt
    # Plotting
    species_names = ['Electron', 'Muon', 'Pion', 'Kaon', 'Proton']
    for i in range(theta_c_hyps.shape[1]):  # Loop over the number of species
        plt.plot(p.flatten(), theta_c_hyps[:, i], label=species_names[i])

    plt.xlabel('Momentum (GeV/c)')
    plt.ylabel('Cherenkov Angle (rad)')
    plt.title('Cherenkov Angle vs. Momentum for Various Particles')
    plt.legend()
    plt.grid(True)
    plt.show()
