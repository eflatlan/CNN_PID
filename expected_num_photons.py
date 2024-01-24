import numpy as np


def expected_num_photons(n, p, pdg):

    """
        filtrer ut tracks som ikke har nok photoner i relevant Ckov-sone fra pdg

        Args:
            n : refractive index
            p : momentum
            pdg : pdg-code

        Returns :
            Limit for evaluation in whole circle and in short side of circle
            limit = landa - sqrt(landa)
            limit_half = landa/2 - sqrt(landa/2)

    """
    if abs(pdg) == 211:
        MASS = 0.1396
    if abs(pdg) == 321:
        MASS = 0.4937
    if abs(pdg) == 2212:
        MASS = 0.9383


    p_lim = MASS/np.sqrt(n**2 - 1)
    p_squared = p ** 2
    # The broadcasting here is correct, and will result in a (1000, 5) shape array
    p_squared_plus_m_squared = p_squared + MASS ** 2

    p_mask = p > p_lim  # Broadcasting comparison
    cos_theta_c = np.sqrt((p_squared_plus_m_squared)) / (p * n)
    sin_theta_c = np.sin(np.arccos(cos_theta_c))


    # 15 photons @ sin^2 = 0.4
    NUM_EXPECTED_PHOTONS_SATURATION = 15
    num_expected_photons = sin_theta_c ** 2 * NUM_EXPECTED_PHOTONS_SATURATION / 0.4

    limit_evaluate = np.floor(num_expected_photons - np.sqrt(num_expected_photons)) # poisson distribution mean - std
    limit_evaluate_half = np.floor(num_expected_photons/2 - np.sqrt(num_expected_photons/2)) # poisson distribution mean - std

    return limit_evaluate, limit_evaluate_half
