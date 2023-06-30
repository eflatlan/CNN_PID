import numpy as np
import math
import matplotlib.pyplot as plt
# Constants
t_gap = 8  # example values, please replace with actual constants
r_w = 1
L = 0.5
q_w = 0.5
n_f = 1.3
n_q = 1.4
n_g = 1
mass_Pion = 0.1396
mass_Kaon = 0.4937
mass_Proton = 0.938


# Define the function to get the index of refraction
def get_freon_index_of_refraction(photon_energy):
    x = photon_energy
    k = 1.177 + (0.0172)*x
    return k


def calc_ckov_from_mass(p_values, n, m_values):
    def ckov(p, m):
        p_sq = p**2
        cos_ckov_denom = p*n

        # sanity check
        if p_sq + m*m < 0:
            return 0

        cos_ckov = np.sqrt(p_sq + m*m) / cos_ckov_denom

        # sanity check
        if cos_ckov > 1 or cos_ckov < -1:
            return 0

        ckov_an_gle = np.arccos(cos_ckov)
        #return ckov_an_gle*180/np.pi  # convert to degrees
        return ckov_an_gle#*180/np.pi  # convert to degrees

    # Create labels for plot
    labels = ['Pion', 'Kaon', 'Proton']

    # Create a new plot
    plt.figure(figsize=(10, 5))
    
    # Calculate the Cherenkov an_gles for each species and plot
    for m, label in zip(m_values, labels):
        ckov_an_gles = [ckov(p, m) for p in p_values]
        plt.plot(p_values, ckov_an_gles, label=label)

    # Set labels and legend
    plt.xlabel('Momentum (GeV/c)')
    plt.ylabel('Cherenkov an_gle (degrees)')
    plt.title('Cherenkov an_gle as a function of momentum')
    plt.legend()
    plt.grid(True)
    plt.show()

def get_theta_p(momentum):
    if momentum < 0.4:
        deg_theta_P = 50;
    else:
        if momentum < 0.5:
            pH, pL = 0.5, 0.4
            deg_theta_H, deg_theta_L = 42.5, 50
        elif momentum < 0.7:
            pH, pL = 0.7, 0.5
            deg_theta_H, deg_theta_L = 27.5, 42.5, 
        elif momentum < 1:
            pH, pL = 1, 0.7
            deg_theta_H, deg_theta_L = 22.5, 27.5
        elif momentum < 1.5:
            pH, pL = 1.5, 1
            deg_theta_H, deg_theta_L = 15, 22.5
        elif momentum < 2.5:
            pH, pL = 2.5, 1.5
            deg_theta_H, deg_theta_L = 10, 15
        else:
            pH, pL = 5, 2.5
            deg_theta_H, deg_theta_L = 8, 10
        deg_theta_P = deg_theta_L + (deg_theta_H - deg_theta_L) / (pH - pL) * (momentum - pL)

    return deg_theta_P * np.pi / 180

# The previously defined functions...

def calculate_theta_p_values(momentum_values):
    return np.array([get_theta_p(p) for p in momentum_values])


def calcCkovFromMass(p, m):
    p_sq = p**2
    cos_ckov_denom = p*n_f

    # sanity check
    if p_sq + m*m < 0:
        return 0

    cos_ckov = np.sqrt(p_sq + m*m)/cos_ckov_denom

    # sanity check
    if cos_ckov > 1 or cos_ckov < -1:
        return 0

    ckovAn_gle = np.arccos(cos_ckov)

    return ckovAn_gle

# Function to compute t_z
def compute_t_z(t_gap, theta_p, phi_prime, theta_c, r_w, L, q_w, n_f, n_q, n_g):
    sin_theta_c = np.sin(theta_c)
    cos_phi_prime = np.cos(phi_prime)
    tan_theta_p = np.tan(theta_p)

    term1 = (r_w - L) / (np.sqrt(1 - sin_theta_c**2))
    term2 = (q_w * n_f) / (np.sqrt(n_q**2 - (n_f * sin_theta_c)**2))

    numerator = t_gap + tan_theta_p * cos_phi_prime * sin_theta_c * (term1 + term2)
    denominator = 1 - ((tan_theta_p * cos_phi_prime * sin_theta_c * n_f) / (np.sqrt(n_g**2 - (n_f * sin_theta_c)**2)))

    return numerator / denominator

def calculate_tan(sin_value):
    # Calculate cosine value usin_g the identity cos(x) = sqrt(1 - sin(x)^2)
    cos_value = np.sqrt(1 - sin_value**2)

    # Calculate tan_gent value usin_g the identity tan(x) = sin(x) / cos(x)
    tan_value = sin_value / cos_value

    return tan_value

def compute_R2(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap):

    tanThetaP = np.tan(theta_p)

    cosPhiL = np.cos(phi_prime)
    sinPhiL = np.sin(phi_prime)  # Assumed to be sin

    cosEtaC = np.cos(theta_c)
    sinEtaC = np.sin(theta_c)

    cosThetaP = np.cos(theta_p)
    tanThetaP = np.tan(theta_p)


    r_wlDeltaR = (r_w - L)/(np.sqrt(1-sinEtaC*sinEtaC))
    q_wDeltaR = (q_w * n_f)/(np.sqrt(n_q*n_q - sinEtaC*sinEtaC*n_f*n_f))
    num = (t_gap + tanThetaP * cosPhiL * sinEtaC * (r_wlDeltaR + q_wDeltaR))

    denum = 1 - (tanThetaP * cosPhiL * sinEtaC * n_f)/(np.sqrt(n_g*n_g - sinEtaC*sinEtaC*n_f*n_f))

    tZ = num/denum

    Lz = (r_w - L) + q_w + tZ

    t_gapDeltaR = (tZ * n_f)/(np.sqrt(n_g*n_g - sinEtaC*sinEtaC*n_f*n_f))

    R = sinEtaC*(r_wlDeltaR+q_wDeltaR+t_gapDeltaR) / cosThetaP

    return R


# Function to compute |R|
def termf1(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap):
    t_z = compute_t_z(t_gap, theta_p, phi_prime, theta_c, r_w, L, q_w, n_f, n_q, n_g)
    sin_theta_c = np.sin(theta_c)

    # stat test
    sin_theta_c = np.sin(theta_c)
    sin_theta_qz = n_f/n_q*sin_theta_c
    sin_theta_0 = n_f/n_g*sin_theta_c
    term1 = (r_w - L) * calculate_tan(sin_theta_c)
    term2 = q_w  * calculate_tan(sin_theta_qz)
    term3 = t_z * calculate_tan(sin_theta_0)
    #R = sin_theta_c * (term1 + term2 + term3)
    # end test

    term1 = (r_w - L) / (np.sqrt(1 - sin_theta_c**2))*sin_theta_c
    term2 = (q_w * n_f) / (np.sqrt(n_q**2 - (n_f * sin_theta_c)**2))*sin_theta_c
    term3 = (t_z * n_f) / (np.sqrt(n_g**2 - (n_f * sin_theta_c)**2))*sin_theta_c
    return term1

def termf2(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap):
    t_z = compute_t_z(t_gap, theta_p, phi_prime, theta_c, r_w, L, q_w, n_f, n_q, n_g)
    sin_theta_c = np.sin(theta_c)

    # stat test
    sin_theta_c = np.sin(theta_c)
    sin_theta_qz = n_f/n_q*sin_theta_c
    sin_theta_0 = n_f/n_g*sin_theta_c
    term1 = (r_w - L) * calculate_tan(sin_theta_c)
    term2 = q_w  * calculate_tan(sin_theta_qz)
    term3 = t_z * calculate_tan(sin_theta_0)
    #R = sin_theta_c * (term1 + term2 + term3)
    # end test

    term1 = (r_w - L) / (np.sqrt(1 - sin_theta_c**2))*sin_theta_c
    term2 = (q_w * n_f) / (np.sqrt(n_q**2 - (n_f * sin_theta_c)**2))*sin_theta_c
    term3 = (t_z * n_f) / (np.sqrt(n_g**2 - (n_f * sin_theta_c)**2))*sin_theta_c
    return term2



def termf3(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap):
    t_z = compute_t_z(t_gap, theta_p, phi_prime, theta_c, r_w, L, q_w, n_f, n_q, n_g)
    sin_theta_c = np.sin(theta_c)

    # stat test
    sin_theta_c = np.sin(theta_c)
    sin_theta_qz = n_f/n_q*sin_theta_c
    sin_theta_0 = n_f/n_g*sin_theta_c
    term1 = (r_w - L) * calculate_tan(sin_theta_c)
    term2 = q_w  * calculate_tan(sin_theta_qz)
    term3 = t_z * calculate_tan(sin_theta_0)
    #R = sin_theta_c * (term1 + term2 + term3)
    # end test

    term1 = (r_w - L) / (np.sqrt(1 - sin_theta_c**2))*sin_theta_c
    term2 = (q_w * n_f) / (np.sqrt(n_q**2 - (n_f * sin_theta_c)**2))*sin_theta_c
    term3 = (t_z * n_f) / (np.sqrt(n_g**2 - (n_f * sin_theta_c)**2))*sin_theta_c
    return term3

# Function to compute |R|
def compute_R(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap):
    t_z = compute_t_z(t_gap, theta_p, phi_prime, theta_c, r_w, L, q_w, n_f, n_q, n_g)
    sin_theta_c = np.sin(theta_c)

    # stat test
    sin_theta_c = np.sin(theta_c)
    sin_theta_qz = n_f/n_q*sin_theta_c
    sin_theta_0 = n_f/n_g*sin_theta_c
    term1 = (r_w - L) * calculate_tan(sin_theta_c)
    term2 = q_w  * calculate_tan(sin_theta_qz)
    term3 = t_z * calculate_tan(sin_theta_0)
    #R = sin_theta_c * (term1 + term2 + term3)
    # end test

    term1 = (r_w - L) / (np.sqrt(1 - sin_theta_c**2))*sin_theta_c
    term2 = (q_w * n_f) / (np.sqrt(n_q**2 - (n_f * sin_theta_c)**2))*sin_theta_c
    term3 = (t_z * n_f) / (np.sqrt(n_g**2 - (n_f * sin_theta_c)**2))*sin_theta_c
    R = (term1 + term2 + term3)
    return R / np.cos(theta_p)

# Momentum values
momentum_values = np.arange(.51, 5.01, 0.01)

# Cerenkov an_gles for each momentum value
ckov_pion = [calcCkovFromMass(p, mass_Pion) for p in momentum_values]
ckov_kaon = [calcCkovFromMass(p, mass_Kaon) for p in momentum_values]
ckov_proton = [calcCkovFromMass(p, mass_Proton) for p in momentum_values]
theta_p_values = calculate_theta_p_values(momentum_values)



# Plot theta_p as a function of momentum
plt.figure()  # create a new figure
plt.plot(momentum_values, theta_p_values*180/np.pi, label='theta_p')
plt.xlabel('Momentum')
plt.ylabel('theta_p')
plt.title('theta_p as a function of Momentum')
plt.legend()
plt.grid(True)
plt.show()


def plot_values_and_areas(momentum_values):  
    phi_prime_values = [0.3, np.pi, np.pi/2]


    fig, axs = plt.subplots(len(phi_prime_values), 1, figsize=(10, 15))  # creates a 3x1 grid of subplots for |R|
    
    # Compute |R| for each particle type and for each phi_prime value
    R_pion_phi = [[compute_R(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap) 
                   for theta_c, theta_p in zip(ckov_pion, theta_p_values)] for phi_prime in phi_prime_values]
    
    term1 = [[termf1(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap) 
                   for theta_c, theta_p in zip(ckov_pion, theta_p_values)] for phi_prime in phi_prime_values]
    
    term2 = [[termf2(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap) 
                for theta_c, theta_p in zip(ckov_pion, theta_p_values)] for phi_prime in phi_prime_values]
    
    term3 = [[termf3(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap) 
                for theta_c, theta_p in zip(ckov_pion, theta_p_values)] for phi_prime in phi_prime_values]
    
    #fig, axs = plt.subplots(len(phi_prime_values), 1, figsize=(10, 15))  # creates a 3x1 grid of subplots for |R|



    plt.figure(figsize=(10, 5))
    plt.plot(momentum_values, term1, label='Term 1')
    plt.plot(momentum_values, term2, label='Term 2')
    plt.plot(momentum_values, term3, label='Term 3')
    plt.xlabel('Momentum')
    plt.ylabel('Value')
    plt.title('Individual terms as functions of Momentum')
    plt.legend()
    plt.grid(True)
    plt.show()


    R_kaon_phi = [[compute_R(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap) 
                   for theta_c, theta_p in zip(ckov_kaon, theta_p_values)] for phi_prime in phi_prime_values]
    R_proton_phi = [[compute_R(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap) 
                     for theta_c, theta_p in zip(ckov_proton, theta_p_values)] for phi_prime in phi_prime_values]

    # Calculate area for each momentum value
    area_pion = [2 * R_pion_phi[2][i] * (R_pion_phi[0][i] + R_pion_phi[1][i]) for i in range(len(momentum_values))]
    area_kaon = [2 * R_kaon_phi[2][i] * (R_kaon_phi[0][i] + R_kaon_phi[1][i]) for i in range(len(momentum_values))]
    area_proton = [2 * R_proton_phi[2][i] * (R_proton_phi[0][i] + R_proton_phi[1][i]) for i in range(len(momentum_values))]
        
    # Plot |R| as a function of momentum for each particle type in separate subplots
    for i, phi_prime in enumerate(phi_prime_values):  
        axs[i].plot(momentum_values, R_pion_phi[i], label='Pion')
        axs[i].plot(momentum_values, R_kaon_phi[i], label='Kaon')
        axs[i].plot(momentum_values, R_proton_phi[i], label='Proton')

        # Labels and legend for each subplot
        axs[i].set_xlabel('Momentum')
        axs[i].set_ylabel('|R|')
        axs[i].set_title(f'|R| as a function of Momentum; phi_prime = {phi_prime}')
        axs[i].legend()
        axs[i].grid(True)

    plt.tight_layout()
    plt.show()

    # Create a new figure for the area plot
    plt.figure(figsize=(10, 5))
    # Plot areas as a function of momentum
    plt.plot(momentum_values, area_pion, label='Pion')
    plt.plot(momentum_values, area_kaon, label='Kaon')
    plt.plot(momentum_values, area_proton, label='Proton')
    plt.xlabel('Momentum')
    plt.ylabel('Area')
    plt.title('Area as a function of Momentum for each Particle')
    plt.legend()
    plt.grid(True)
    plt.show()

# Call the function
plot_values_and_areas(momentum_values)


photon_energy = 5.5

mass_values = [0.139570, 0.493677, 0.938272]  # masses for pion, kaon, and proton in GeV/c^2

calc_ckov_from_mass(momentum_values, get_freon_index_of_refraction(photon_energy=photon_energy), m_values=mass_values)