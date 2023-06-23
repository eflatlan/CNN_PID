import numpy as np
import math
import matplotlib.pyplot as plt
# Constants
t_gap = 8  # example values, please replace with actual constants
phi_prime = 3.14
theta_p = 1
r_w = 1
L = 0.5
q_w = 0.5
n_f = 1.3
n_q = 1.4
n_g = 1
mass_Pion = 0.1396
mass_Kaon = 0.4937
mass_Proton = 0.938

def get_theta_p(momentum):
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

    ckovAngle = np.arccos(cos_ckov)

    return ckovAngle

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
    # Calculate cosine value using the identity cos(x) = sqrt(1 - sin(x)^2)
    cos_value = np.sqrt(1 - sin_value**2)

    # Calculate tangent value using the identity tan(x) = sin(x) / cos(x)
    tan_value = sin_value / cos_value

    return tan_value

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
momentum_values = np.arange(0.4, 5.01, 0.01)

# Cerenkov angles for each momentum value
ckov_pion = [calcCkovFromMass(p, mass_Pion) for p in momentum_values]
ckov_kaon = [calcCkovFromMass(p, mass_Kaon) for p in momentum_values]
ckov_proton = [calcCkovFromMass(p, mass_Proton) for p in momentum_values]
theta_p_values = calculate_theta_p_values(momentum_values)


# Compute |R| for each particle type and for each theta_p
R_pion = [compute_R(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap) 
          for theta_c, theta_p in zip(ckov_pion, theta_p_values)]
R_kaon = [compute_R(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap) 
          for theta_c, theta_p in zip(ckov_kaon, theta_p_values)]
R_proton = [compute_R(theta_c, phi_prime, theta_p, r_w, L, q_w, n_f, n_q, n_g, t_gap) 
            for theta_c, theta_p in zip(ckov_proton, theta_p_values)]

# Plot |R| as a function of momentum for each particle type
plt.plot(momentum_values, R_pion, label='Pion')
plt.plot(momentum_values, R_kaon, label='Kaon')
plt.plot(momentum_values, R_proton, label='Proton')

# Labels and legend
plt.xlabel('Momentum')
plt.ylabel('|R|')
plt.title('|R| as a function of Momentum')
plt.legend()
plt.grid(True)
plt.show()



#... previous code

# Plot theta_p as a function of momentum
plt.figure()  # create a new figure
plt.plot(momentum_values, theta_p_values*180/np.pi, label='theta_p')
plt.xlabel('Momentum')
plt.ylabel('theta_p')
plt.title('theta_p as a function of Momentum')
plt.legend()
plt.grid(True)
plt.show()
