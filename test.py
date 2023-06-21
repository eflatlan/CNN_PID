
# 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# IP : * theta_p, phi_p : MIP Azimuthal and polar an_gles
#      * x_MIP, y_MIP : impact points of MIP on projected on Photocathode       
#       

# Constants : 
class Constants:
    def __init__(self):
        self.defaultPhotonEnergy = 6.75
        self.CH4GapWidth = 8.0                         # t_gap
        self.RadiatorWidth = 1.0                       # r_w
        self.QuartzWindowWidth = 0.5                   # q_w
        self.EmissionLen_gth = self.RadiatorWidth / 2.0 # L 

class DataPointGenerator:
    def __init__(self, phi_p = None, theta_p = None, x_p = None, y_p = None, eta_c = None, photon_energy = None):
        self.r_w =  1
        self.L = self.r_w/2
        self.q_w =  0.5
        self.t_gap = 8
        self.theta_p = theta_p
        self.phi_p = phi_p
        self.x_p = y_p
        self.y_p = y_p
        self.eta_c = eta_c
        self.cos_tehta_p = np.cos(theta_p)
        self.sin_tehta_p = np.sin(theta_p)
        self.tan_tehta_p = np.sin(theta_p)
        self.n_f = self.GetFreonIndexOfRefraction(photon_energy)
        self.n_q = self.GetQuartzIndexOfRefraction(photon_energy)
        # n_f , n_q , n_g = freon, quartz and methane refraction indices.
        self.n_g = 1 # Giacomo's macro
        self.a = ((self.r_w-self.L) + self.q_w + self.t_gap) * self.tan_tehta_p


    # 4.4 cos etac = sin theta p cos(phi - phi_p) + cos theta_p cos theta  
    # return tehta
    def phi_2_cos(self, phi = 0):
        valid_cos_theta = False
        while not valid_cos_theta:
            cos_theta = (np.cos(self.eta_c) - self.sin_tehta_p * np.cos(phi - self.phi_p)) / self.cos_tehta_p
            print(f"phi_2_cos : cos_theta = {cos_theta} | phi = {phi}")

            if -1 <= valid_cos_theta <= 1:
                valid_cos_theta = True
            else:
                phi += 0.1  # Increase phi by a small increment for the next iteration
        # also have to decrease phi to check other ran_ge?
        return phi, np.arccos(cos_theta)
    

    # eq 4.3 : b = (r_w - L) * tan_theta \
    # + n_f * sin_theta * (q_w /(sqrt(n_q^2-n_f^2*(sin_theta)^2)) + t_gap /(sqrt(n_g^2-n_f^2*(sin_theta)^2)))
    def b_from_theta(self, theta):
        sin_theta = np.sin(theta)
        tan_theta = np.tan(theta)
        b = (self.r_w - self.L) * tan_theta \
        + self.n_f * sin_theta * (self.q_w /(np.sqrt(self.n_q**2-self.n_f**2*(sin_theta)**2)) + self.t_gap /(np.sqrt(self.n_g**2-self.n_f**2*(sin_theta)**2)))
        return b


    # 4.2 : R^2 = [a cos_theta_p - b * cos_theta]^2 + [a*sin_theta_p - b*sin_theta]^2 ; R =  |R| * [cos_phi', sin_phi']
    #           
    # solve for b : 
    def b_from_phi(self, phi, phi_local):
        tan_phi_local = np.tan(phi_local)
        cos_phi = np.cos(phi)
        b = self.a* (tan_phi_local * self.cos_tehta_p - self.sin_tehta_p) /(tan_phi_local * cos_phi - np.sqrt(1-cos_phi * cos_phi))
        return b

    def GetFreonIndexOfRefraction(self, photon_energy = None):
        n_f = 1.177 + (0.0172) * photon_energy
        return n_f

    def GetQuartzIndexOfRefraction(self, photon_energy = None):
        n_q = np.sqrt(1 + 46.411 / (113.763556 - photon_energy) + 228.71 / (328.51563 - photon_energy))
        return n_q

    def b_equals(self, phi, phi_local):

        phi_values = []
        theta_values = []
        b_phi_values = []
        b_theta_values = []

        def f(phi):
            phi, theta = self.phi_2_cos(phi)
            b_theta = self.b_from_theta(theta)
            b_phi = self.b_from_phi(phi, phi_local)
            return b_theta - b_phi

        def df(phi):
            h = 1e-4  # small number for finite difference
            return (f(phi+h) - f(phi-h)) / (2*h)

        # initial guess for phi
        phi = 0
        for _ in range(1000):  # iterate up to 1000 times
            phi, theta = self.phi_2_cos(phi)
            b_theta = self.b_from_theta(theta)
            b_phi = self.b_from_phi(phi, phi_local)
            phi_corrected = phi % (2*np.pi)

            phi_values.append(phi_corrected)
            theta_values.append(theta)
            b_phi_values.append(b_phi)
            b_theta_values.append(b_theta)

            f_phi = f(phi)
            df_phi = df(phi)

            if abs(f_phi) < 1e-6:  # if close enough to 0, then stop
                break

            phi = phi - f_phi / df_phi  # Newton's method update

        # Convert lists to numpy arrays
        phi_values = np.array(phi_values)
        theta_values = np.array(theta_values)
        b_phi_values = np.array(b_phi_values)
        b_theta_values = np.array(b_theta_values)

        # Create plots
        fig, axs = plt.subplots(2)
        axs[0].plot(phi_values, label='phi')
        axs[0].plot(theta_values, label='theta')
        axs[0].legend()

        axs[1].plot(b_phi_values, label='b_phi')
        axs[1].plot(b_theta_values, label='b_theta')
        axs[1].legend()

        plt.show()

        return theta, phi

    def plot_vector(self,theta, phi, theta_p, phi_p):
        fig = plt.figure()

        # Create 3D axes
        ax = fig.add_subplot(111, projection='3d')

        # Define the photon vector
        x_photon = np.sin(theta) * np.cos(phi)
        y_photon = np.sin(theta) * np.sin(phi)
        z_photon = np.cos(theta)

        # Define the mip vector
        x_mip = np.sin(theta_p) * np.cos(phi_p)
        y_mip = np.sin(theta_p) * np.sin(phi_p)
        z_mip = np.cos(theta_p)

        # Plot the photon vector
        ax.quiver(0, 0, 0, x_photon, y_photon, z_photon, color='b', label='Photon')

        # Plot the mip vector
        ax.quiver(0, 0, 0, x_mip, y_mip, z_mip, color='r', label='MIP')

        # Set the limits
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.set_zlim([0, 1])

        # Set the labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        # Add a legend
        ax.legend()

        plt.show()




eta_c = 0.67
photon_energy = 6.75
theta_p = 5*np.pi/180
phi_p = 5*np.pi/180
x_p, y_p = 2, 2
photon = DataPointGenerator(phi_p = phi_p, theta_p = theta_p, x_p = x_p, y_p = y_p, eta_c = eta_c, photon_energy = photon_energy)
phi = 0.1
phi_local = 0.2
theta, phi = photon.b_equals(phi, phi_local)

photon.plot_vector(theta, phi, theta_p, phi_p)