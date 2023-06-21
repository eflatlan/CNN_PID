# Get Ckov angle for every photon cluster candidate: 
# segmented into regions by mass-hypothesis

# IP : * theta_p, phi_p : MIP Azimuthal and polar angles
#      * x_MIP, y_MIP : impact points of MIP on projected on Photocathode
       
#       x, y : photon impact coords projected on Photocathode plane
#       

# Constants : 
class Constants:
    def __init__(self):
        self.defaultPhotonEnergy = 6.75
        self.CH4GapWidth = 8.0                         # t_gap
        self.RadiatorWidth = 1.0                       # r_w
        self.QuartzWindowWidth = 0.5                   # q_w
        self.EmissionLength = self.RadiatorWidth / 2.0 # L 


def FromCkov(phi_p, theta_p):

    constants = Constants()
    # a = [(r_w-L)+q_w + t_gap] * tan(theta_p)
    a = ((constants.RadiatorWidth-constants.EmissionLength) + constants.QuartzWindowWidth + constants.CH4GapWidth) * tan(theta_p)

# 4.1 : tan_phi = (y - L * tan_theta_p * sin_phi_p ) / (x - L * tan_theta_p * cos_phi_p )
#phi = np.atan((y-L *  tan_theta_p * sin_phi_p ) / (x - L * tan_theta_p * cos_phi_p ))

    # phi of photons distributed around MIP 0..2pi
    phi = np.random(0, 2*np.pi)

    # eta_c normally distributed around theta_c w 0.08 stddev
    eta_c = Gauss(theta_c, 0.08)

    # 4.4 cos etac = sin theta p cos(phi - phi_p) + cos tetha_p cos theta  
    theta = np.acos((np.cos(eta_c) -np.sin(theta_p) * cos(phi-phi_p))/np.cos(theta_p)) 

    # theta is known, get b (4.3): 

    # 4.2 : R^2 = [a cos_theta_p - b * cos_theta]^2 + [a*sin_theta_p - b*sin_theta]^2
    #           = (x_p - x)^2 + (y_p - y)^2
    # solve for b : 


    # eq 4.3 : b = (r_w - L) * tan_theta \
    # + n_f * sin_theta * (q_w /(sqrt(n_q^2-n_f^2*(sin_theta)^2)) + t_gap /(sqrt(n_g^2-n_f^2*(sin_theta)^2)))
    b = (constants.RadiatorWidth-constants.EmissionLength) * tan_theta \
        + n_f * sin_theta * (constants.QuartzWindowWidth /(np.sqrt(n_q * n_q - n_f*n_f * sin_theta * sin_theta))\
        + constants.CH4GapWidth  /(np.sqrt(n_q * n_q - n_f*n_f * sin_theta * sin_theta)))

    # a, b, phi, phi_p known: get x, y


