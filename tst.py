pH, pL = 0.5, 0.4
deg_theta_H, deg_theta_L = 42.5, 50


momentum = 0.2


deg_theta_P = deg_theta_L + (deg_theta_H - deg_theta_L) / (pH - pL) * (momentum - pL)
deg_theta_P

print(deg_theta_P)

