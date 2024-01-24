
# import multiprocessing
# def process_track(params):
#     recon, cluX, cluY, thetaP, phiP = params
#     # Call findPhotCkov and receive updated thetaCer and phiCer
#     updated_thetaCer, updated_phiCer = recon.findPhotCkov(cluX, cluY)



#     Alisigma2_ = cppyy.gbl.Alisigma2_(1.2904)

#     #Double_t aliSigma2(Double_t trkTheta,Double_t trkPhi, Double_t ckovTh, Double_t ckovPh)
#     # Call the sigLoc method

#     print(f"updated_thetaCer {updated_thetaCer} updated_phiCer {updated_phiCer}")

#     sigma2 = Alisigma2_.sigma2(thetaP, phiP,  updated_thetaCer, updated_phiCer)
#     print(f"sigma2 {sigma2}")

#     return updated_thetaCer, updated_phiCer, sigma2


# #theta_cer_padded, phi_cer_padded = parallel_process(recon, x_padded, y_padded, self.ThetaP[0], self.PhiP[0])

# def parallel_process(recon, x_padded, y_padded, thetaP, phiP):
#     params_list = [(recon, cluX, cluY, thetaP, phiP) for cluX, cluY in zip(x_padded, y_padded)]

#     with multiprocessing.Pool() as pool:
#         updated_values = pool.map(process_track, params_list)

#     # Unpack and update the theta_cer_padded and phi_cer_padded arrays
#     theta_cer_padded = []
#     phi_cer_padded = []
#     for updated_thetaCer, updated_phiCer in updated_values:
#         theta_cer_padded.append(updated_thetaCer)
#         phi_cer_padded.append(updated_phiCer)

#     return theta_cer_padded, phi_cer_padded

# def get_other_tracks(event_data_dict, x_padded, y_padded):

#     print(f"called get_other_tracks(event_data_dict, x_padded, y_padded)")



#     class ClassForPhotonProbability:

#         def __init__(self):
#             dict = ["Momentum", "xRad", "yRad", "xMip", "yMip", "ThetaP", "PhiP"]

#             # Initialize all attributes to None or a default value
#             self.Momentum = None
#             self.xRad = None
#             self.yRad = None
#             self.xMip = None
#             self.yMip = None
#             self.ThetaP = None
#             self.PhiP = None

#             self.fTrkPos = [0, 0]
#             self.fMipPos = [0, 0]
#             #self.fTrkDir = [self.xMip, self.yMip]



#         def populate_attributes(self):
#             # Iterate over each key and fill the class attributes with corresponding values

#             setattr(self, "Momentum", event_data_dict["Momentum"])
#             setattr(self, "xRad", event_data_dict["xRad"])
#             setattr(self, "yRad", event_data_dict["yRad"])
#             setattr(self, "xMip", event_data_dict["xMip"])
#             setattr(self, "yMip", event_data_dict["yMip"])
#             setattr(self, "ThetaP", event_data_dict["ThetaP"])
#             setattr(self, "PhiP", event_data_dict["PhiP"])




#             self.fTrkPos = [self.xRad[0], self.yRad[1]]
#             self.fMipPos = [self.xMip[0], self.yMip[1]]

#             recon = Recon(self.fTrkPos, self.fMipPos, self.ThetaP[0])

#             phiCer, thetaCer = 0, 0



#             # Example usage
#             theta_cer_padded, phi_cer_padded = parallel_process(recon, x_padded, y_padded, self.ThetaP[0], self.PhiP[0])


#             # x_padded, y_padded

#             #Double_t aliSigma2(Double_t trkTheta,Double_t trkPhi, Double_t ckovTh, Double_t ckovPh)
#             # Call the sigLoc method


#             #result = Alisigma2_.sigLoc(1.0, 0.5, 0.3, 0.4, 0.9)
#             print(f"shape theta_cer_padded {theta_cer_padded.shape}")


#     photon_prob = ClassForPhotonProbability()
#     print("photon_prob = ClassForPhotonProbability()")

#     # Populate the attributes using event_data_dict
#     photon_prob.populate_attributes()



#     #return theta_c_hyps
