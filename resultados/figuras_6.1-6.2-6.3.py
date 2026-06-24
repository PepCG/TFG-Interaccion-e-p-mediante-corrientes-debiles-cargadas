import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','grid','notebook'])

# int EM
angulos_lepton_500 = np.load("angulos_lepton_500.npy")
angulos_lepton_1000 = np.load("angulos_lepton_1000.npy")
angulos_lepton_2000 = np.load("angulos_lepton_2000.npy")
angulos_nucleon_500 = np.load("angulos_nucleon_500.npy")
angulos_nucleon_1000 = np.load("angulos_nucleon_1000.npy")
angulos_nucleon_2000 = np.load("angulos_nucleon_2000.npy")
omega_lepton_500 = np.load("omega_lepton_500.npy")
omega_lepton_1000 = np.load("omega_lepton_1000.npy")
omega_lepton_2000 = np.load("omega_lepton_2000.npy")
omega_nucleon_500 = np.load("omega_nucleon_500.npy")
omega_nucleon_1000 = np.load("omega_nucleon_1000.npy")
omega_nucleon_2000 = np.load("omega_nucleon_2000.npy")
sigma_lepton_500 = np.load("sigma_lepton_500.npy")
sigma_lepton_1000 = np.load("sigma_lepton_1000.npy")
sigma_lepton_2000 = np.load("sigma_lepton_2000.npy")
sigma_nucleon_500 = np.load("sigma_nucleon_500.npy")
sigma_nucleon_1000 = np.load("sigma_nucleon_1000.npy")
sigma_nucleon_2000 = np.load("sigma_nucleon_2000.npy")
sigma_energia_lepton_500 = np.load("sigma_energia_lepton_500.npy")
sigma_energia_lepton_1000 = np.load("sigma_energia_lepton_1000.npy")
sigma_energia_lepton_2000 = np.load("sigma_energia_lepton_2000.npy")
sigma_energia_nucleon_500 = np.load("sigma_energia_nucleon_500.npy")
sigma_energia_nucleon_1000 = np.load("sigma_energia_nucleon_1000.npy")
sigma_energia_nucleon_2000 = np.load("sigma_energia_nucleon_2000.npy")

# int débil

angulos_lepton_500_weak = np.load("angulos_lepton_500_weak.npy")
angulos_lepton_1000_weak = np.load("angulos_lepton_1000_weak.npy")
angulos_lepton_2000_weak = np.load("angulos_lepton_2000_weak.npy")

angulos_nucleon_500_weak = np.load("angulos_nucleon_500_weak.npy")
angulos_nucleon_1000_weak = np.load("angulos_nucleon_1000_weak.npy")
angulos_nucleon_2000_weak = np.load("angulos_nucleon_2000_weak.npy")

omega_lepton_500_weak = np.load("omega_lepton_500_weak.npy")
omega_lepton_1000_weak = np.load("omega_lepton_1000_weak.npy")
omega_lepton_2000_weak = np.load("omega_lepton_2000_weak.npy")

omega_nucleon_500_weak = np.load("omega_nucleon_500_weak.npy")
omega_nucleon_1000_weak = np.load("omega_nucleon_1000_weak.npy")
omega_nucleon_2000_weak = np.load("omega_nucleon_2000_weak.npy")

sigma_lepton_500_weak = np.load("sigma_lepton_500_weak.npy")
sigma_lepton_1000_weak = np.load("sigma_lepton_1000_weak.npy")
sigma_lepton_2000_weak = np.load("sigma_lepton_2000_weak.npy")

sigma_nucleon_500_weak = np.load("sigma_nucleon_500_weak.npy")
sigma_nucleon_1000_weak = np.load("sigma_nucleon_1000_weak.npy")
sigma_nucleon_2000_weak = np.load("sigma_nucleon_2000_weak.npy")

sigma_energia_lepton_500_weak = np.load("sigma_energia_lepton_500_weak.npy")
sigma_energia_lepton_1000_weak = np.load("sigma_energia_lepton_1000_weak.npy")
sigma_energia_lepton_2000_weak = np.load("sigma_energia_lepton_2000_weak.npy")

sigma_energia_nucleon_500_weak = np.load("sigma_energia_nucleon_500_weak.npy")
sigma_energia_nucleon_1000_weak = np.load("sigma_energia_nucleon_1000_weak.npy")
sigma_energia_nucleon_2000_weak = np.load("sigma_energia_nucleon_2000_weak.npy")



# graficamos

plt.figure(figsize=(12,8))
    
plt.subplot(1, 2, 1)
plt.plot(angulos_nucleon_500,sigma_nucleon_500,'-',color='blue',label=r"$\epsilon_i = 500 MeV$")
plt.plot(angulos_nucleon_1000,sigma_nucleon_1000,'-',color='green',label=r"$\epsilon_i = 1000 MeV$")
plt.plot(angulos_nucleon_2000,sigma_nucleon_2000,'-',color='orange',label=r"$\epsilon_i = 2000 MeV$")
plt.yscale('log')
plt.xlabel(r"$θ$ (grados)")
plt.ylabel("dσ/dcosθ  (fm²)")
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(omega_nucleon_500,sigma_energia_nucleon_500,'-',color='blue',label=r"$\epsilon_i = 1000 MeV$")
plt.plot(omega_nucleon_1000,sigma_energia_nucleon_1000,'-',color='green',label=r"$\epsilon_i = 1000 MeV$")
plt.plot(omega_nucleon_2000,sigma_energia_nucleon_2000,'-',color='orange',label=r"$\epsilon_i = 1000 MeV$")
plt.yscale('log')
plt.xlabel("ω (MeV)")
plt.ylabel("dσ/dω  (fm²/MeV)")
plt.grid(True)

plt.tight_layout()
plt.show()