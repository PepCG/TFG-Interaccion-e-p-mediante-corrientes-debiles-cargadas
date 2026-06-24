# comparacion de las secciones eficaces para la interaccion elecrtomagnetica y debil
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','grid','notebook'])

# cargamos los datos
angulos_electron = np.load("angulos_electron.npy")
angulos_proton = np.load("angulos_proton.npy")
angulos_neutron = np.load("angulos_neutron.npy")
angulos_neutrino = np.load("angulos_neutrino.npy")
omega_electron = np.load("omega_electron.npy")
omega_proton = np.load("omega_proton.npy")
omega_neutron = np.load("omega_neutron.npy")
omega_neutrino = np.load("omega_neutrino.npy")
sigma_electron = np.load("sigma_electron.npy")
sigma_proton = np.load("sigma_proton.npy")
sigma_neutron = np.load("sigma_neutron.npy")
sigma_neutrino = np.load("sigma_neutrino.npy")
sigma_energia_electron = np.load("sigma_energia_electron.npy")
sigma_energia_proton = np.load("sigma_energia_proton.npy")
sigma_energia_neutron = np.load("sigma_energia_neutron.npy")
sigma_energia_neutrino = np.load("sigma_energia_neutrino.npy")

# cargamos los datos para el propagador EM en el caso débil
angulos_electron_weak = np.load("angulos_electron_weak.npy")
angulos_proton_weak = np.load("angulos_proton_weak.npy")
omega_electron_weak = np.load("omega_electron_weak.npy")
omega_proton_weak = np.load("omega_proton_weak.npy")
sigma_electron_weak = np.load("sigma_electron_weak.npy")
sigma_proton_weak = np.load("sigma_proton_weak.npy")
sigma_energia_electron_weak = np.load("sigma_energia_electron_weak.npy")
sigma_energia_proton_weak = np.load("sigma_energia_proton_weak.npy")

# graficamos

plt.figure(figsize=(12,8))
    
plt.subplot(2, 2, 1)
plt.plot(angulos_electron,sigma_electron,'-',color='blue',label="Int. EM")
plt.plot(angulos_electron_weak,sigma_electron_weak,'-',color='green',label=r"Int. débil modificada")
plt.yscale('log')
plt.xlabel(r"$θ_l$ (grados)")
plt.ylabel("dσ/dcosθ  (fm²)")
#plt.title("Comparacion seccion eficaz electromagnética y electrodébil")
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 2)
plt.plot(omega_electron,sigma_energia_electron,'-',color='blue',label="Int. EM")
plt.plot(omega_electron_weak,sigma_energia_electron_weak,'-',color='green',label=r"Int. débil modificada")
plt.yscale('log')
plt.xlabel("ω (MeV)")
plt.ylabel("dσ/dω  (fm²)")
#plt.title("Sección eficaz diferencial para el lepton")
plt.grid(True)
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(angulos_proton,sigma_proton,'-',color='blue',label="Int. EM")
plt.plot(angulos_proton_weak,sigma_proton_weak,'-',color='green',label="Int. débil modificada")
plt.yscale('log')
plt.xlabel(r"$θ_n$ (grados)")
plt.ylabel("dσ/dcosθ  (fm²)")
#plt.title("Comparacion seccion eficaz electromagnética y electrodébil")
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 4)
plt.plot(omega_proton,sigma_energia_proton,'-',color='blue',label="Int. EM")
plt.plot(omega_proton_weak,sigma_energia_proton_weak,'-',color='green',label="Int. débil modificada")
plt.yscale('log')
plt.xlabel("ω (MeV)")
plt.ylabel("dσ/dω  (fm²/MeV)")
#plt.title("Sección eficaz diferencial para el lepton")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()