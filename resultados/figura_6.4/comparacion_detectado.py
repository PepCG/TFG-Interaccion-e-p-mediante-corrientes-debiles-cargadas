# comparacion de las secciones eficaces para la interaccion elecrtomagnetica y debil
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','grid','notebook'])

# cargamos los datos
angulos_electron = np.load("angulos_electron.npy")
angulos_proton = np.load("angulos_proton.npy")
omega_electron = np.load("omega_electron.npy")
omega_proton = np.load("omega_proton.npy")
sigma_electron = np.load("sigma_electron.npy")
sigma_proton = np.load("sigma_proton.npy")
sigma_energia_electron = np.load("sigma_energia_electron.npy")
sigma_energia_proton = np.load("sigma_energia_proton.npy")

# graficamos

plt.figure(figsize=(12,8))
    
# PANEL IZQUIERDO
ax1 = plt.subplot(1, 2, 1)
ax1.plot(angulos_electron, sigma_electron, '-', color='blue', label="$e^-$")
ax1.plot(angulos_proton, sigma_proton, '-', color='green', label="$p$")
ax1.set_yscale('log')
ax1.set_xlabel(r"$θ_e$ (grados)")
ax1_top = ax1.secondary_xaxis('top')
ax1_top.set_xlabel(r"$θ_p$ (grados)")
ax1.set_ylabel("dσ/dcosθ  (fm²)")
ax1.legend()
ax1.grid(True)

# PANEL DERECHO
ax2 = plt.subplot(1, 2, 2)
ax2.plot(omega_electron, sigma_energia_electron, '-', color='blue', label="$e^-$")
ax2.plot(omega_proton, sigma_energia_proton, '-', color='green', label="$p$")
ax2.set_yscale('log')
ax2.set_xlabel("ω (MeV)")
ax2.set_ylabel("dσ/dω (fm²/MeV)")
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.show()
