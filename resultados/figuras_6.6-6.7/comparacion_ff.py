# comparamos segun la parametrizacion

import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','grid','notebook'])

# ============================
# CARGA DE DATOS
# ============================

angulos_neutron = np.load("angulos_neutron.npy")
angulos_neutron_galster = np.load("angulos_neutron_galster.npy")
angulos_neutron_kelly = np.load("angulos_neutron_kelly.npy")
angulos_neutron_as = np.load("angulos_neutron_ka.npy")
angulos_neutron_gk = np.load("angulos_neutron_gk.npy")
angulos_neutron_gkex = np.load("angulos_neutron_gkex.npy")

angulos_proton = np.load("angulos_proton.npy")
angulos_proton_galster = np.load("angulos_proton_galster.npy")
angulos_proton_kelly = np.load("angulos_proton_kelly.npy")
angulos_proton_as = np.load("angulos_proton_ka.npy")
angulos_proton_gk = np.load("angulos_proton_gk.npy")
angulos_proton_gkex = np.load("angulos_proton_gkex.npy")

omega_neutron = np.load("omega_neutron.npy")
omega_neutron_galster = np.load("omega_neutron_galster.npy")
omega_neutron_kelly = np.load("omega_neutron_kelly.npy")
omega_neutron_as = np.load("omega_neutron_ka.npy")
omega_neutron_gk = np.load("omega_neutron_gk.npy")
omega_neutron_gkex = np.load("omega_neutron_gkex.npy")

omega_proton = np.load("omega_proton.npy")
omega_proton_galster = np.load("omega_proton_galster.npy")
omega_proton_kelly = np.load("omega_proton_kelly.npy")
omega_proton_as = np.load("omega_proton_ka.npy")
omega_proton_gk = np.load("omega_proton_gk.npy")
omega_proton_gkex = np.load("omega_proton_gkex.npy")

sigma_neutron = np.load("sigma_neutron.npy")
sigma_neutron_galster = np.load("sigma_neutron_galster.npy")
sigma_neutron_kelly = np.load("sigma_neutron_kelly.npy")
sigma_neutron_as = np.load("sigma_neutron_ka.npy")
sigma_neutron_gk = np.load("sigma_neutron_gk.npy")
sigma_neutron_gkex = np.load("sigma_neutron_gkex.npy")

sigma_proton = np.load("sigma_proton.npy")
sigma_proton_galster = np.load("sigma_proton_galster.npy")
sigma_proton_kelly = np.load("sigma_proton_kelly.npy")
sigma_proton_as = np.load("sigma_proton_ka.npy")
sigma_proton_gk = np.load("sigma_proton_gk.npy")
sigma_proton_gkex = np.load("sigma_proton_gkex.npy")

sigma_energia_neutron = np.load("sigma_energia_neutron.npy")
sigma_energia_neutron_galster = np.load("sigma_energia_neutron_galster.npy")
sigma_energia_neutron_kelly = np.load("sigma_energia_neutron_kelly.npy")
sigma_energia_neutron_as = np.load("sigma_energia_neutron_as.npy")
sigma_energia_neutron_gk = np.load("sigma_energia_neutron_gk.npy")
sigma_energia_neutron_gkex = np.load("sigma_energia_neutron_gkex.npy")

sigma_energia_proton = np.load("sigma_energia_proton.npy")
sigma_energia_proton_galster = np.load("sigma_energia_proton_galster.npy")
sigma_energia_proton_kelly = np.load("sigma_energia_proton_kelly.npy")
sigma_energia_proton_as = np.load("sigma_energia_proton_as.npy")
sigma_energia_proton_gk = np.load("sigma_energia_proton_gk.npy")
sigma_energia_proton_gkex = np.load("sigma_energia_proton_gkex.npy")



# ============================
# COCIENTES
# ============================

cociente = sigma_neutron / sigma_neutron
cociente_galster = sigma_neutron / sigma_neutron_galster
cociente_kelly = sigma_neutron / sigma_neutron_kelly
cociente_as = sigma_neutron / sigma_neutron_as
cociente_gk = sigma_neutron / sigma_neutron_gk
cociente_gkex = sigma_neutron / sigma_neutron_gkex

cociente_energia = sigma_energia_neutron / sigma_energia_neutron
cociente_energia_galster = sigma_energia_neutron / sigma_energia_neutron_galster
cociente_energia_kelly = sigma_energia_neutron / sigma_energia_neutron_kelly
cociente_energia_as = sigma_energia_neutron / sigma_energia_neutron_as
cociente_energia_gk = sigma_energia_neutron / sigma_energia_neutron_gk
cociente_energia_gkex = sigma_energia_neutron / sigma_energia_neutron_gkex

cociente_proton = sigma_proton / sigma_proton
cociente_proton_galster = sigma_proton / sigma_proton_galster
cociente_proton_kelly = sigma_proton / sigma_proton_kelly
cociente_proton_as = sigma_proton / sigma_proton_as
cociente_proton_gk = sigma_proton / sigma_proton_gk
cociente_proton_gkex = sigma_proton / sigma_proton_gkex

cociente_proton_energia = sigma_energia_proton / sigma_energia_proton
cociente_proton_energia_galster = sigma_energia_proton / sigma_energia_proton_galster
cociente_proton_energia_kelly = sigma_energia_proton / sigma_energia_proton_kelly
cociente_proton_energia_as = sigma_energia_proton / sigma_energia_proton_as
cociente_proton_energia_gk = sigma_energia_proton / sigma_energia_proton_gk
cociente_proton_energia_gkex = sigma_energia_proton / sigma_energia_proton_gkex



# ============================
# REPRESENTACIÓN DE COCIENTES
# ============================

plt.figure(figsize=(12,8))

plt.subplot(2, 2, 1)
plt.plot(angulos_neutron, cociente, '-', color='blue', label="Galster")
plt.plot(angulos_neutron_kelly, cociente_kelly, '-', color='red', label="Kelly")
plt.plot(angulos_neutron_as, cociente_as, '-', color='orange', label="Arrington-Sick")
plt.plot(angulos_neutron_gk, cociente_gk, '--', color='purple', label="GK")
plt.plot(angulos_neutron_gkex, cociente_gkex, '--', color='green', label="GKex")
plt.xlabel("θ (grados)")
plt.ylabel("R")
plt.grid(True)

plt.subplot(2, 2, 2)
plt.plot(omega_neutron, cociente_energia, '-', color='blue', label="Galster")
plt.plot(omega_neutron_kelly, cociente_energia_kelly, '-', color='red', label="Kelly")
plt.plot(omega_neutron_as, cociente_energia_as, '-', color='orange', label="Arrington-Sick")
plt.plot(omega_neutron_gk, cociente_energia_gk, '--', color='purple', label="GK")
plt.plot(omega_neutron_gkex, cociente_energia_gkex, '--', color='green', label="GKex")
plt.xlabel(r"$Q^2$ $(GeV^2)$")
plt.legend(loc="best", fontsize=10)
plt.grid(True)


plt.subplot(2, 2, 3)
plt.plot(angulos_proton, cociente_proton, '-', color='blue', label="Galster")
plt.plot(angulos_proton_kelly, cociente_proton_kelly, '-', color='red', label="Kelly")
plt.plot(angulos_proton_as, cociente_proton_as, '-', color='orange', label="Arrington-Sick")
plt.plot(angulos_proton_gk, cociente_proton_gk, '--', color='purple', label="GK")
plt.plot(angulos_proton_gkex, cociente_proton_gkex, '--', color='green', label="GKex")
plt.xlabel("θ (grados)")
plt.ylabel("R")
plt.grid(True)

plt.subplot(2, 2, 4)
plt.plot(omega_proton, cociente_proton_energia, '-', color='blue', label="Galster")
plt.plot(omega_proton_kelly, cociente_proton_energia_kelly, '-', color='red', label="Kelly")
plt.plot(omega_proton_as, cociente_proton_energia_as, '-', color='orange', label="Arrington-Sick")
plt.plot(omega_proton_gk, cociente_proton_energia_gk, '--', color='purple', label="GK")
plt.plot(omega_proton_gkex, cociente_proton_energia_gkex, '--', color='green', label="GKex")
plt.xlabel("$Q^2$ $(GeV^2)$")
plt.grid(True)



# ============================
# REPRESENTACIÓN DE SIGMAS
# ============================

plt.figure(figsize=(12,8))
    
plt.subplot(2, 2, 1)
plt.plot(angulos_neutron, sigma_neutron, '-', color='blue', label="Galster")
plt.plot(angulos_neutron_kelly, sigma_neutron_kelly, '-', color='red', label="Kelly")
plt.plot(angulos_neutron_as, sigma_neutron_as, '-', color='orange', label="Arrington-Sick")
plt.plot(angulos_neutron_gk, sigma_neutron_gk, '--', color='purple', label="GK")
plt.plot(angulos_neutron_gkex, sigma_neutron_gkex, '--', color='green', label="GKex")
plt.yscale('log')
plt.xlabel("θ (grados)")
plt.ylabel("dσ/dcosθ  (fm²)")
plt.grid(True)

plt.subplot(2, 2, 2)
plt.plot(omega_neutron, sigma_neutron, '-', color='blue', label="Galster")
plt.plot(omega_neutron_kelly, sigma_neutron_kelly, '-', color='red', label="Kelly")
plt.plot(omega_neutron_as, sigma_neutron_as, '-', color='orange', label="Arrington-Sick")
plt.plot(omega_neutron_gk, sigma_neutron_gk, '--', color='purple', label="GK")
plt.plot(omega_neutron_gkex, sigma_neutron_gkex, '--', color='green', label="GKex")
plt.yscale('log')
plt.xlabel("$Q^2$ $(GeV^2)$")
plt.legend(loc="best")
plt.grid(True)


plt.subplot(2, 2, 3)
plt.plot(angulos_proton, sigma_proton, '-', color='blue', label="Galster")
plt.plot(angulos_proton_kelly, sigma_proton_kelly, '-', color='red', label="Kelly")
plt.plot(angulos_proton_as, sigma_proton_as, '-', color='orange', label="Arrington-Sick")
plt.plot(angulos_proton_gk, sigma_proton_gk, '--', color='purple', label="GK")
plt.plot(angulos_proton_gkex, sigma_proton_gkex, '--', color='green', label="GKex")
plt.yscale('log')
plt.xlabel("θ (grados)")
plt.ylabel("dσ/dcosθ  (fm²)")
plt.grid(True)

plt.subplot(2, 2, 4)
plt.plot(omega_proton, sigma_proton, '-', color='blue', label="Galster")
plt.plot(omega_proton_kelly, sigma_proton_kelly, '-', color='red', label="Kelly")
plt.plot(omega_proton_as, sigma_proton_as, '-', color='orange', label="Arrington-Sick")
plt.plot(omega_proton_gk, sigma_proton_gk, '--', color='purple', label="GK")
plt.plot(omega_proton_gkex, sigma_proton_gkex, '--', color='green', label="GKex")
plt.yscale('log')
plt.xlabel("$Q^2$ $(GeV^2)$")
plt.grid(True)

plt.tight_layout()
plt.show()





