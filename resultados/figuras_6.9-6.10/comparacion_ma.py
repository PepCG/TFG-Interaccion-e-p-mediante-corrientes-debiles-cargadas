import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','grid','notebook'])

# cargamos los datos
angulos_lepton = np.load("angulos_lepton.npy")
angulos_lepton_ma1270 = np.load("angulos_lepton_ma1270.npy")
angulos_lepton_ma1350 = np.load("angulos_lepton_ma1350.npy")

angulos_nucleon = np.load("angulos_nucleon.npy")
angulos_nucleon_ma1270 = np.load("angulos_nucleon_ma1270.npy")
angulos_nucleon_ma1350 = np.load("angulos_nucleon_ma1350.npy")

omega_lepton = np.load("omega_lepton.npy")
omega_lepton_ma1270 = np.load("omega_lepton_ma1270.npy")
omega_lepton_ma1350 = np.load("omega_lepton_ma1350.npy")

omega_nucleon = np.load("omega_nucleon.npy")
omega_nucleon_1270 = np.load("omega_nucleon_1270.npy")
omega_nucleon_1350 = np.load("omega_nucleon_1350.npy")

sigma_lepton = np.load("sigma_lepton.npy")
sigma_lepton_ma1270 = np.load("sigma_lepton_ma1270.npy")
sigma_lepton_ma1350 = np.load("sigma_lepton_ma1350.npy")

sigma_nucleon = np.load("sigma_nucleon.npy")
sigma_nucleon_ma1270 = np.load("sigma_nucleon_ma1270.npy")
sigma_nucleon_ma1350 = np.load("sigma_nucleon_ma1350.npy")

sigma_energia_lepton = np.load("sigma_energia_lepton.npy")
sigma_energia_lepton_ma1270 = np.load("sigma_energia_lepton_ma1270.npy")
sigma_energia_lepton_ma1350 = np.load("sigma_energia_lepton_ma1350.npy")

sigma_energia_nucleon = np.load("sigma_energia_nucleon.npy")
sigma_energia_nucleon_ma1270 = np.load("sigma_energia_nucleon_ma1270.npy")
sigma_energia_nucleon_ma1350 = np.load("sigma_energia_nucleon_ma1350.npy")

# definimos los cocientes

cociente_sigma_1270 = sigma_lepton_ma1270/sigma_lepton
cociente_sigma_1350 = sigma_lepton_ma1350/sigma_lepton
cociente_sigma_nucleon_1270 = sigma_nucleon_ma1270/sigma_nucleon
cociente_sigma_nucleon_1350 = sigma_nucleon_ma1350/sigma_nucleon
cociente_energia_1270 = sigma_energia_lepton_ma1270/sigma_energia_lepton
cociente_energia_1350 = sigma_energia_lepton_ma1350/sigma_energia_lepton
cociente_energia_nucleon_1270 = sigma_energia_nucleon_ma1270/sigma_energia_nucleon
cociente_energia_nucleon_1350 = sigma_energia_nucleon_ma1350/sigma_energia_nucleon


# graficamos

plt.figure(figsize=(12,8))
    
plt.subplot(2, 2, 1)
plt.plot(angulos_lepton,sigma_lepton,'-',color='red',label=r"$M_A = 1,03 \; GeV$")
plt.plot(angulos_lepton_ma1270,sigma_lepton_ma1270,'--',color='orange',label=r"$M_A = 1,27 \; GeV$")
plt.plot(angulos_lepton_ma1350,sigma_lepton_ma1350,':',color='m',label=r"$M_A = 1,35 \; GeV$")
plt.yscale('log')
plt.xlabel(r"$θ_l$ (grados)")
plt.ylabel("dσ/dcosθ  (fm²)")
#plt.title("Comparacion seccion eficaz electromagnética y electrodébil")
plt.legend()
plt.grid(True)

plt.subplot(2, 2, 2)
plt.plot(omega_lepton,sigma_energia_lepton,'-',color='red',label=r"$M_A = 1,03 GeV$")
plt.plot(omega_lepton_ma1270,sigma_energia_lepton_ma1270,'--',color='orange',label=r"$M_A = 1,27 GeV$")
plt.plot(omega_lepton_ma1350,sigma_energia_lepton_ma1350,':',color='m',label=r"$M_A = 1,35 GeV$")
plt.yscale('log')
plt.xlabel("ω (MeV)")
plt.ylabel("dσ/dω  (fm²/MeV)")
#plt.title("Sección eficaz diferencial para el lepton")
plt.grid(True)

plt.subplot(2, 2, 3)
plt.plot(angulos_nucleon,sigma_nucleon,'-',color='red',label=r"$M_A = 1,03 GeV$")
plt.plot(angulos_nucleon_ma1270,sigma_nucleon_ma1270,'--',color='orange',label=r"$M_A = 1,27 GeV$")
plt.plot(angulos_nucleon_ma1350,sigma_nucleon_ma1350,':',color='m',label=r"$M_A = 1,35 GeV$")
plt.yscale('log')
plt.xlabel(r"$θ_n$ (grados)")
plt.ylabel("dσ/dcosθ  (fm²)")
#plt.title("Comparacion seccion eficaz electromagnética y electrodébil")
plt.grid(True)

plt.subplot(2, 2, 4)
plt.plot(omega_nucleon, sigma_energia_nucleon, '-', color='red', label=r"$M_A = 1,03 GeV$")
plt.plot(omega_nucleon_1270, sigma_energia_nucleon_ma1270, '--', color='orange', label=r"$M_A = 1,27 GeV$")
plt.plot(omega_nucleon_1350, sigma_energia_nucleon_ma1350, ':', color='m', label=r"$M_A = 1,35 GeV$")

plt.yscale('log')
plt.xlabel("ω (MeV)")
plt.ylabel("dσ/dω  (fm²/MeV)")
#plt.title("Sección eficaz diferencial para el lepton")
plt.grid(True)

plt.tight_layout()
plt.show()

# graficamos los cocientes

plt.figure(figsize=(12,8))

plt.subplot(2, 2, 1)
plt.plot(angulos_lepton,cociente_sigma_1270,'-',color='orange',label=r"$M_A = 1,27 GeV$")
plt.plot(angulos_lepton,cociente_sigma_1350,'--',color='m',label=r"$M_A = 1,35 GeV$")
#plt.yscale('log')
plt.xlabel("$θ_l$ (grados)")
plt.ylabel(r"$R_A$")
plt.legend(loc="best", fontsize=10)
plt.grid(True)

plt.subplot(2, 2, 2)
plt.plot(omega_lepton,cociente_energia_1270,'-',color='orange',label="Galster")
plt.plot(omega_lepton,cociente_energia_1350,'--',color='m',label="Galster <<")
#plt.yscale('log')
plt.xlabel("ω (MeV)")
#plt.ylabel("R")
#plt.title("Sección eficaz diferencial para el lepton")
plt.grid(True)


plt.subplot(2, 2, 3)
plt.plot(angulos_nucleon, cociente_sigma_nucleon_1270, '-', color='orange', label="Galster")
plt.plot(angulos_nucleon, cociente_sigma_nucleon_1350, '--', color='m', label="Galster <<")
#plt.yscale('log')
plt.xlabel("$θ_n$ (grados)")
plt.ylabel(r"$R_A$")
plt.grid(True)

plt.subplot(2, 2, 4)
plt.plot(omega_nucleon, cociente_energia_nucleon_1270, '-', color='orange', label="Galster")
plt.plot(omega_nucleon, cociente_energia_nucleon_1350, '--', color='m', label="Galster <<")
#plt.yscale('log')
plt.xlabel("ω (MeV)")
#plt.ylabel("R")
plt.grid(True)

plt.show()
