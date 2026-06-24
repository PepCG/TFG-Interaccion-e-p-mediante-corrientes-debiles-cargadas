import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','grid','notebook'])

# ÁNGULOS
angulos_lepton = np.load("angulos_lepton.npy")
angulos_antilepton = np.load("angulos_antilepton.npy")
angulos_neutron = np.load("angulos_neutron.npy")
angulos_proton = np.load("angulos_proton.npy")

# OMEGA
omega_lepton = np.load("omega_lepton.npy")
omega_antilepton = np.load("omega_antilepton.npy")
omega_neutron = np.load("omega_neutron.npy")
omega_proton = np.load("omega_proton.npy")

# SIGMA
sigma_lepton = np.load("sigma_lepton.npy")
sigma_antilepton = np.load("sigma_antilepton.npy")
sigma_neutron = np.load("sigma_neutron.npy")
sigma_proton = np.load("sigma_proton.npy")

# SIGMA ENERGÍA
sigma_energia_lepton = np.load("sigma_energia_lepton.npy")
sigma_energia_antilepton = np.load("sigma_energia_antilepton.npy")
sigma_energia_neutron = np.load("sigma_energia_neutron.npy")
sigma_energia_proton = np.load("sigma_energia_proton.npy")

# cargamos los datos a diferentes masas axiales

angulos_103_electron = np.load("angulos_electron_103.npy")
angulos_127_electron = np.load("angulos_electron_127.npy")
angulos_135_electron = np.load("angulos_electron_135.npy")

sigma_103_electron = np.load("sigma_electron_103.npy")
sigma_127_electron = np.load("sigma_electron_127.npy")
sigma_135_electron = np.load("sigma_electron_135.npy")

q2_103_electron = np.load("q2_electron_103.npy")
q2_127_electron = np.load("q2_electron_127.npy")
q2_135_electron = np.load("q2_electron_135.npy")

angulos_103_positron = np.load("angulos_positron_103.npy")
angulos_127_positron = np.load("angulos_positron_127.npy")
angulos_135_positron = np.load("angulos_positron_135.npy")

sigma_103_positron = np.load("sigma_positron_103.npy")
sigma_127_positron = np.load("sigma_positron_127.npy")
sigma_135_positron = np.load("sigma_positron_135.npy")

q2_103_positron = np.load("q2_positron_103.npy")
q2_127_positron = np.load("q2_positron_127.npy")
q2_135_positron = np.load("q2_positron_135.npy")

# calculamos las restas

resta_103 = sigma_103_electron - sigma_103_positron
resta_127  = sigma_127_electron - sigma_127_positron
resta_135  = sigma_135_electron - sigma_135_positron

# graficamos la resta

plt.figure(figsize=(12,8))

plt.subplot(1, 2, 1)
plt.plot(angulos_103_electron,resta_103,'-',color='red',label=r"$M_A = 1,03 GeV$")
plt.plot(angulos_127_electron,resta_127,'-',color='orange',label=r"$M_A = 1,27 GeV$")
plt.plot(angulos_135_electron,resta_135,'-',color='purple',label=r"$M_A = 1,35 GeV$")
plt.xlabel("$θ_n$ (grados)")
plt.ylabel(r"$\Delta\sigma$")
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(q2_103_electron,resta_103,'-',color='red',label=r"$M_A = 1,03 GeV$")
plt.plot(q2_127_electron,resta_127,'-',color='orange',label=r"$M_A = 1,27 GeV$")
plt.plot(q2_135_electron,resta_135,'-',color='purple',label=r"$M_A = 1,35 GeV$")
plt.xlabel(r"$Q^2$ $(MeV^2)$")
plt.ylabel(r"$\Delta\sigma$")
plt.grid(True)


plt.show()
 
# graficamos las secciones eficaces

plt.figure(figsize=(12,10))

plt.subplot(2,2,1)
plt.plot(angulos_lepton, sigma_lepton, '-', color='green', label=r'$\nu$')
plt.plot(angulos_antilepton, sigma_antilepton, '-', color='blue', label=r'$\bar{\nu}$')
plt.xlabel(r'$\theta$ (°)')
plt.ylabel(r'$d\sigma/d\cos\theta$ (fm$^2$)')
plt.yscale('log')
plt.legend()
plt.grid(True)

plt.subplot(2,2,2)
plt.plot(omega_lepton, sigma_energia_lepton, '-', color='green', label=r'$\nu$')
plt.plot(omega_antilepton, sigma_energia_antilepton, '-', color='blue', label=r'$\bar{\nu}$')
plt.xlabel(r'$\omega$ (MeV)')
plt.ylabel(r'$d\sigma/d\omega$ (fm$^2$/MeV)')
plt.yscale('log')
plt.legend()
plt.grid(True)

plt.subplot(2,2,3)
plt.plot(angulos_neutron, sigma_neutron, '-', color='green', label='n')
plt.plot(angulos_proton, sigma_proton, '-', color='blue', label='p')
plt.xlabel(r'$\theta$ (°)')
plt.ylabel(r'$d\sigma/d\cos\theta$ (fm$^2$)')
plt.yscale('log')
plt.legend()
plt.grid(True)

plt.subplot(2,2,4)
plt.plot(omega_neutron, sigma_energia_neutron, '-', color='green', label='n')
plt.plot(omega_proton, sigma_energia_proton, '-', color='blue', label='p')
plt.xlabel(r'$\omega$ (MeV)')
plt.ylabel(r'$d\sigma/d\omega$ (fm$^2$/MeV)')
plt.yscale('log')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
