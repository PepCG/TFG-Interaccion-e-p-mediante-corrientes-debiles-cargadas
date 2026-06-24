import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','grid','notebook'])

# normal
angulos_nucleon = np.load("angulos_nucleon.npy")
sigma_nucleon = np.load("sigma_nucleon.npy")
q2_nucleon = np.load("q2_nucleon.npy")


# GP = 0
angulos_nucleon_nulo = np.load("angulos_nucleon_nulo.npy")
sigma_nucleon_nulo = np.load("sigma_nucleon_nulo.npy")
q2_nucleon_nulo = np.load("q2_nucleon_nulo.npy")

# cocientes
R = sigma_nucleon_nulo / sigma_nucleon

# graficamos

plt.figure(figsize=(10,4))
    
plt.subplot(1, 2, 1)
plt.plot(angulos_nucleon,R,color='blue')
#plt.yscale('log')
plt.xlabel(r"$θ$ (grados)")
plt.ylabel(r"$R_P$")
#plt.title("Comparacion seccion eficaz electromagnética y electrodébil")
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(q2_nucleon,R,'-',color='blue')
#plt.yscale('log')
plt.xlabel(r"$Q^2$ $(GeV^2)$")
#plt.title("Sección eficaz diferencial para el lepton")
plt.grid(True)

plt.tight_layout()
plt.show()