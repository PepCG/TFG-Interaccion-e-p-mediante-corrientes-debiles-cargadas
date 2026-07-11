# ================================================================
#  Código para el cáculo de las secciones eficaces detectando el nucleón
# ================================================================

import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','grid','notebook'])

# ================================================================
#  CONSTANTES 
# ================================================================

I = 1j

Mp  = 938.27203
MN  = 938.27203
MN2 = 881568.315821

Mpi  = 138.0389867
Mpi2 = 19054.761849

electronmass = 0.511
muonmass     = 105.658369

Pi  = 3.141592654
hbc = 197.3270

M_W     = 80385.0
M_Z     = 91187.6
Cabibbo = 0.974
G_Fermi = 1.16637e-11
Alpha   = 0.007297353
Sin2W   = 0.23122
QWeak   = 0.07512
gA      = 1.2695   
fpi     = 93

gmunu = np.array([
    [ 1.0,  0.0,  0.0,  0.0],
    [ 0.0, -1.0,  0.0,  0.0],
    [ 0.0,  0.0, -1.0,  0.0],
    [ 0.0,  0.0,  0.0, -1.0]
])

# ================================================================
#  FORM FACTORS
# ================================================================

def formfactors(Qsq, param):
    QsqGeV = Qsq / 1e6
    xmu = 2 * MN
    tau = Qsq / (xmu * xmu)

    MA = 1.032
    MA2 = MA * MA
    MV = 0.843
    MV2 = MV * MV

    DipV = 1.0 / (1.0 + QsqGeV / MV2)**2
    DipA = 1.0 / (1.0 + QsqGeV / MA2)**2

    rmup = 2.79285
    rmun = -1.91304

    # ============================================================
    # GALSTER
    # ============================================================
    if param == 1:
        GEp = DipV
        GMp = rmup * GEp
        GEn = -rmun * tau * GEp / (1.0 + 5.6 * tau)
        GMn = rmun * GEp

    # ============================================================
    # KELLY
    # ============================================================
    elif param == 3:
        a_m_neutron = [1, 2.33]
        a_m_proton = [1, 0.12]
        a_e_proton = [1, -0.24]
        b_e_proton = [10.98, 12.82, 21.97]
        b_m_proton = [10.97, 18.86, 6.55]
        b_m_neutron = [14.72, 24.2, 84.1]
        Omega = 0.71
        A = 1.7
        B = 3.3

        def G(a, b, t):
            num = sum(a[k] * t**k for k in range(len(a)))
            den = 1 + sum(b[k] * t**(k+1) for k in range(len(b)))
            return num / den

        GD = 1 / (1 + Qsq / Omega**2)**2
        GEp = G(a_e_proton, b_e_proton, tau)
        GMp = G(a_m_proton, b_m_proton, tau)
        GEn = (A * tau) / (1 + B * tau) * GD
        GMn = G(a_m_neutron, b_m_neutron, tau)

    # ============================================================
    # Arrington–Sick 
    # ============================================================
    elif param == 4:

        
        b_GEp  = [3.440, -0.178, -1.212, 1.176, -0.284]
        b_GMp  = [3.173, -0.314, -1.165, 5.619, -1.087]
        b_GEn  = [0.977, -20.82, 22.02]
        b_GMn  = [3.297, -0.258, 0.001]

        Q = np.sqrt(max(QsqGeV, 0.0))

        # Continued Fraction
        def CF(Q2, b):
            val = 1.0
            for bi in reversed(b):
                val = 1.0 + bi * Q2 / val
            return 1.0 / val

        # Proton
        GEp = CF(QsqGeV, b_GEp)
        GMp = rmup * CF(QsqGeV, b_GMp)

        # Neutron
        GEn = 0.484 * QsqGeV * CF(QsqGeV, b_GEn)
        GMn = rmun * CF(QsqGeV, b_GMn)
        
    F1p = (GEp + tau * GMp) / (1.0 + tau)
    F2p = (GMp - GEp) / (1.0 + tau)
    F1n = (GEn + tau * GMn) / (1.0 + tau)
    F2n = (GMn - GEn) / (1.0 + tau)

    GA = -gA * DipA
    GP = GA * 4 * MN / (Qsq + Mpi2)

    return F1p, F2p, F1n, F2n, GA, GP, 0.0, 0.0, 0.0

# ================================================================
#  TENSOR LEPTÓNICO
# ================================================================

def Lmunu(process, Helicity, Ki, Kf):
    Ki = np.asarray(Ki, dtype=float)
    Kf = np.asarray(Kf, dtype=float)

    Ki_Kf = Ki[0]*Kf[0] - Ki[1]*Kf[1] - Ki[2]*Kf[2] - Ki[3]*Kf[3]

    L = np.zeros((4, 4), dtype=complex)
    for i in range(4):
        for j in range(4):
            L[i, j] = 0.5 * (Ki[i]*Kf[j] + Ki[j]*Kf[i] - gmunu[i, j]*Ki_Kf)

    if process == 0:
        return L

    LA = np.zeros((4, 4), dtype=complex)

    LA[0,1] = (-I)*(Ki[2]*Kf[3] - Kf[2]*Ki[3])
    LA[0,2] = (-I)*(-Ki[1]*Kf[3] + Kf[1]*Ki[3])
    LA[0,3] = (-I)*(Ki[1]*Kf[2] - Kf[1]*Ki[2])

    LA[1,0] = -LA[0,1]
    LA[1,2] = (I)*(Ki[0]*Kf[3] - Kf[0]*Ki[3])
    LA[1,3] = (I)*(-Ki[0]*Kf[2] + Kf[0]*Ki[2])

    LA[2,0] = -LA[0,2]
    LA[2,1] = -LA[1,2]
    LA[2,3] = (I)*(Ki[0]*Kf[1] - Kf[0]*Ki[1])

    LA[3,0] = -LA[0,3]
    LA[3,1] = -LA[1,3]
    LA[3,2] = -LA[2,3]

    return 4.0 * L + 2.0 * (-Helicity) * LA

# ================================================================
#  TENSOR HADRÓNICO
# ================================================================

def Hmunu(process, nucleon, Pm, PN,param):
    Pm = np.asarray(Pm, dtype=float)
    PN = np.asarray(PN, dtype=float)

    Q = PN - Pm
    Qsq = - (Q[0]**2 - Q[1]**2 - Q[2]**2 - Q[3]**2)

    F1p, F2p, F1n, F2n, GA, FP, F1s, F2s, GAs = formfactors(Qsq,param)

    H = np.zeros((4,4), dtype=complex)

    if process == 0:
        PmPN = Pm[0]*PN[0] - Pm[1]*PN[1] - Pm[2]*PN[2] - Pm[3]*PN[3]

        if nucleon == 1:
            F1, F2 = F1p, F2p
        else:
            F1, F2 = F1n, F2n

        for i in range(4):
            for j in range(4):
                term1 = (F1+F2)**2 * (Pm[i]*PN[j] + Pm[j]*PN[i] + (MN2-PmPN)*gmunu[i,j])
                term2 = ((F2/(2*MN))**2 * (PmPN + MN2) - F2*(F1+F2)) * (Pm[i]+PN[i])*(Pm[j]+PN[j])
                H[i,j] = (term1 + term2) / (2*MN2)

        return H

    if process == 1:
        Pm_c = np.array([Pm[0], -Pm[1], -Pm[2], -Pm[3]])
        Q_c  = np.array([Q[0], -Q[1], -Q[2], -Q[3]])

        F1V = F1p - F1n
        F2V = F2p - F2n

        GA2 = GA*GA
        FP2 = FP*FP
        tau = Qsq / (4*MN2)

        W1 = tau*((F1V+F2V)**2 + GA2) + GA2
        W2 = F1V**2 + tau*F2V**2 + GA2
        W3 = 2*GA*(F1V+F2V)
        W4 = (F2V**2)/4*(tau-1) - F1V*F2V/2 - FP*GA*MN + MN2*tau*FP2
        W5 = W2

        Hs = np.zeros((4,4), dtype=float)
        for i in range(4):
            for j in range(4):
                Hs[i,j] = (
                    -W1*MN2*gmunu[i,j]
                    + W2*Pm[i]*Pm[j]
                    + W4*Q[i]*Q[j]
                    + W5/2*(Pm[i]*Q[j] + Q[i]*Pm[j])
                )

        Ha = np.zeros((4,4), dtype=float)

        Ha[0,1] = (-1*Pm_c[2]*Q_c[3] + 1*Pm_c[3]*Q_c[2]) * W3/2
        Ha[1,0] = -Ha[0,1]

        Ha[0,2] = (1*Pm_c[1]*Q_c[3] + -1*Pm_c[3]*Q_c[1]) * W3/2
        Ha[2,0] = -Ha[0,2]

        Ha[0,3] = (-1*Pm_c[1]*Q_c[2] + 1*Pm_c[2]*Q_c[1]) * W3/2
        Ha[3,0] = -Ha[0,3]

        Ha[1,2] = (-1*Pm_c[0]*Q_c[3] + 1*Pm_c[3]*Q_c[0]) * W3/2
        Ha[2,1] = -Ha[1,2]

        Ha[1,3] = (1*Pm_c[0]*Q_c[2] + -1*Pm_c[2]*Q_c[0]) * W3/2
        Ha[3,1] = -Ha[1,3]

        Ha[2,3] = (-1*Pm_c[0]*Q_c[1] + 1*Pm_c[1]*Q_c[0]) * W3/2
        Ha[3,2] = -Ha[2,3]

        return (Hs + I*Ha) / MN2

    if process == 2:
        Pm_c = np.array([Pm[0], -Pm[1], -Pm[2], -Pm[3]])
        PN_c = np.array([PN[0], -PN[1], -PN[2], -PN[3]])

        PmPN = Pm[0]*PN[0] - Pm[1]*PN[1] - Pm[2]*PN[2] - Pm[3]*PN[3]

        tau3 = 1 if nucleon == 1 else -1

        wF1 = (0.5-Sin2W)*(F1p-F1n)*tau3 - Sin2W*(F1p+F1n) - 0.5*F1s
        wF2 = (0.5-Sin2W)*(F2p-F2n)*tau3 - Sin2W*(F2p+F2n) - 0.5*F2s
        wGA = 0.5*(tau3*GA + GAs)

        HVV = np.zeros((4,4), dtype=complex)
        HAA = np.zeros((4,4), dtype=complex)
        HVA = np.zeros((4,4), dtype=complex)

        for i in range(4):
            for j in range(4):
                termVV = (wF1+wF2)**2 * (
                    Pm[i]*PN[j] + Pm[j]*PN[i] + (MN2-PmPN)*gmunu[i,j]
                )
                termVV += ((wF2/(2*MN))**2*(PmPN+MN2) - wF2*(wF1+wF2)) * \
                          (Pm[i]+PN[i])*(Pm[j]+PN[j])
                HVV[i,j] = termVV / (2*MN2)

                tau = Qsq/(4*MN2)
                HAA[i,j] = (
                    -gmunu[i,j]*(tau+1)
                    + Pm[i]*Pm[j]/MN2
                    + (Pm[i]*Q[j] + Pm[j]*Q[i])/(2*MN2)
                ) * (wGA*wGA)

        temp = 2*I*wGA*(wF1+wF2)/(2*MN2)

        HVA[0,1] = temp*(-1*Pm_c[2]*PN_c[3] + 1*Pm_c[3]*PN_c[2])
        HVA[1,0] = -HVA[0,1]

        HVA[0,2] = temp*(1*Pm_c[1]*PN_c[3] + -1*Pm_c[3]*PN_c[1])
        HVA[2,0] = -HVA[0,2]

        HVA[0,3] = temp*(-1*Pm_c[1]*PN_c[2] + 1*Pm_c[2]*PN_c[1])
        HVA[3,0] = -HVA[0,3]

        HVA[1,2] = temp*(-1*Pm_c[0]*PN_c[3] + 1*Pm_c[3]*PN_c[0])
        HVA[2,1] = -HVA[1,2]

        HVA[1,3] = temp*(1*Pm_c[0]*PN_c[2] + -1*Pm_c[2]*PN_c[0])
        HVA[3,1] = -HVA[1,3]

        HVA[2,3] = temp*(-1*Pm_c[0]*PN_c[1] + 1*Pm_c[1]*PN_c[0])
        HVA[3,2] = -HVA[2,3]

        return HVV + HAA + HVA

# ================================================================
#  CINEMÁTICA
# ================================================================

def main():
    with open("Input.txt") as f:
        process  = int(f.readline())
        nucleon  = int(f.readline())
        Ei       = float(f.readline())
        Helicity = int(f.readline())
        param = int(f.readline())

    mf = electronmass
    if process == 1:
        mf = 0.001

    ultrarel = 1

    P = np.array([MN, 0, 0, 0], dtype=float)

    ki = Ei
    Ki = np.array([Ei, 0, 0, ki], dtype=float)
    Ki_c = np.array([Ei, 0, 0, -ki], dtype=float)

    out = open(f"Ei{int(Ei)}.out", "w")

    angles = []
    dsig_costh = []
    omega_vals = []
    dsig_dQ2_vals = []
    dsig_omega_vals = []

    for thetaf_deg in range(1, 181, 2):
        thetaf = thetaf_deg * Pi / 180
        costheta_f = np.cos(thetaf)
        sintheta_f = np.sin(thetaf)

        Ef = Ei / (1 + Ei/MN * (1 - costheta_f))

        if ultrarel == 1:
            mi = 0.0
            AA = Ei*MN + (mi*mi + mf*mf)/2
            BB = Ei + MN
            CC = ki * costheta_f

            determ = (AA*BB)**2 - (AA*AA + CC*CC*mf*mf)*(BB*BB - CC*CC)
            if determ < 0:
                continue

            if costheta_f < 0:
                Ef = (AA*BB - np.sqrt(determ)) / (BB*BB - CC*CC)
            else:
                Ef = (AA*BB + np.sqrt(determ)) / (BB*BB - CC*CC)

            if Ef >= Ei or Ef < 0:
                continue
        kf = np.sqrt(Ef*Ef - mf*mf)
        Kf = np.array([Ef, kf*sintheta_f, 0.0, kf*costheta_f], dtype=float)
        Kf_c = np.array([Kf[0], -Kf[1], -Kf[2], -Kf[3]], dtype=float)

        Q = Ki - Kf
        if Q[0] < 0:
            continue

        Q2 = Q[0]**2 - Q[1]**2 - Q[2]**2 - Q[3]**2
        if Q2 >= 0:
            continue

        PN = P + Q
        if PN[0] < MN:
            continue

        Lmn = Lmunu(process, Helicity, Ki_c, Kf_c)
        Hmn = Hmunu(process, nucleon, P, PN,param)

        LmnHmn = np.sum(Lmn * Hmn)

        if process == 2:
            FF = (G_Fermi/np.sqrt(2.0))**2
        elif process == 1:
            FF = (G_Fermi*Cabibbo/np.sqrt(2.0))**2
        else:
            FF = (4*Pi*Alpha/Q2)**2

        Ki_P = Ki[0]*P[0] - Ki[1]*P[1] - Ki[2]*P[2] - Ki[3]*P[3]
        KK = MN2 / (abs(Ki_P) * Kf[0] * PN[0]) * FF

        frec = 1 + Kf[0] * (kf - ki*costheta_f) / (kf * PN[0])

        d2sig = KK * (kf/(2*Pi))**2 / abs(frec) * abs(LmnHmn)

        dsig = 2*Pi * hbc**2 * d2sig

        # Jacobiano correcto dω/d(cosθe) = Ef²/MN
        # derivado de Ef = Ei·MN / (MN + Ei·(1-cosθe))
        jac_e = Ef**2 / MN
        dsig_omega = dsig / jac_e

        jacQ2 = 2*Ei*Ef * (1 + Ef/MN * (costheta_f - 1))
        dsigQ2 = dsig / abs(jacQ2)

        angles.append(thetaf_deg)
        dsig_costh.append(dsig)
        omega_vals.append(Q[0])
        dsig_dQ2_vals.append(dsigQ2)
        dsig_omega_vals.append(dsig_omega)

        out.write(f"{Q[0]} {Kf[0]} {PN[0]} {-Q2/1e6} {dsigQ2} {thetaf_deg} {dsig}\n")

    out.close()

    # ============================
    #  PLOTS SECCIONES EFICACES
    # ============================

    angles_arr = np.array(angles)
    dsig_costh_arr = np.array(dsig_costh)
    omega_arr = np.array(omega_vals)
    dsig_dQ2_arr = np.array(dsig_dQ2_vals)
    dsigma_omega_arr = np.array(dsig_omega_vals)
    
    # guardamos variables para plotear aparte
    #np.save("datos/comparacion_energia/angulos_lepton_500.npy",angles_arr)
    #np.save("datos/comparacion_energia/sigma_lepton_500.npy",dsig_costh_arr)
    #np.save("datos/comparacion_energia/sigma_energia_lepton_500.npy",dsigma_omega_arr)
    #np.save("datos/comparacion_energia/omega_lepton_500.npy",omega_arr) 


    plt.figure(figsize=(14,6))
    
    plt.subplot(1, 2, 1)
    plt.plot(angles_arr, dsig_costh_arr, '-', color='blue')
    plt.yscale('log')
    plt.xlabel(r"$θ_\nu$ (grados)")
    plt.ylabel("dσ/dcosθ  (fm²)")
    #plt.title("Sección eficaz diferencial para el lepton")
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(omega_arr, dsigma_omega_arr, '-', color='red')
    plt.yscale('log')
    plt.xlabel("ω (MeV)")
    plt.ylabel("dσ/dω  (fm²/MeV)")
    #plt.title("Sección eficaz diferencial para el lepton")
    plt.grid(True)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
    