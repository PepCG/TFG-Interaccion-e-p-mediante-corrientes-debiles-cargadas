# ================================================================
#  Sección eficaz detectando el nucleón
#  Traducción del código C++ de Raúl González
# ================================================================

import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','grid','notebook'])

# ================================================================
#  CONSTANTES  Y MÉTRICA
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
gA      = 1.265
fpi     = 93   

#  Parámetros de la parametrización z-expansion
zexp_tc = 9.0 * 0.134**2   
zexp_t0 = -0.50            
zexp_coeffs = np.array([
    0.72115656, -1.70104640,  0.26324902,
    1.53433681,  0.01061114, -1.49893610,
    0.67062898
])
zexp_gA_beta = 1.2754
zexp_a1a2_mean = np.array([-1.70104640, 0.26324902])
zexp_a1a2_cov = np.array([
    [ 0.00175929, -0.00294651],
    [-0.00294651,  0.00807158]
])

# métrica
gmunu = np.array([
    [ 1.0,  0.0,  0.0,  0.0],
    [ 0.0, -1.0,  0.0,  0.0],
    [ 0.0,  0.0, -1.0,  0.0],
    [ 0.0,  0.0,  0.0, -1.0]
])

# ================================================================
#  FORM FACTORS
# ================================================================

def z_expansion_FA(QsqGeV, tc=zexp_tc, t0=zexp_t0, coeffs=zexp_coeffs):
    sqrt_tcQ2 = np.sqrt(tc + QsqGeV)
    sqrt_tct0 = np.sqrt(tc - t0)
    z = (sqrt_tcQ2 - sqrt_tct0) / (sqrt_tcQ2 + sqrt_tct0)

    FA = 0.0
    for k, ak in enumerate(coeffs):
        FA += ak * z**k
    return FA


def formfactors(Qsq, param, axial_param=1):
    QsqGeV = Qsq / 1e6
    xmu = 2 * MN
    tau = Qsq / (xmu * xmu)

    MA = 1.23
    MA2 = MA * MA
    MV = 0.843
    MV2 = MV * MV

    # forma dipolar 
    DipV = 1.0 / (1.0 + QsqGeV / MV2)**2
    DipA = 1.0 / (1.0 + QsqGeV / MA2)**2
    
    # momentos magnéticos
    rmup = 2.79285
    rmun = -1.91304

    # ============================================================
    # 1. GALSTER
    # ============================================================
    if param == 1:
        GEp = DipV
        GMp = rmup * GEp
        GEn = -rmun * tau * GEp / (1.0 + 5.6 * tau)
        GMn = rmun * GEp

    # ============================================================
    # 2. GALSTER Q<< (no utilizado al final)
    # ============================================================
    elif param == 2:
        GD = 1 / (1 + 4.97 * tau)**2
        epsilon = 1 / (1 + 4.97 * tau)
        GEp = GD
        GMp = rmup * GD
        GEn = -rmun * tau * GD * epsilon
        GMn = rmun * GD

    # ============================================================
    # 3. KELLY
    # ============================================================
    elif param == 3:
     # coeficientes de la tabla del artículo original
     a_m_neutron = [1, 2.33]
     a_m_proton = [1, 0.12]
     a_e_proton = [1, -0.24]
     b_e_proton = [10.98, 12.82, 21.97]
     b_m_proton = [10.97, 18.86, 6.55]
     b_m_neutron = [14.72, 24.2, 84.1]
     A = 1.7
     B = 3.3
     
     # función de la parametrización
     def G(a, b, t):
        num = sum(a[k] * t**k for k in range(len(a)))
        den = 1 + sum(b[k] * t**(k+1) for k in range(len(b)))
        return num / den

     GEp = G(a_e_proton, b_e_proton, tau)
     GMp = rmup * G(a_m_proton, b_m_proton, tau)
     GEn = (A * tau) / (1 + B * tau) * DipV
     GMn = rmun * G(a_m_neutron, b_m_neutron, tau)

    # ============================================================
    # 4. ARRINGTON-SICK
    # ============================================================
    elif param == 4:

        # Coeficientes de la tabla del artículo orginal
        b_GEp  = [3.440, -0.178, -1.212, 1.176, -0.284]
        b_GMp  = [3.173, -0.314, -1.165, 5.619, -1.087]
        b_GEn  = [0.977, -20.82, 22.02]
        b_GMn  = [3.297, -0.258, 0.001]

        Q = np.sqrt(max(QsqGeV, 0.0))

        # funcion de fraccion continua CF
        def CF(Q2, b):
            val = 1.0
            for bi in reversed(b):
                val = 1.0 + bi * Q2 / val
            return 1.0 / val

        GEp = CF(QsqGeV, b_GEp)
        GMp = rmup * CF(QsqGeV, b_GMp)
        GEn = 0.484 * QsqGeV * CF(QsqGeV, b_GEn)
        GMn = rmun * CF(QsqGeV, b_GMn)

    # ============================================================
    # 5. GK
    # ============================================================
    elif param == 5:
        # parametros
        g_rho = 0.377      
        K_rho = 6.62       
        g_omega = 0.411    
        K_omega = 0.163    

        kappa_v = 3.706    
        kappa_s = -0.12    

        m_rho = 0.776
        m_omega = 0.784

        Lambda1 = 0.795
        Lambda2 = 2.27
        LambdaQCD = 0.29
        
        # funcion
        Q2 = Qsq / 1e6
        logterm = np.log((Lambda2**2 + Q2) / LambdaQCD**2)
        log0 = np.log(Lambda2**2 / LambdaQCD**2)
        Qtilde2 = Q2 * log0 / logterm

        pole1 = Lambda1**2 / (Lambda1**2 + Qtilde2)
        pole2 = Lambda2**2 / (Lambda2**2 + Qtilde2)

        # factores de Dirac-Pauli
        F1 = pole1 * pole2        
        F2 = pole1 * pole2**2     

        F1V = (g_rho*(m_rho**2/(m_rho**2+Q2)) + (1-g_rho)) * F1
        F2V = (K_rho*g_rho*(m_rho**2/(m_rho**2+Q2)) + (kappa_v - K_rho*g_rho)) * F2

        F1S = (g_omega*(m_omega**2/(m_omega**2+Q2)) + (1-g_omega)) * F1
        F2S = (K_omega*g_omega*(m_omega**2/(m_omega**2+Q2)) + (kappa_s - K_omega*g_omega)) * F2

        F1p = 0.5*(F1S + F1V)
        F1n = 0.5*(F1S - F1V)
        F2p = 0.5*(F2S + F2V)
        F2n = 0.5*(F2S - F2V)

        if axial_param == 2:
            GA = -z_expansion_FA(QsqGeV)
        else:
            GA = -gA * DipA
        GP = GA * 4 * MN / (Qsq + Mpi2)

        return F1p, F2p, F1n, F2n, GA, GP, 0.0, 0.0, 0.0

    # ============================================================
    # 6. GKex
    # ============================================================
    elif param == 6:
        # parametros
        m_rho = 0.776
        m_omega = 0.784
        m_phi = 1.019

        g_rho, K_rho = 0.4466, 4.3472      
        g_omega, K_omega = 0.4713, 21.762  
        g_phi, K_phi = -0.8461, 11.849     
        mu_phi = 1.1498

        kappa_v, kappa_s = 3.706, -0.12    

        Lambda1 = 0.9006    
        LambdaD = 1.7038     
        Lambda2 = 1.1336     
        LambdaQCD = 0.0312

        # funcion
        Q2 = Qsq / 1e6
        logterm = np.log((Lambda2**2 + Q2) / LambdaQCD**2)
        log0 = np.log(Lambda2**2 / LambdaQCD**2)
        Qtilde2 = Q2 * log0 / logterm   # Q2 efectivo de pQCD

        pole2 = Lambda2**2 / (Lambda2**2 + Qtilde2)

        # factores de Dirac-Pauli
        F1m = (Lambda1**2 / (Lambda1**2 + Qtilde2)) * pole2
        F2m = (Lambda1**2 / (Lambda1**2 + Qtilde2))**2 * pole2
        F1D = (LambdaD**2 / (LambdaD**2 + Qtilde2)) * pole2
        F2D = (LambdaD**2 / (LambdaD**2 + Qtilde2))**2 * pole2

        F1phi = F1m * (Q2 / (Lambda1**2 + Q2))**1.5
        F2phi = F2m * (Lambda1**2/mu_phi**2 * (Q2 + mu_phi**2)/(Lambda1**2 + Q2))**1.5

        F1V = g_rho*(m_rho**2/(m_rho**2+Q2))*F1m + (1-g_rho)*F1D
        F2V = K_rho*g_rho*(m_rho**2/(m_rho**2+Q2))*F2m + (kappa_v - K_rho*g_rho)*F2D

        F1S = g_omega*(m_omega**2/(m_omega**2+Q2))*F1m \
            + g_phi*(m_phi**2/(m_phi**2+Q2))*F1phi \
            + (1 - g_omega)*F1D

        F2S = K_omega*g_omega*(m_omega**2/(m_omega**2+Q2))*F2m \
            + K_phi*g_phi*(m_phi**2/(m_phi**2+Q2))*F2phi \
            + (kappa_s - K_omega*g_omega - K_phi*g_phi)*F2D

        F1p = 0.5*(F1S + F1V)
        F1n = 0.5*(F1S - F1V)
        F2p = 0.5*(F2S + F2V)
        F2n = 0.5*(F2S - F2V)

        if axial_param == 2:
            GA = -z_expansion_FA(QsqGeV)
        else:
            GA = -gA * DipA
        GP = GA * 4 * MN / (Qsq + Mpi2)

        return F1p, F2p, F1n, F2n, GA, GP, 0.0, 0.0, 0.0

    # ============================================================
    # CONVERSIÓN A F1, F2 
    # ============================================================
    F1p = (GEp + tau * GMp) / (1.0 + tau)
    F2p = (GMp - GEp) / (1.0 + tau)
    F1n = (GEn + tau * GMn) / (1.0 + tau)
    F2n = (GMn - GEn) / (1.0 + tau)

    if axial_param == 2:
        GA = -z_expansion_FA(QsqGeV)
    else:
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

    if process == 0: # proceso EM
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

def Hmunu(process, nucleon, Pm, PN, param, axial_param=1):
    Pm = np.asarray(Pm, dtype=float)
    PN = np.asarray(PN, dtype=float)

    Q = PN - Pm
    Qsq = - (Q[0]**2 - Q[1]**2 - Q[2]**2 - Q[3]**2)

    F1p, F2p, F1n, F2n, GA, FP, F1s, F2s, GAs = formfactors(Qsq, param, axial_param)

    H = np.zeros((4,4), dtype=complex)

    # EM
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

    # CC
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

    # NC
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

# leemos el archivo Input.txt que nos permite seleccionar el tipo de proceso que queremos medir
def main():
    with open("Input.txt") as f:
        process  = int(f.readline())
        nucleon  = int(f.readline())
        Ei       = float(f.readline())
        Helicity = int(f.readline())
        param = int(f.readline())
        axial_param = int(f.readline())


    mf = electronmass
    if process == 1:
        mf = 0.001  # valor de masa arbitrario para el neutrino

    P = np.array([MN, 0, 0, 0], dtype=float)

    ki = Ei
    Ki = np.array([Ei, 0, 0, ki], dtype=float)

    out = open(f"Ei{int(Ei)}.out", "w")
    
    # creamos las variables vacías necesarias para el ploteo posterior
    angles = []
    dsig_vals = []
    pN_vals = []
    Q2_vals = []
    dsig_omega_vals = []
    omega_vals = []

    for thetaN_deg in range(1, 181, 2):

        thetaN = thetaN_deg * Pi / 180
        cosN = np.cos(thetaN)
        sinN = np.sin(thetaN)

        denom = (Ei + MN)**2 - (ki**2) * (cosN**2)
        if denom <= 0:
            continue

        pN = (2 * MN * (Ei + MN) * ki * cosN) / denom
        if pN <= 0:
            continue

        EN = np.sqrt(pN**2 + MN2)

        Ef = Ei + MN - EN
        if Ef <= 0:
            continue

        kf = np.sqrt(max(Ef**2 - mf**2, 0.0))

        PN = np.array([EN, pN*sinN, 0.0, pN*cosN], dtype=float)

        Q = P - PN
        Kf = Q - Ki
        Q2 = Q[0]**2 - Q[1]**2 - Q[2]**2 - Q[3]**2
        if Q2 >= 0:
            continue

        frec = (pN / EN) + (pN - ki*cosN) / Ef
        if frec == 0:
            continue

        Lmn = Lmunu(process, Helicity, Ki, Kf)
        Hmn = Hmunu(process, nucleon, P, PN, param, axial_param)

        LmnHmn = np.sum(Lmn * Hmn)

        if process == 2:
            FF = (G_Fermi/np.sqrt(2.0))**2
        elif process == 1:
            FF = (G_Fermi*Cabibbo/np.sqrt(2.0))**2
        else:
            FF = (4*Pi*Alpha/Q2)**2

        Ki_P = Ki[0]*P[0] - Ki[1]*P[1] - Ki[2]*P[2] - Ki[3]*P[3]
        Kpref = MN2 / (np.sqrt(Ki_P**2 - mf**2*MN**2) * EN * Ef) * FF

        d2sig = Kpref * (pN/(2*Pi))**2 / abs(frec) * abs(LmnHmn)
        dsig = hbc**2 * d2sig
        
        jac_e = Ef**2 / MN
        dsig_omega = dsig / jac_e

        angles.append(thetaN_deg)
        dsig_vals.append(dsig)
        pN_vals.append(pN)
        omega_vals.append(np.abs(Q[0]))
        Q2_vals.append(-Q2/1e6)
        dsig_omega_vals.append(dsig_omega)

        out.write(f"{thetaN_deg} {pN} {EN} {Ef} {Q2_vals[-1]} {dsig}\n")

    out.close()

    angles_arr = np.array(angles)
    dsig_arr = np.array(dsig_vals)
    pN_arr = np.array(pN_vals)
    Q2_arr = np.array(Q2_vals)
    dsigma_omega_arr = np.array(dsig_omega_vals)
    omega_arr = np.array(omega_vals)
    
    # guardamos variables para plotear aparte
    #np.save("datos/comparacion_ma/angulos_nucleon_ma1230.npy",angles_arr)
    #np.save("datos/comparacion_ma/sigma_nucleon_ma1230.npy",dsig_arr)
    #np.save("datos/comparacion_detectado/sigma_energia_proton.npy",dsigma_omega_arr)
    #np.save("datos/asimetria_positron/omega_electron_103.npy",omega_arr)
    #np.save("datos/comparacion_ma/q2_nucleon_ma1230.npy",Q2_arr) 
    #np.save("datos/asimetria_positron/sigma.npy",dsig_arr) 
    
    # ================================================================
    #  GRÁFICAS
    # ================================================================

    plt.figure(figsize=(14,6))
    
    plt.subplot(1, 2, 1)
    plt.plot(angles_arr, dsig_arr, '-', color='blue')
    plt.yscale('log')
    plt.xlabel("$θ_n$ (grados)")
    plt.ylabel("dσ/dcosθ  (fm²)")
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(omega_arr, dsigma_omega_arr, '-', color='red')
    plt.yscale('log')
    plt.xlabel("ω (MeV)")
    plt.ylabel("dσ/dω  (fm²/MeV)")
    plt.grid(True)

    plt.tight_layout()

    plt.figure(figsize=(8,6))
    plt.plot(Q2_arr, dsig_arr, color='green')
    plt.yscale('log')
    plt.xlabel(r"$Q^2$  (GeV$^2$)")
    plt.ylabel("dσ/dΩ_N  (fm²)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# ================================================================
#  EJECUCIÓN
# ================================================================

if __name__ == "__main__":
    main()
