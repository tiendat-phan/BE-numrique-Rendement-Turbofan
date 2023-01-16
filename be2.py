import numpy as np
import matplotlib.pyplot as plt

# Donnees initiales 

M0 =	0.78
P0	= 22700
T0	= 217
eps_e = 0.98
n_c	= 0.9
n_f	= 0.92
n_comp = 0.99
eps_cc = 0.95
n_m	= 0.98  # rendement de l'arbre
n_thp = 0.89
n_tbp = 0.9
eps_tuy = 0.98
Tt4 = 1600
pi_c = 40
pi_f = 1.4

# Question 2: l'evolution du rendement thermique 

pi_cph = 22

r = 287
r_etoil	= 291.6
gamma	= 1.4
gamma_etoil	= 1.33
Pk	= 42800000

lda = 11
Fn = 21000

def returnRendement(pi_f, Tt4, pi_c, lda):
    #Start
    cp = (gamma/(gamma-1))*r
    
    cp_hot = (gamma_etoil/(gamma_etoil-1))*r_etoil
  
    # print('\n------- Inlet-------\n')
    Pt0 = P0*np.power((1+0.5*(gamma-1)*M0*M0),gamma/(gamma-1))
   
    Tt0 = T0*np.power((1+0.5*(gamma-1)*M0*M0),1)
    
    V0 = M0*np.sqrt(r*gamma*T0)

    Tt2 = Tt0
   
    Pt2 = eps_e*Pt0
  
    # print('\n------- Rotor oulet-------\n')
    Pt21 = pi_f*Pt2
    
    a1 = (gamma-1)/(gamma*n_f)
    Tt21 = Tt2*np.power(pi_f,a1)
    
    # print('\n------- LPC oulet-------\n')
    pi_25 = pi_c/(pi_cph*pi_f)

    Pt25 = pi_25*Pt21
    
    a2 = (gamma-1)/(gamma*n_c)
    Tt25 = Tt21*np.power(pi_25,a2)
   
    # print('\n------- HPC oulet-------\n')

    Pt3 = pi_c*Pt2 
    a3 = (gamma-1)/(gamma*n_c)
    Tt3 =  Tt25*np.power(pi_cph,a3)
    
    # print('\n------- Combustor oulet-------\n')

    Pt4 = eps_cc*Pt3   
    alpha = (cp_hot*Tt4 - cp*Tt3)/(n_comp*Pk - cp_hot*Tt4) 
    

    # print('\n------- HPT oulet-------\n')
    wu_hpt = cp*(Tt3 - Tt25)
    wut = -wu_hpt


    Tt45 = Tt4 + wut/(n_m*(1+alpha)*cp_hot)
    Pt45 = Pt4*np.power(Tt45/Tt4, gamma_etoil/((gamma_etoil-1)*n_thp))
  
    # print('\n------- LPT oulet-------\n')

    wu_lpt = cp*(Tt25-Tt21) + cp*(1+lda)*(Tt21-Tt2)
    Tt5 = Tt45 - wu_lpt/(n_m*cp_hot*(1+alpha))
    Pt5 = Pt45*np.power(Tt5/Tt45, gamma_etoil/((gamma_etoil-1)*n_tbp))

 
    # print('\n------- Nozzle oulet core-------\n')
    Pt9 = eps_tuy*Pt5
    Tt9 = Tt5
    P9 = P0 # tuyere adaptee

    M9 = np.sqrt(2*(np.power(Pt9/P9,(gamma_etoil-1)/gamma_etoil)-1)/(gamma_etoil-1))
    T9 = Tt9/(1+ 0.5*(gamma_etoil-1)*M9*M9)
    V9 = M9*np.sqrt(r_etoil*gamma_etoil*T9)
    F9 = ((1+alpha)*V9 - V0)
    F9_total = ((1+alpha)*V9 - V0)*1/(lda+1)
   

    # print('\n------- Nozzle oulet bypass-------\n')
    Pt17 = Pt21
    Tt17 = Tt21
    Tt19 = Tt17
    Pt19 = eps_tuy*Pt17 #tuyere adapte
    P19 = P0
    M19 = np.sqrt(2*(np.power(Pt19/P19,(gamma-1)/gamma)-1)/(gamma-1))
    T19 = Tt19/(1+ 0.5*(gamma-1)*M19*M19)
    V19 = M19*np.sqrt(r*gamma*T19)
    F19 = (V19 - V0)
    F19_total = (V19 - V0)*lda/(lda+1) 

 
    # print('\n------- Profulsif cooefficents-------\n')

    wchim = alpha*Pk/(lda+1)

    wcy9 = 0.5*((1+alpha)*V9*V9 - V0*V0)
    wcy9_total = (1/(lda+1))*0.5*((1+alpha)*V9*V9 - V0*V0)

    wcy19 = 0.5*(V19*V19 - V0*V0)
    wcy19_total = (lda/(lda+1))*0.5*(V19*V19 - V0*V0)

    wpr19 = F19*V0
    wpr9 = F9*V0

    wpr19_total = F19_total*V0
    wpr9_total = F9_total*V0

    n_19 = wpr19/wcy19
    n_9 = wpr9/wcy9


    Ftot = F19_total + F9_total

    n_th = (wcy19_total+wcy9_total)/wchim
    n_pr = (wpr19_total+wpr9_total)/(wcy19_total+wcy9_total)
    n_global = n_th*n_pr

    # print('\n------- Sizing-------\n')
    k_total =  lda*(V19-V0)/(lda+1) + ((1+alpha)*V9-V0)/(1+lda)

    m_total = Fn/k_total
    m_core = (1/(lda+1))*m_total
    m_bypass = m_total - m_core

   
    rho0 = P0/(r*T0)
    S_2 = m_total/(V0*rho0)

    D_fan = np.sqrt(4*S_2/np.pi)

    return n_th,n_global 

   

if __name__ == '__main__':
    
    #Question 2: l'evolution du rendement thermique en fonction du taux de compression
    pi_c = np.arange(30,70,1).tolist()
    a = np.array(pi_c)
    
    lamda = np.arange(1,30,1).tolist()
    lamda = np.array(lamda)


    nth1, nglo1 = returnRendement(pi_f = 1, Tt4 = 1600, pi_c = 40, lda = lamda)
    nth2, nglo2 = returnRendement(pi_f = 1.2, Tt4 = 1600, pi_c = 40, lda = lamda)
    nth3, nglo3 = returnRendement(pi_f = 1.3, Tt4 = 1600, pi_c = 40, lda = lamda)
    nth4, nglo4 = returnRendement(pi_f = 1.4, Tt4 = 1600, pi_c = 40, lda = lamda)
    nth5, nglo5 = returnRendement(pi_f = 1.6, Tt4 = 1600, pi_c = 40, lda = lamda)
    nth6, nglo6 = returnRendement(pi_f = 2, Tt4 = 1600, pi_c = 40, lda = lamda)

    
    # print(round(nth1,5))
    plt.plot(lamda, nglo1, label=r'$\pi_f = 1 $')
    plt.plot(lamda, nglo2, label=r'$\pi_f = 1.2$')
    plt.plot(lamda, nglo3, label=r'$\pi_f = 1.4$')
    plt.plot(lamda, nglo4, label=r'$\pi_f = 1.4$')
    plt.plot(lamda, nglo5, label=r'$\pi_f = 1.8$')
    plt.plot(lamda, nglo6, label=r'$\pi_f = 2$')
    # plt.plot(a, nglo1, label=r'$T_{t4} = 1600 C, global$', color= 'g' )
    # plt.title(r'Le rendement thermique - le OPR')
    plt.xlabel(r'$\lambda$')   
    plt.ylabel(r'$n_{gl}$')
    plt.legend()
    plt.show()

