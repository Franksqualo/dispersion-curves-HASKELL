import numpy as np
import os
import matplotlib.pyplot as plt
import cmath
from solidlayer import *
from fluidlayer import *
from math import sqrt


def dispersion_curves_gen (dir_principale,n,Vp,Vs,rho,D,fmax,directory_freq,directoy_vphase,nomefilef,nomefilev,PLOT=True):
    """Function to calculate dispersion curves for a 1D medium, with vertical heterogeneity.
    The method adopthed is the Haskell method.
    We suppose that: the medium is is totally elastic, monodimensional and the waves are plane waves.
    NB || The operator must take in account that the Haskell method assume that the curvature of the wavefronts
    is negligible, so that the distance source-receiver must be large enough or an higher frequency have to be chosen.
    
    Args:
        n (j)   : numero di strati considerando anche il semispazio
        Vp (j)  : vel delle onde P nello strato j in (m/s)
        Vs (j)  : vel delle onde S nello strato j in (m/s)
        rho (j) : densita' dello strato j in (g/cm^3)
        D (j)   : spessore dello strato j in (m)
        fmax (j): frequenza max da considerare in (Hz)
        
        dir_principale (j) : Storage funzioni solid-fluid layers
        directory_freq (j) : where to save the fobserved
        directoy_vphase (j): where to save the vphase observed
        nomefilef (j): how to save fobserved
        nomefilev (j): how to save vphase observed
    """
       
    os.chdir(dir_principale)
   
    # Inizializzazione degli altri parametri:
    
    # vsmin=min(Vs);
    vels=Vs[0]
    velp=Vp[0]
    if (Vs[0]==0):
        vels=velp
    TH=D[0]
    vp=Vp/vels
    vs=Vs/vels
    rho=rho/rho[0]
    d=D/TH
    k0=0.0001
    kf=(2*fmax*2*pi)/np.min(Vs)
    # kf=15
    dk=0.001
    nk=round((kf-k0)/dk)
    k=np.linspace(k0,kf,nk)
    k=k*TH
    print(k.shape)
    # Applicazione funzioni
    if (Vs[0]==0):
        T,croot=fluidlayer(n,k,d,rho,vp,vs,vels,TH,fmax)
    else:
        T,croot=solidlayer(n,k,d,rho,vp,vs,vels,TH,fmax)
        
    # --------------------------------------------------------------------------------- #  
    #                    Rappresentazione della curva di dispersione                    #    
    # --------------------------------------------------------------------------------- #  
    freq=1/T            
    vphase=croot*vels   
    
    if PLOT:    
        fig, ax = plt.subplots()
        ax.plot(freq, vphase, linewidth=2.0)
        ax.set_xlabel('f(Hz)')
        ax.set_ylabel('c(m/s)')
        plt.show()
    
    
    freqobs=freq
    vphaseobs=vphase
    
    os.chdir(directory_freq)
    np.save(nomefilef,freqobs)

    os.chdir (directoy_vphase)
    np.save(nomefilev,vphaseobs)
    
    os.chdir(dir_principale)
    
    # return freq,vphase
    