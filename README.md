# Dispersion-curves-HASKELL

To launch the functions follow the steps below (the help of the functions will also help):

##### Set:
n       = number of layers

Vp,Vs   = Velocities as numpy array

rho     = density as an array

D       = depth of each interface as an array

freqmax = maximum frequency 
##### To set the maximum frequency if unknown:
freqmax=min(Vs)/(depth[len(depth)-2]-depth[len(depth)-3]) 


##### Once the initialization parameters are fixed, use the "dispersion_curves_gen" to generate the curves.
 ________________
 Here an exemple:
 ________________ 

import numpy as np

from funzione_heskell import *

n=1 

Vp=[600, 1400]

Vp=np.array(Vp)

Vs=[300, 600]

Vs=np.array(Vs)

rho=[1.82, 1.91]

rho=np.array(rho)

D=[8.0, 0]

D=np.array(D)

fmax=50

dir_principale=  ............. #where to find the functions

directoy_vphase= ............. #where to find the vphase array

directory_freq=  ............. #where to find the frequency array

nomefilef='f1'

nomefilev='v1'

dispersion_curves_gen (dir_principale,n,Vp,Vs,rho,D,fmax,directory_freq,directoy_vphase,nomefilef,nomefilev, PLOT=True) 

![image](https://user-images.githubusercontent.com/108676675/200411255-301a2ebd-3cec-47b6-83f2-454d1cc86755.png)

