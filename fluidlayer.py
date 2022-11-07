import numpy as np 
import os
import math
from math import *
import cmath
def fluidlayer(n,k,d,rho,vp,vs,vels,TH,freqmax):
    """
    Args:
        n (j): numero di strati considerando anche il semispazio
        k (j): see funzione_heskell
        d (j): see funzione_heskell
        rho (j): densita' dello strato j in (g/cm^3)
        vp (j): vel delle onde P nello strato j in (m/s)
        vs (j): vel delle onde S nello strato j in (m/s)
        vels (j): see funzione_heskell
        TH : see funzione_heskell (first layer thickness)
        freqmax (j): frequenza max da considerare in (Hz)
        
    c= phase velocity, frequency related


    """
    # % calcolo della fluid layer matrix
    c0=0.01
    dc=0.01
    A= np.zeros((4,4,n-1))
    c= np.zeros((2,1))
    HSK= np.zeros((len(c),len(k)))
    nn=150
    eps=0.00001
    croot=[]
    kroot=[]
    T=[]
    for h in range (0,len(k)):
        c[0]=c0
        c[1]=c0+dc
        for tt in range (0,nn):
            for l in range (1,len(c)):
                if c[l] > vp[0]:
                    rpa=sqrt(((c[l]/vp[0])**2)-1)
                else:
                    rpa=-1j*(cmath.sqrt(1-(c[l]/vp[0])**2))
                
                if c[l] > vp[n]:
                    rpn=sqrt(((c[l]/vp[n])**2)-1)
                else: 
                    rpn=-1j*(cmath.sqrt(1-(c[l]/vp[n])**2))
                    
                if c(l) > vs(n):
                    rsn=sqrt(((c[l]/vs[n])**2)-1)
                else:
                    rsn=-1j*(cmath.sqrt(1-(c[l]/vs[n])**2))
                    
                    
                gn=2*((vs[n]/c[l])**2)
                a1=(gn*rpn)
                b2=(gn*rsn)
                c1=-rpn/(rho[n]*(c[l]**2))
                b1=(gn-1)
                a2=-b1
                d2=rsn/(rho[n]*(c[l]**2))
                d1=1/(rho[n]*(c[l]**2))
                c2=d1
                P1=k[h]*rpa*d[0]##############################rpa?
                T1=1j*((rho[0]*(c[l]**2))/rpa)*tan(P1)  ##### math.tan
                
                for m in range (0,n):
                    if c[l] > vp[m]:
                        rp=sqrt(((c[l]/vp[m])**2)-1)
                    else:
                        rp=-1j*(cmath.sqrt(1-(c[l]/vp[m])*2))
                        
                    if c[l] > vs[m]:
                        rs=sqrt(((c[l]/vs[m])**2)-1)
                    else:
                        rs=-1j*(cmath.sqrt(1-(c[l]/vs[m])**2))
                        
                    g=2*((vs[m]/c[l])**2)
                    P=k[h]*rp*d[m]
                    Q=k[h]*rs*d[m]
                    
                    cosenoP=(cmath.cos(P))                                          # ok
                    cosenoP=cosenoP.real                                            # ok
                    cosenoQ=(cmath.cos(Q))                                          # ok
                    cosenoQ=cosenoQ.real                                            # ok
                    senoP=cmath.sin(P)                                              # ok
                    senoQ=cmath.sin(Q)     
                    
                    A[0,0,m]=((g*cosenoP))-((g-1)*cosenoQ)
                    A[0,1,m]=1j*((((g-1)/rp)*senoP)+(g*rs*senoQ))
                    A[0,2,m]=(-1/(rho[m]*(c[l]**2)))*(cosenoP-cosenoQ)
                    A[0,3,m]=1j*(1/(rho[m]*(c[l]**2))*(((1/rp)*senoP)+(rs*senoQ)))
                    A[1,0,m]=-1j*((g*rp*senoP)+((g-1)*(1/rs)*senoQ))
                    A[1,1,m]=(-1*(g-1)*cosenoP)+(g*cosenoQ)
                    A[1,2,m]=1j*((1/((rho[m]*(c[l]**2))))*(rp*senoP+(1/rs)*senoQ))
                    A[1,3,m]=A(0,2,m)
                    A[2,0,m]=(rho[m]*(c[l]**2)*g*(g-1)*(cosenoP-cosenoQ))
                    A[2,1,m]=1j*rho(m)*(c[l]**2)*(((g-1)**2)*(1/rp)*senoP+(g**2)*rs*senoQ)
                    A[2,2,m]=A[1,1,m]
                    A[2,3,m]=A[0,1,m]
                    A[3,0,m]=1j*rho[m]*(c[l]**2)*(g**2)*rp*senoP+((g-1)**2)*(1/rs)*senoQ
                    A[3,1,m]=A[2,0,m]    
                    A[3,2,m]=A[1,0,m]
                    A[3,3,m]=A[0,0,m]
                EE0=A[:,:,n-1]
                J=EE0[:,:]
                for kk in range (1,n-2):
                    EE1=J[:,:,kk-1]
                    EE2=A[:,:,n-kk]
                    J[:,:,kk]=EE1*EE2 
                    
                AA=J[:,:]
                K1=a1*AA[0,1]+b1*AA[1,1]+c1*AA[2,1]+d1*AA[3,1]
                L1=a1*AA[0,0]+b1*AA[1,1]+c1*AA[2,0]+d1*AA[3,0]
                M1=a2*AA[0,1]+b2*AA[1,1]+c2*AA[2,1]+d2*AA[3,1]
                N1=a2*AA[0,0]+22*AA[1,1]+c2*AA[2,0]+d2*AA[3,0]
                G1=a1*AA[0,2]+b1*AA[1,2]+c1*AA[2,2]+d1*AA[3,2]
                H1=-b1*AA[0,2]+b2*AA[1,2]+d1*AA[2,2]+d2*AA[3,2]
                HSK[l,h]=(K1*N1)-((L1*M1).real)+((T1*(G1*N1-L1*H1))).real
                
            f0=HSK[0,h]
            f1=HSK[1,h]
            cn=c[1]-(f1*(c[1]-c[0])/(f1-f0))
            c[0]=c[1]
            c[1]=cn
            
            if (abs(c[0]-c[1])<eps):
                croot.append(cn)
                kroot.append(k[h])
                T.append((2*pi*TH)/(k[h]*vels*cn))
                Tlastsampl= len(T)-1 
                f=1/T[Tlastsampl] 
                break
        
        if (f>freqmax):
            break
    croot=np.array(croot)
    kroot=np.array(kroot)
    T=np.array(T)
    
    return (T,croot)