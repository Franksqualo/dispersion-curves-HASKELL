import numpy as np 
import os
import math
from math import *
import cmath

def solidlayer(n,k,d,rho,vp,vs,vels,TH,freqmax):
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
    # calcolo della solid layer matrix
    c0=0.01
    dc=0.01
    A=np.zeros((4,4,n),complex)

    c=np.zeros((2,1)) 
    HSK=np.zeros((len(c),len(k)))
    nn=150
    eps=0.000001
    croot=[]
    kroot=[]
    T=[]
    for h in range (0,len(k)):
        c[0]=c0
        c[1]=c0+dc
        for tt in range (0,nn):
            for l in range (0,len(c)):
                if c[l] > vp[n]:
                    rpn=sqrt(((c[l]/vp[n])**2)-1)                               #ok
                    # rpn=complex(np.round(rpn.real,4),np.round(rpn.imag,4))

                else:
                    rpn=-1j*(cmath.sqrt(1-(c[l]/vp[n])**2))                     #ok
                    # rpn=complex(np.round(rpn.real,4),np.round(rpn.imag,4))

                if c[l] > vs[n]:
                    rsn=sqrt(((c[l]/vs[n])**2)-1)                               #ok
                    # rsn=complex(np.round(rsn.real,4),np.round(rsn.imag,4))

                else:
                    rsn=-1j*(cmath.sqrt(1-(c[l]/vs[n])**2))                     #ok
                    # rsn=complex(np.round(rsn.real,4),np.round(rsn.imag,4))

                
                gn=2*((vs[n]/c[l])**2)     #ok
                a1=(gn*rpn)                #ok
                b2=(gn*rsn)                #ok
                c1=-rpn/(rho[n]*(c[l]**2)) #ok
                # c1=complex(np.round(c1.real,1),np.round(c1.imag,1)) #ok

                b1=(gn-1)                  #ok
                a2=-b1                     #ok
                d2=rsn/(rho[n]*(c[l]**2))  #ok
                # d2=complex(np.round(d2.real,1),np.round(d2.imag,1)) #ok

                d1=1/(rho[n]*(c[l]**2))    #ok
                # d1=complex(np.round(d1.real,1),np.round(d1.imag,1)) #ok

                c2=d1
                for m in range (0,n):
                    if c[l] > vp[m]:
                        rp=sqrt(((c[l]/vp[m])**2)-1)           ############# ok 
                        # rp=complex(np.round(rp.real,4),np.round(rp.imag,4))
                    else:
                        rp=-1j*(cmath.sqrt(1-(c[l]/vp[m])**2)) ############# ok 
                        # rp=complex(np.round(rp.real,4),np.round(rp.imag,4)) 
                        
                        # ** fixed **
                        #_______________________________________________
                        # PROBLEMA DI APPROSSIMAZIONE rp
                        # c[0],vp[0]= 0.01,2 per entrambi matlab e python
                        
                        # *** Python version ***
                        # print(-1j*sqrt(1-(0.01/2)**2))
                        # -0.999987499921874j                  #va fatto round a 4?
                        
                        # *** Matlab version ***
                        # disp(-i*sqrt(1-(0.01/2)^2))
                        # 0.0000 - 1.0000i
                        #________________________________________________

                    if c[l] > vs[m]:
                        rs=sqrt(((c[l]/vs[m])**2)-1)            ############# ok 
                        # rs=complex(np.round(rs.real,4),np.round(rs.imag,4))
                    else:                       
                        rs=-1j*(cmath.sqrt(1-(c[l]/vs[m])**2))  ############# ok 
                        # rs=complex(np.round(rs.real,4),np.round(rs.imag,4))

                    g=2*((vs[m]/c[l])**2)   # ok
                    P=k[h]*rp*d[m]          # ok
                    # P=np.round(P,8)         # ok
                    Q=k[h]*rs*d[m]          # ok
                    # Q=np.round(Q,8)         # ok   
                    cosenoP=(cmath.cos(P))                                          # ok
                    cosenoP=cosenoP.real                                            # ok
                    # cosenoP=np.round(cosenoP,8)                                     # ok
                    cosenoQ=(cmath.cos(Q))                                          # ok
                    cosenoQ=cosenoQ.real                                            # ok
                    # cosenoQ=np.round(cosenoQ,8)                                     # ok
                    senoP=cmath.sin(P)                                              # ok
                    # senoP=complex(np.round(senoP.real,9),np.round(senoP.imag,8))    # ok
                    senoQ=cmath.sin(Q)                                              # ok
                    # senoQ=complex(np.round(senoQ.real,8),np.round(senoQ.imag,8))    # ok
                    
                    # **fixed**
                    #______________________________________________________________________________________________
                    # Python version
                    #______________________________________________________________________________________________
                    # P=(-0.0007999899999374992j), Q=(-0.0007999599989999501j)
                    # sin(P)= -0-0.0007999900852676353j , sin(Q)= -0-0.0007999600843204864j   
                    # usando P arrotondato:
                    # sin(P)= -0-0.0007999900853301361j , sin(Q)= -0-0.0007999600853205366j
                    # g:[20000.]
                    # cosP=1.000000319992017 cosQ=1.000000319968017 [1.0000008]
                    
                    #_______________________________________________________________________________________________
                    # Matlab version
                    #_______________________________________________________________________________________________
                    # P=(0.0000e+00 - 7.9999e-04i), Q=(0.0000e+00 - 7.9996e-04i)
                    # g:[20000.]
                    # cosP=1.0000 cosQ=1.0000 [1.0000]
                    # sin(P)=0.0000e+00 - 7.9999e-04i, sin(Q)=0.0000e+00 - 7.9996e-04i
                    
                    
                    #risoluzione: P,Q,sinP,sinQ= round a 8
                    #risoluzione: cosP,cosQ round a 4
                    #______________________________________________________________________________________________
                    
                    A[0,0,m]=((g*cosenoP))-((g-1)*cosenoQ).real                          #ok   #rimane complex perchè la matrice allocata è complex  
                    # A[0,0,m]=complex(np.round(A[0,0,m].real,4),np.round(A[0,0,m].imag,4))     
                            
                    A[0,1,m]=1j*((((g-1)/rp)*senoP)+(g*rs*senoQ))                        #ok
                    # A[0,1,m]=complex(np.round(A[0,1,m].real,4),np.round(A[0,1,m].imag,4))
                    
                    A[0,2,m]=((-1/(rho[m]*(c[l]**2)))*(cosenoP-cosenoQ)).real   ############################## Non torna (=0)
                    # A[0,2,m]=complex(np.round(A[0,2,m].real,4),np.round(A[0,2,m].imag,4))
                    
                    A[0,3,m]=1j*(1/(rho[m]*(c[l]**2))*(((1/rp)*senoP)+(rs*senoQ)))       #ok
                    # A[0,3,m]=complex(np.round(A[0,3,m].real,4),np.round(A[0,3,m].imag,4))
                    
                    A[1,0,m]=-1j*((g*rp*senoP)+((g-1)*(1/rs)*senoQ))            ############################## Non torna (=4e-08j invece di= 4e-04j)
                    # A[1,0,m]=complex(np.round(A[1,0,m].real,4),np.round(A[1,0,m].imag,4)) 
                    A[1,1,m]=(-1*(g-1)*cosenoP)+(g*cosenoQ)                              #ok
                    # A[1,1,m]=complex(np.round(A[1,1,m].real,4),np.round(A[1,1,m].imag,4))
                    A[1,2,m]=1j*((1/((rho[m]*(c[l]**2))))*(rp*senoP+(1/rs)*senoQ))   ######################### Non torna (=0.0004j invece di= 2e-04j)
                    # A[1,2,m]=complex(np.round(A[1,2,m].real,4),np.round(A[1,2,m].imag,4))
                    A[1,3,m]=A[0,2,m]                                           ############################## Non torna di riflesso a A[0,2,m]
                    # A[1,3,m]=complex(np.round(A[1,3,m].real,4),np.round(A[1,3,m].imag,4))
                    
                    A[2,0,m]=(rho[m]*(c[l]**2)*g*(g-1)*(cosenoP-cosenoQ))
                    A[2,1,m]=1j*rho[m]*(c[l]**2)*(((g-1)**2)*(1/rp)*senoP+(g**2)*rs*senoQ)
                    A[2,2,m]=A[1,1,m]
                    A[2,3,m]=A[0,1,m]
                    A[3,0,m]=1j*rho[m]*(c[l]**2)*((g**2)*rp*senoP+((g-1)**2)*(1/rs)*senoQ)
                    A[3,1,m]=A[2,0,m]
                    A[3,2,m]=A[1,0,m]
                    A[3,3,m]=A[0,0,m]
                EE0=A[:,:,n-1]
                J=EE0[:,:]     
                for kk in range (1,n-1):
                    EE1=J[:,:,kk-1]
                    EE2=A[:,:,n-kk]
                    J[:,:,kk]=EE1*EE2 #EE0
                AA=J[:,:]#,n-1]
                K=(a1*AA[0,1]+b1*AA[1,1]+c1*AA[2,1]+d1*AA[3,1]).real #come devono essere: AA[0,1]= re, AA(2,1) =img, AA(3,1)=re , AA(4,1)=img
                # K=K.astype(float)
                L=a1*AA[0,0]+b1*AA[1,0]+c1*AA[2,0]+d1*AA[3,0]
                # L=L.astype(float)
                M=a2*AA[0,1]+b2*AA[1,1]+c2*AA[2,1]+d2*AA[3,1]
                # M=M.astype(float)
                N=(a2*AA[0,0]+b2*AA[1,0]+c2*AA[2,0]+d2*AA[3,0]).real
                # N=N.astype(float)
                HSK[l,h]=(K*N)-((L*M).real)
                
            f0=HSK[0,h]
            f1=HSK[1,h]
            cn=c[1]-(f1*(c[1]-c[0])/(f1-f0))
            c[0]=c[1]    
            c[1]=cn
            
            if (abs(c[0]-c[1])<eps):
                
                croot.append(cn)
                # croot[h]=cn
                kroot.append(k[h])
                # kroot[h]=k[h]
                T.append((2*pi*TH)/(k[h]*vels*cn))
                # T[h]=(2*pi*TH)/(k[h]*vels*cn)
                
                #problema, in quanto h supera la lunghezza di T ( se T è lungo 43, andrà da 0 a 42 e la posizione 43 non esisterà)
                #problema che si ripete ogni volta che non entro nell'if
                Tlastsampl= len(T)-1 #uso per risolvere il problema
                f=1/T[Tlastsampl] 
                break
        if (f>freqmax):
            break
    croot=np.array(croot)
    kroot=np.array(kroot)
    T=np.array(T)
    return(T,croot)