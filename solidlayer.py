import numpy as np 

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
    A=np.zeros((n,4,4),complex) 
######################################
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
                    rpn=sqrt(((c[l]/vp[n])**2)-1)                               
                else:
                    rpn=-1j*(cmath.sqrt(1-(c[l]/vp[n])**2))                     
                if c[l] > vs[n]:
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
                             
                for m in range (0,n):
                    if c[l] > vp[m]:
                        rp=sqrt(((c[l]/vp[m])**2)-1) 
                               
                    else:
                        rp=-1j*(cmath.sqrt(1-(c[l]/vp[m])**2)) 
                        
                        #________________________________________________

                    if c[l] > vs[m]:
                        rs=sqrt(((c[l]/vs[m])**2)-1)    
                               
                    else:                       
                        rs=-1j*(cmath.sqrt(1-(c[l]/vs[m])**2))  
                       
                        
                    g=2*((vs[m]/c[l])**2)   
 
                    P=k[h]*rp*d[m]          
                    Q=k[h]*rs*d[m]          
                    
                    cosenoP=(cmath.cos(P))   
                    cosenoP=cosenoP.real                                            
                    cosenoQ=(cmath.cos(Q))                                          
                    cosenoQ=cosenoQ.real                                           
                    senoP=cmath.sin(P)                                              
                    senoQ=cmath.sin(Q)                                              

                    #______________________________________________________________________________________________                
                    A[m,0,0]=((g*cosenoP))-((g-1)*cosenoQ)                #.real                           
                    A[m,0,1]=1j*((((g-1)/rp)*senoP)+(g*rs*senoQ)) #ok     #img                  
                    A[m,0,2]=((-1/(rho[m]*(c[l]**2)))*(cosenoP-cosenoQ))  #.real                 
                    A[m,0,3]=1j*(1/(rho[m]*(c[l]**2))*(((1/rp)*senoP)+(rs*senoQ))) #ok    #img
                                    
                    A[m,1,0]=-1j*((g*rp*senoP)+((g-1)*(1/rs)*senoQ)) #ok    #img      
                    A[m,1,1]=(-1*(g-1)*cosenoP)+(g*cosenoQ)                 #re
                    A[m,1,2]=1j*((1/((rho[m]*(c[l]**2))))*(rp*senoP+(1/rs)*senoQ))        #img
                    A[m,1,3]=A[m,0,2]                                       #re
                                                                                 
                    A[m,2,0]=(rho[m]*(c[l]**2)*g*(g-1)*(cosenoP-cosenoQ))   #re
                    A[m,2,1]=1j*rho[m]*(c[l]**2)*(((g-1)**2)*(1/rp)*senoP+(g**2)*rs*senoQ) #img
                    A[m,2,2]=A[m,1,1]                                       #re
                    A[m,2,3]=A[m,0,1]                                       #img
                    
                    A[m,3,0]=1j*rho[m]*(c[l]**2)*((g**2)*rp*senoP+((g-1)**2)*(1/rs)*senoQ) #img
                    A[m,3,1]=A[m,2,0]                                                      #re
                    A[m,3,2]=A[m,1,0]                                                      #img
                    A[m,3,3]=A[m,0,0]                                                      #re
                    
                EE0=A[n-1,:,:] #n=1 avrà solo 1 dimensione (0)
                
                                           
                J=EE0[:,:] # se n=1 non entro nel loop e J rimane così
                #_______________________________________________________________________ 
                for kk in range (1,n):# se n>=2 entro nel loop
                    if kk==1:
                        EE1=J[:,:]
                    else:
                        EE1=J[kk-1,:,:]

                    EE2=A[n-(kk+1),:,:]#+1 poichè se n=4, kk arriverà max a 3 nel ciclo
                    ET=np.dot(EE1,EE2)
                    if kk==1:          #reshape di J perchè solo 2D
                        J=np.reshape(J,(1,4,4))
                    J=np.insert(J,obj=len(J),values=ET,axis=0) 
                
                if n==1:
                    AA=J[:,:] #nel caso di n=1: 3d=0
                else:
                    AA=J[n-1,:,:] #nel caso di n>=2
                #________________________________________________________________________
                
                K=(a1*AA[0,1]+b1*AA[1,1]+c1*AA[2,1]+d1*AA[3,1]).real 
                L=a1*AA[0,0]+b1*AA[1,0]+c1*AA[2,0]+d1*AA[3,0]
                M=a2*AA[0,1]+b2*AA[1,1]+c2*AA[2,1]+d2*AA[3,1]
                N=(a2*AA[0,0]+b2*AA[1,0]+c2*AA[2,0]+d2*AA[3,0]).real
                HSK[l,h]=(K*N)-((L*M).real)
                
            f0=HSK[0,h]                       
            f1=HSK[1,h]                       
            ##################
            cn=c[1]-(f1*(c[1]-c[0])/(f1-f0))
            c[0]=c[1]    
            c[1]=cn
 
            if (abs(c[0]-c[1])<eps):
                
                croot.append(cn)
                kroot.append(k[h])
                T.append((2*pi*TH)/(k[h]*vels*cn)) 
                Tlastsampl= len(T)-1 #uso per risolvere il problema
                f=1/T[Tlastsampl] 
                break
       
        if (f>freqmax):
            break
    croot=np.array(croot)
    kroot=np.array(kroot)
    T=np.array(T)
    
    return(T,croot)
