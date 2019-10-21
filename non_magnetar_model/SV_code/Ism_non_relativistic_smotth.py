import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize




t=np.logspace(0,3.,400)





#NU=np.logspace(-2.,8.,100)



nuc_init=3.5e14
nua_init=1e6
M_sun=2e33
num_init=1e8 #GHz
fm_init=28.5e-3

beta0=0.5
c=3e10 ## speed of light
s=10.
p=2.2
y=num_init/nua_init
mp=1.67e-24
me=9.1e-28
epse=0.3
n=0.02
sm=2.
sc=2.
st=5.

Mej=0.006*M_sun
T_init=25.
####-------------------------------------------------------------------------------------------------------------------########



def beta_fun(T):
 beta=beta0*(1.+(T/T_init)**(3.*s/5.))**(-1./s)
 return beta

def r_fun(T):
 return beta_fun(T)*c*T




def nuc(T):
 vc=nuc_init*(T/T_init)**-2.*(beta_fun(T)/beta_fun(T_init))**-3.
 return np.log10(vc)


def fun_root2(T):
 return nua_init*(r_fun(T)/r_fun(T_init))**(3./5.)*(beta_fun(T)/beta_fun(T_init))**(-8./5.)-num_init*(beta_fun(T)/beta_fun(T_init))**5.
 
  
teq=optimize.brentq(fun_root2,1e-2,5000.)    ##### time where nua = num for nua>num case

'''  
plt.loglog(t,fun_root(t))  ####just to check the function
plt.show()
'''


#print teq

def num(T):
 vm=num_init*(beta_fun(T)/beta_fun(T_init))**5.
 return np.log10(vm)




def nua(T):
 if T<=teq:
  va=nua_init*(r_fun(T)/r_fun(T_init))**(3./5.)*(beta_fun(T)/beta_fun(T_init))**(-8./5.)
  return np.log10(va)
 elif T>teq:
  nua_init1=nua_init*(r_fun(teq)/r_fun(T_init))**(3./5.)*(beta_fun(teq)/beta_fun(T_init))**(-8./5.)
  va=nua_init1*(r_fun(T)/r_fun(teq))**(2./(p+4.))*(beta_fun(T)/beta_fun(teq))**((5.*p-2.)/(p+4.))
  return np.log10(va)



def fm(T):
 return np.log10(fm_init*(beta_fun(T)/beta_fun(T_init))*(r_fun(T)/r_fun(T_init))**3.)



def flux(lgnu,T):
 beta1 = 1./3.

 if num(T) >= nua(T):
  tp = 10**(-5.*(num(T)-nua(T))/3.0)
 elif num(T) < nua(T):
  tp = 10**(-(p+4.)*0.5*(num(T)-nua(T)))

 taunu = tp*10**((lgnu-num(T))*(-5./3.))*(1.+10.**(sm*((p/2.)+(1./3.))*(lgnu-num(T))))**(-1./sm)*(1.+10.**(sc*0.5*(lgnu-nuc(T)))) **(-1./sc)
 
	
#================= calculating snu =================

 if num(T) >= nua(T):
  smt = fm(T)+(5.0*(num(T)-nua(T))/3.0)  #S_nu(nu = nu_p)
 elif num(T) < nua(T):
  smt = fm(T)+(p+4.)*0.5*(num(T)-nua(T))
 
 #snu = smt+(log10(10**(2.0*sm*(lgnu-numt))+ &
 #              10**(5.0*sm*(lgnu-numt)*0.5)))/sm

 a1 = 10**(2.0*sm*(lgnu-num(T)))
 a2= 10**(5.0*sm*(lgnu-num(T))*0.5)
 snu = smt+(np.log10(a1+a2))/sm

 #fnu = snu
 if taunu <= 5e-8:
  lgfl = snu+np.log10(taunu)
 else:
  lgfl = snu+np.log10(1.-np.exp(-taunu))
 
 return lgfl


'''

f=np.zeros(np.size(t))
f1=np.zeros(np.size(t))
f2=np.zeros(np.size(t))
f3=np.zeros(np.size(t))
f4=np.zeros(np.size(t))
#f5=np.zeros(np.size(NU))




for i in range(np.size(t)):
 
 f[i]=fm(t[i])
 f1[i]=10**flux(np.log10(350e9),t[i])
 f2[i]=nua(t[i])
 f3[i]=num(t[i])
 f4[i]=nuc(t[i])




dta=np.loadtxt('350GHz.mdl')
fl1=10**(dta[:,5])
fl=np.zeros(226)
for i in range(226):
 fl[i]=fl1[i]+f1[i]


d=np.loadtxt('350GHz')
td=np.log10(d[:,0])
Fl=np.log10(d[:,1]*1e-3)
Z=d[:,2]*1e-3
print Z
err=np.zeros(len(td))
for i in range(len(td)):
 err[i]=0.434*Z[i]/10**Fl[i]
'''

'''

d2=np.loadtxt('1.4GHz')
td2=np.log10(d2[1:,0])
fl2=np.log10(d2[1:,1]*1e-3)
Z2=0.434*d2[1:,2]*1e-3
err2=np.zeros(len(td2))
for i in range(len(td2)):
 err2[i]=Z2[i]/10**fl2[i]

plt.errorbar(td2,fl2,err2,fmt='go')


d=np.loadtxt('350GHz')
td=np.log10(d[:,0])
Fl=np.log10(d[:,1]*1e-3)
Z=0.434*d[:,2]*1e-3
err=np.zeros(len(td))
for i in range(len(td)):
 err[i]=Z[i]/10**Fl[i]
'''



# '#176F5A' 45 GHz, #0F8AF1 0.6, 6 5BDD1A, 15 4D4BF6
'''
plt.plot(t,fl1,color='#E74C37',linestyle='--',label='jet',linewidth=3.)
plt.errorbar(10**td,10**Fl,Z,fmt='o',color='#EA4D0D',markersize=10.0)
plt.plot(t,f1,color='#626567',linestyle='--',label='cocoon',linewidth=3.)
plt.plot(t,fl,color='#EA4D0D',linestyle='-',label='combined',linewidth=5.)
#plt.ylim(-8,-1)
#plt.title('90GHz',size=20.)
plt.ylabel('Flux (Jy)',size=40.0)
plt.xlabel('Time (days)',size=40.0)
plt.yscale('log')
plt.xscale('log')
plt.tick_params(axis='both',which='major',length=20.0,labelsize=25.0)
plt.tick_params(axis='both',which='minor',length=14.0,labelsize=20.0)
plt.legend(loc='best',prop={'size': 25})
plt.show()


plt.plot(np.log10(t),f3)
plt.plot(np.log10(t),f2)
#plt.plot(np.log10(t),f)
#plt.loglog(t,f4)
#plt.errorbar(td,Fl,err,fmt='ro')
plt.show()
'''






'''
plt.ylabel('Flux (Jy)',size=30.0)
plt.xlabel('Time (days)',size=30.0)
plt.tick_params(axis='both',which='major',length=18.0,labelsize=20.0)
plt.tick_params(axis='both',which='minor',length=12.0,labelsize=15.0)
plt.legend(loc='best',prop={'size': 20})

'''






