#This is to compare the equations of three different papers
#assuming magnetar energized merger ejecta model
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import subprocess
fontP=FontProperties()
 
#---------------------------------------
def MetzBower(E52,n0,b0,dl26):
    t0=1.2*((E52)**(1./3.))/( (n0**(1./3.))*(b0**(5./3.)) )
    dl28 = dl26/100.
    f0=3.*(E52)*(n0**0.83)*(b0**2.3)/(dl28*dl28)
    #This is done for 1.4gHz
    tf=(t0,f0)
    return tf

 
def Horesh16(E52,n0,nuobs,dl26):
    #Optically thin model
    #Assumed epse = epsb = 0.1
    p=2.2
    dl28 = dl26/100.
    a = (p+2)/(3.*p+2)
    b=4./(3.*p+2)
    d =2.*(p+4)/(3.*p+2.)
 
    f = (2.*p+3.)/(3.*p+2)
    g = (8-3.*p)/(6.*p+4.)
    h = (5*p-5.)/(3.*p+2)
 
    t0 = (170./365.)*((E52/3.)**(a))*(n0**b)/((nuobs/3.)**d)  #in yrs
    f0 = 5.*(E52**f)*(n0**g)*((nuobs/3.)**h)/(dl28*dl28)
    tf=(t0,f0)
    return tf
 
#def Horesh16_1(E52,n0,nuobs):
#    #Optically thin model
#    #Assumed epse = epsb = 0.1
#    #Here, t0 is the point where obs frequency crosses num. This is useless because num is almost always below nua.
#    t0 = (120./365.)*((E52/3.)**(1./3.))/((nuobs/3.)**(2./3.))  #in yrs
#
#    f0 = 1.
#    tf=(t0,f0)
#    return tf
 
def Fong16(E52,mej,n0,dl26,nuobs):
    dl27 = dl26/10.
    t0 = (300./365.)*(E52**-0.5)*((mej/0.01)**(5./6.))/(n0**(1./3.))
    p=2.2
    a=(5.*p-3)*0.25;b=0.25*(5.*p-7.); d=0.25*(p+1.); g=0.5*(p-1)
    c1 = 3.e5*( 1.1**(0.5*(5.*p-7.)) )*( 4.3**(-0.5*(p-1)))
#    print("Bar=\t",E52)
    #print 'c1 is', c1
    f0 = c1*(E52**a)*((mej/0.01)**-b)*(n0**d)/(dl27*dl27*((nuobs/6.)**g))
 #   print ("Flux=\t",f0)
    tf=(t0,f0)
    return tf
#----------------------------Main code-----------------------
#Give nuobs in GHz
 
#GRB160821b
#dl_26 = 771.e6*3.08e18/1.e26; epoch=1.5
 
#GRB140903a
#dl_26=1874e6*3.08e18/1.e26;  epoch=3.5
 
#GRB100625a
#dl_26=2536.6e6*3.08e18/1.e26;  epoch=7.5
 
#GRB100625a
dl_26=771.1e6*3.08e18/1.e26;  epoch=7.5
 
mej1=0.1
p = 2.2
y = -(15.*p-21.)/10.
 
#plt.text(20.,50.,'GRB160821b, $M_{ej} = 0.01$',horizontalalignment='center',verticalalignment='center')
 
#plt.text(20.,25.,'GRB050709, $M_{ej} = 0.1$',horizontalalignment='center',verticalalignment='center')
 
#td1= MetzBower(10.,1.,0.2,dl_26)[0]
#fp1= MetzBower(10.,1.,0.2,dl_26)[1]
#print td1,fp1
#tin=np.logspace(np.log10(1.5),np.log10(td1),20)
#plt.loglog(tin,10**fout,'g-',label='E52=10, n0 = 1.')
 
 
 
#----- 610 upper limit ----
#plt.scatter(11.,0.12,marker='v',color='k',s=100)
 
#
#
ag47=open('datafile.txt','w')
ag47.write('#Energy\t n_density\t time\t flux\n')
#ism=np.array([0.001,0.01,0.1,1.0])
scale=0.001
number_density=0.001
#Fong16(E52,mej,n0,dl26,nuobs)

#for j in range(0,len(ism)):      
while scale <= 1.0 : 
      number_density,counter=0.001,0
#      print ("Scale=\t",(scale))
      while number_density <= 1.0 :
             td1= Fong16(scale,mej1,number_density,dl_26,0.6)[0]
             print ("Time=\t",td1)
             tin=np.logspace(np.log10(1.),1.+np.log10(td1),20)
#      print(tin) 
             fp1= Fong16(scale,mej1,number_density,dl_26,0.6)[1]
             fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
             fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
             fout=np.where(tin < td1,fout1,fout2)
#      print(fout)
             print ('\nscale=\t{:5.5f}'.format(scale),'number_density=\t{:5.5f}'.format(number_density),'tin=\t{:5.5f}'.format(tin[counter]))
#             for sh in range(0,len(fout)):
                 
             ag47.write('%5.5f\t%5.5f\t%5.5f\t\t%3.5e\n' % (scale,number_density,tin[counter],10**(fout[counter]-3)))
             
             plt.loglog(tin,10**(fout-3),'r--',lw=2.0,label='E52=0.2, n0 = 0.1')
 #            print ("Number_density=\t",number_density)
             number_density,counter=number_density*10,counter+1
      ag47.write('#-------------------------------------------------\n')
      scale=scale+0.01
        
ag47.close()      
#plt.loglog(tin,10**(fout-3),'r--',lw=3.0,label='E52=0.5, n0 = 0.1')
#
 
#subprocess.call("gedit datafile.txt &",shell=True)
#fp1= Fong16(5.e-1,mej1,1.e-2,dl_26,1.4)[1]
#fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
#fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
#fout=np.where(tin < td1,fout1,fout2)
#plt.loglog(tin,10**(fout-3),'r-',lw=3.0,label='E52=0.5, n0 = 1.e-2, 1.4GHz')
 
 
#
#
 
#Horesh16(E52,n0,nuobs,dl26)
#td1 = Horesh16(1.,1.e-4,0.6,dl_26)[0]
#f1 = Horesh16(1.,1.e-4,0.6,dl_26)[1]
#print td1,fp1
#tin=np.logspace(np.log10(1.),1.+np.log10(td1),10)
#fout = np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
#fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
#fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
#fout=np.where(tin < td1,fout1,fout2)
#plt.loglog(tin,10**fout,'r--',label='E52=1, n0 = 1.e-4, Horesh')
 
#
#-------------------------------changing number density------------
#
#td1= MetzBower(1.,0.01,0.2,dl_26)[0]
#fp1= MetzBower(1.,0.01,0.2,dl_26)[1]
#print td1,fp1
#tin=np.logspace(np.log10(1.5),np.log10(td1),10)
#fout = np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
#plt.loglog(tin,10**fout,'r-',label='E52=10, n0 =0.01')
 
#
#
#Fong16(E52,mej,n0,dl26,nuobs)
#td1= Fong16(2.e-1,mej1,1.0,dl_26,0.6)[0]
#tin=np.logspace(np.log10(1.),1.+np.log10(td1),20)
#
#fp1= Fong16(2.e-1,mej1,1.0,dl_26,0.6)[1]
#fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
#fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
#fout=np.where(tin < td1,fout1,fout2)
#plt.loglog(tin,10**(fout-3),'g--',lw=4.0,label='E52=0.2, n0 = 1.0')
#plt.loglog(tin,10**(fout-3),'g--',lw=4.0,label='E52=0.5, n0 = 0.1, 610MHz')
#
#fp1= Fong16(5.e-1,mej1,0.1,dl_26,1.4)[1]
#fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
#fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
#fout=np.where(tin < td1,fout1,fout2)
#plt.loglog(tin,10**(fout-3),'g-',lw=4.0,label='E52=0.5, n0 = 0.1, 1.4GHz')
 
#------------------------------change E52---------------
 
#Fong16(E52,mej,n0,dl26,nuobs)
#td1= Fong16(1.,mej1,0.1,dl_26,0.6)[0]
#tin=np.logspace(np.log10(1.),1.+np.log10(td1),20)
#
#fp1= Fong16(1.,mej1,0.1,dl_26,0.6)[1]
#fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
#fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
#fout=np.where(tin < td1,fout1,fout2)
#plt.loglog(tin,10**(fout-3),'m--',lw=4.0,label='E52=1.0, n0 = 0.1')
#
#fp1= Fong16(1.,mej1,0.01,dl_26,1.4)[1]
#fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
#fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
#fout=np.where(tin < td1,fout1,fout2)
#plt.loglog(tin,10**(fout-3),'m-',lw=4.0,label='E52=1., n0 = 1.e-2, 1.4GHz')
 
#------------------------------change E52---------------
 
#Fong16(E52,mej,n0,dl26,nuobs)
#td1= Fong16(1.,mej1,1,dl_26,0.6)[0]
#tin=np.logspace(np.log10(1.),1.+np.log10(td1),20)
#
#fp1= Fong16(1.,mej1,1.0,dl_26,0.6)[1]
#fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
#fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
#fout=np.where(tin < td1,fout1,fout2)
#plt.loglog(tin,10**(fout-3),'b--',lw=4.0,label='E52=1.0, n0 = 1.0')
#
#fp1= Fong16(1.,mej1,0.01,dl_26,1.4)[1]
#fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
#fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
#fout=np.where(tin < td1,fout1,fout2)
#plt.loglog(tin,10**(fout-3),'m-',lw=4.0,label='E52=1., n0 = 1.e-2, 1.4GHz')
 
 
 
 
#Horesh16(E52,n0,nuobs,dl26)
#td1 = Horesh16(10.,0.01,0.6,dl_26)[0]
#f1 = Horesh16(10.,0.01,0.6,dl_26)[1]
#print td1,fp1
#tin=np.logspace(np.log10(1.),1.+np.log10(td1),10)
#fout = np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
#fout1 = np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
#fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
#fout=np.where(tin < td1,fout1,fout2)
#plt.loglog(tin,10**fout,'r-',label='E52=10, n0 = 0.01, Horesh')
 
plt.axhline(y=0.1,xmin=0,xmax=1,ls='-.',c='k')
plt.axvline(x=epoch,ymin=0,ymax=1,ls='-.',c='k')
fontP.set_size('x-small')
#plt.legend(prop=fontP,ncol=3,loc='best').draw_frame(False)
#pylab.legend(prop=fontP,ncol=3,loc='lower left')
plt.xlabel(r'($t-t_0$)yrs',size=12)
plt.ylabel('Flux (mJy)',size=12)
plt.gca().set_xlim(1.,50.)
plt.gca().set_ylim(1.e-3,1.e2)
#plt.grid()
plt.savefig('grb100625a_610.eps')
plt.show()
 
 
