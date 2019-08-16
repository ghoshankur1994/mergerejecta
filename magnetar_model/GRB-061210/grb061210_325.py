import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
fontP=FontProperties()
 
#---------------------------------------


#plt.subplot(2,1,1)
def Fong16(E52,mej,n0,dl26,nuobs):
    dl27 = dl26/10.
    t0 = (300./365.)*(E52**-0.5)*((mej/0.01)**(5./6.))/(n0**(1./3.))
    p=2.2
    a=(5.*p-3)*0.25;b=0.25*(5.*p-7.); d=0.25*(p+1.); g=0.5*(p-1)
    c1 = 3.e5*( 1.1**(0.5*(5.*p-7.)) )*( 4.3**(-0.5*(p-1)))
    #print 'c1 is', c1
    f0 = c1*(E52**a)*((mej/0.01)**-b)*(n0**d)/(dl27*dl27*((nuobs/6.)**g))
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
plt.subplot(1,2,1) 
#GRB100625a
dl_26=2200.1e6*3.08e18/1.e26;  epoch=7.5
 
mej1=0.1
p = 2.2
y = -(15.*p-21.)/10.
 
#plt.text(20.,50.,'GRB160821b, $M_{ej} = 0.01$',horizontalalignment='center',verticalalignment='center')
 
plt.text(200.,1.,' $M_{ej} = 0.1$',horizontalalignment='center',verticalalignment='center')
 
#td1= MetzBower(10.,1.,0.2,dl_26)[0]
#fp1= MetzBower(10.,1.,0.2,dl_26)[1]
#print td1,fp1
#tin=np.logspace(np.log10(1.5),np.log10(td1),20)
#plt.loglog(tin,10**fout,'g-',label='E52=10, n0 = 1.')
 
 
 
#----- 325 upper limit ----
plt.scatter(9.9,0.36,marker='v',color='k',s=50)



td1= Fong16(1.,mej1,0.001,dl_26,0.325)[0]
print (td1)
tin=np.logspace(np.log10(1.),1.+np.log10(td1),20)
 
fp1= Fong16(1.,mej1,0.001,dl_26,0.325)[1]
fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
fout=np.where(tin < td1,fout1,fout2)
plt.loglog(tin,10**(fout-3),'r--',lw=2.0,label='E52=1.0, n0 = 0.001')



#Fong16(E52,mej,n0,dl26,nuobs)
td1= Fong16(1.,mej1,0.01,dl_26,0.325)[0]
print(td1)
tin=np.logspace(np.log10(1.),1.+np.log10(td1),20)
#
fp1= Fong16(1.,mej1,0.01,dl_26,0.6)[1]
fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
fout=np.where(tin < td1,fout1,fout2)
plt.loglog(tin,10**(fout-3),'g--',lw=2.0,label='E52=1.0, n0 = 0.01')


#Fong16(E52,mej,n0,dl26,nuobs)
td1= Fong16(1.,mej1,0.1,dl_26,0.325)[0]
tin=np.logspace(np.log10(1.),1.+np.log10(td1),20)
#
fp1= Fong16(1.,mej1,0.1,dl_26,0.325)[1]
fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
fout=np.where(tin < td1,fout1,fout2)
plt.loglog(tin,10**(fout-3),'m--',lw=2.0,label='E52=1.0, n0 = 0.1')


#Fong16(E52,mej,n0,dl26,nuobs)
td1= Fong16(1.,mej1,1.,dl_26,0.325)[0]
tin=np.logspace(np.log10(1.),1.+np.log10(td1),20)
#
fp1= Fong16(1.,mej1,1.,dl_26,0.325)[1]
fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
fout=np.where(tin < td1,fout1,fout2)
plt.loglog(tin,10**(fout-3),'b--',lw=2.0,label='E52=1.0, n0 = 1.0')
plt.legend(prop=fontP,ncol=1,loc='lower right').draw_frame(False)
plt.xlabel(r'($t-t_0$)yrs',size=12)
plt.ylabel('Flux (mJy)',size=12)
plt.subplot(1,2,2) 
#GRB100625a
dl_26=2200.1e6*3.08e18/1.e26;  epoch=7.5
 
mej1=0.05
p = 2.2
y = -(15.*p-21.)/10.

plt.text(40.,5.,'$M_{ej} = 0.05$',horizontalalignment='center',verticalalignment='center')

plt.scatter(11.,0.36,marker='v',color='k',s=50)



td1= Fong16(1.,mej1,0.001,dl_26,0.325)[0]
print (td1)
tin=np.logspace(np.log10(1.),1.+np.log10(td1),20)
 
fp1= Fong16(1.,mej1,0.001,dl_26,0.325)[1]
fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
fout=np.where(tin < td1,fout1,fout2)
plt.loglog(tin,10**(fout-3),'r--',lw=2.0,label='E52=1.0, n0 = 0.001')
plt.xlabel(r'($t-t_0$)yrs',size=12)
plt.ylabel('Flux (mJy)',size=12)
#plt.legend(prop=fontP,ncol=1,loc='lower right').draw_frame(False)


#Fong16(E52,mej,n0,dl26,nuobs)
td1= Fong16(1.,mej1,0.01,dl_26,0.325)[0]
tin=np.logspace(np.log10(1.),1.+np.log10(td1),20)
#
fp1= Fong16(1.,mej1,0.01,dl_26,0.325)[1]
fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
fout=np.where(tin < td1,fout1,fout2)
plt.loglog(tin,10**(fout-3),'g--',lw=2.0,label='E52=1.0, n0 = 0.01')


#Fong16(E52,mej,n0,dl26,nuobs)
td1= Fong16(1.,mej1,0.1,dl_26,0.325)[0]
tin=np.logspace(np.log10(1.),1.+np.log10(td1),20)
#
fp1= Fong16(1.,mej1,0.1,dl_26,0.325)[1]
fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
fout=np.where(tin < td1,fout1,fout2)
plt.loglog(tin,10**(fout-3),'m--',lw=2.0,label='E52=1.0, n0 = 0.1')


#Fong16(E52,mej,n0,dl26,nuobs)
td1= Fong16(1.,mej1,1.,dl_26,0.325)[0]
tin=np.logspace(np.log10(1.),1.+np.log10(td1),20)
#
fp1= Fong16(1.,mej1,1.,dl_26,0.325)[1]
fout1= np.log10(fp1)+3.*(np.log10(tin)-np.log10(td1))
fout2= np.log10(fp1)+y*(np.log10(tin)-np.log10(td1))
fout=np.where(tin < td1,fout1,fout2)
plt.loglog(tin,10**(fout-3),'b--',lw=2.0,label='E52=1.0, n0 = 1.0')

#plt.subplot(2,1,1)
#plt.axhline(y=0.1,xmin=0,xmax=1,ls='-.',c='k')
#plt.axvline(x=epoch,ymin=0,ymax=1,ls='-.',c='k')
fontP.set_size('x-small')

#pylab.legend(prop=fontP,ncol=3,loc='lower left')
plt.xlabel(r'($t-t_0$)yrs',size=12)
plt.ylabel('Flux (mJy)',size=12)
plt.gca().set_xlim(1.,100.)
plt.gca().set_ylim(1.e-4,1.e1)
#plt.grid()
plt.subplots_adjust(top=0.8, bottom=0.15, left=0.15, right=0.95, hspace=0.5,
                    wspace=0.5)
plt.savefig('grb061210_325.jpg',dpi=1000)
plt.show()



























