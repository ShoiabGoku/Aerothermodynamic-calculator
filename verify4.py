# -*- coding: utf-8 -*-
# Nonequilibrium (finite-rate, Park) post-shock relaxation verification.
# 1-D steady conservation + Park R1-R5 single-temperature chemistry (same as
# the user's Fluent setup). Must start at the frozen shock state and relax to
# the equilibrium shock state already verified in verify.py.
import math
Ru=8.314462618;kB=1.380649e-23;hPl=6.62607015e-34;NA=6.02214076e23
Rk=8314.462618  # J/kmol/K
SP={'O2':dict(M=31.9988,mol=True,thv=2273.6,thr=2.079,sig=2,ge=3,hf=0.0),
    'N2':dict(M=28.0134,mol=True,thv=3393.5,thr=2.884,sig=2,ge=1,hf=0.0),
    'NO':dict(M=30.0061,mol=True,thv=2739.0,thr=2.452,sig=1,ge=4,hf=90290.0),
    'O':dict(M=15.9994,mol=False,ge=9,hf=249170.0),
    'N':dict(M=14.0067,mol=False,ge=4,hf=472680.0)}   # M in kg/kmol, hf J/mol
Tref=298.15
def hS(s,T):  # J/mol
    d=SP[s]
    if not d['mol']: return d['hf']+2.5*Ru*(T-Tref)
    vib=lambda th:d['thv']/(math.exp(d['thv']/th)-1)
    return d['hf']+3.5*Ru*(T-Tref)+Ru*(vib(T)-vib(Tref))
def cpS(s,T):
    d=SP[s]
    if not d['mol']: return 2.5*Ru
    x=d['thv']/T;ex=math.exp(x);return 3.5*Ru+Ru*x*x*ex/((ex-1)**2)
def qVol(s,T):
    d=SP[s];m=d['M']/1000/NA;qt=(2*math.pi*m*kB*T/(hPl*hPl))**1.5
    if not d['mol']: return qt*d['ge']
    return qt*(T/(d['sig']*d['thr']))*(1/(1-math.exp(-d['thv']/T)))*d['ge']
D_O2=2*SP['O']['hf'];D_N2=2*SP['N']['hf'];D_NO=SP['N']['hf']+SP['O']['hf']-SP['NO']['hf']
def Kp_O2(T):return qVol('O',T)**2/qVol('O2',T)*kB*T*math.exp(-D_O2/(NA*kB*T))
def Kp_N2(T):return qVol('N',T)**2/qVol('N2',T)*kB*T*math.exp(-D_N2/(NA*kB*T))
def Kp_NO(T):return qVol('N',T)*qVol('O',T)/qVol('NO',T)*kB*T*math.exp(-D_NO/(NA*kB*T))
SPL=['O2','N2','NO','O','N']
def mix_h(Y,T):  # J/kg
    return sum(Y[s]*hS(s,T)/(SP[s]['M']/1000)/1000*1 for s in SPL)  # hS J/mol / (kg/mol)
def mix_h2(Y,T): # cleaner: J/kg
    return sum(Y[s]*(hS(s,T)/(SP[s]['M']*1e-3)) for s in SPL)
def mix_Mbar(Y): # kg/kmol
    return 1.0/sum(Y[s]/SP[s]['M'] for s in SPL)
def mix_cp(Y,T):
    return sum(Y[s]*(cpS(s,T)/(SP[s]['M']*1e-3)) for s in SPL)
def solveT(Y,p,htarget):
    lo,hi=200.0,30000.0
    for _ in range(70):
        m=0.5*(lo+hi)
        if mix_h2(Y,m)>htarget: hi=m
        else: lo=m
    return 0.5*(lo+hi)
# ---- Park rates (SI, kmol, from thesis table) ----
# kf = A * T^b * exp(-E/(Rk*T)); [A]: m3/kmol/s (R1-R3 with third body), etc.
RXN=[dict(A=2.0e18,b=-1.5,E=4.95e8,typ='dissO2'),
     dict(A=7.0e18,b=-1.6,E=9.41e8,typ='dissN2'),
     dict(A=5.0e12,b=0.0,E=6.28e8,typ='dissNO'),
     dict(A=6.4e14,b=-1.0,E=3.19e8,typ='exchN2O'),   # N2+O -> NO+N
     dict(A=8.4e9 ,b=0.0,E=1.62e8,typ='exchNOO')]    # NO+O -> O2+N
def thirdBody(c,typ):
    if typ=='dissO2': eff={'O2':1,'N2':1,'NO':1,'O':5,'N':5}
    elif typ=='dissN2': eff={'O2':1,'N2':1,'NO':1,'O':4.3,'N':4.3}
    else: eff={'O2':1,'N2':1,'NO':22,'O':22,'N':22}
    return sum(eff[s]*c[s] for s in SPL)
def rates(Y,rho,T,mode):
    c={s:rho*Y[s]/SP[s]['M'] for s in SPL}      # kmol/m3
    w={s:0.0 for s in SPL}
    act={'s3':[0],'s4':[0,1],'s5':[0,1,2,3,4]}[mode]
    for i in act:
        r=RXN[i];kf=r['A']*T**r['b']*math.exp(-r['E']/(Rk*T))
        if r['typ']=='dissO2':
            Kc=Kp_O2(T)/(Rk*T);kb=kf/Kc;M3=thirdBody(c,r['typ'])
            q=kf*c['O2']*M3-kb*c['O']**2*M3
            w['O2']-=q;w['O']+=2*q
        elif r['typ']=='dissN2':
            Kc=Kp_N2(T)/(Rk*T);kb=kf/Kc;M3=thirdBody(c,r['typ'])
            q=kf*c['N2']*M3-kb*c['N']**2*M3
            w['N2']-=q;w['N']+=2*q
        elif r['typ']=='dissNO':
            Kc=Kp_NO(T)/(Rk*T);kb=kf/Kc;M3=thirdBody(c,r['typ'])
            q=kf*c['NO']*M3-kb*c['N']*c['O']*M3
            w['NO']-=q;w['N']+=q;w['O']+=q
        elif r['typ']=='exchN2O':
            Kc=Kp_N2(T)/Kp_NO(T);kb=kf/Kc
            q=kf*c['N2']*c['O']-kb*c['NO']*c['N']
            w['N2']-=q;w['O']-=q;w['NO']+=q;w['N']+=q
        else:  # NO+O -> O2+N
            Kc=Kp_NO(T)/Kp_O2(T);kb=kf/Kc
            q=kf*c['NO']*c['O']-kb*c['O2']*c['N']
            w['NO']-=q;w['O']-=q;w['O2']+=q;w['N']+=q
    return w  # kmol/m3/s
def stateFromY(Y,mflux,Pmom,Htot,rho_guess):
    rho=rho_guess
    for _ in range(60):
        u=mflux/rho;p=Pmom-mflux*u;h=Htot-0.5*u*u
        T=solveT(Y,p,h);Mb=mix_Mbar(Y);rho_new=p*Mb/(Rk*T)
        if abs(rho_new-rho)<1e-10*rho: rho=rho_new;break
        rho=0.5*(rho+rho_new)
    u=mflux/rho;p=Pmom-mflux*u;T=solveT(Y,p,Htot-0.5*u*u)
    return T,p,rho,u
def relax(T1,p1,V,mode,xmax=0.5):
    # frozen upstream air
    Mb_air=1/(0.233/SP['O2']['M']+0.767/SP['N2']['M'])
    Y={'O2':0.233,'N2':0.767,'NO':0.0,'O':0.0,'N':0.0}
    rho1=p1*mix_Mbar(Y)/(Rk*T1)
    h1=mix_h2(Y,T1)
    mflux=rho1*V;Pmom=p1+rho1*V*V;Htot=h1+0.5*V*V
    # frozen shock jump (no chemistry): iterate eps with fixed Y
    eps=0.1
    for _ in range(200):
        u2=eps*V;p2=p1+rho1*V*V*(1-eps);h2=h1+0.5*V*V*(1-eps*eps)
        T2=solveT(Y,p2,h2);rho2=p2*mix_Mbar(Y)/(Rk*T2)
        d=rho1/rho2-eps;eps+=0.6*d
        if abs(d)<1e-9:break
    Tf,pf,rhof,uf=T2,p2,rho2,eps*V
    # integrate dY/dx = w_i*M_i/(rho*u)
    x=0.0;T,p,rho,u=Tf,pf,rhof,uf
    hist=[(x,T,dict(Y))]
    dx=1e-7
    while x<xmax:
        w=rates(Y,rho,T,mode)
        dY={s:w[s]*SP[s]['M']/(rho*u) for s in SPL}
        mx=max(abs(dY[s]) for s in SPL)
        # adaptive dx: limit dY per step
        if mx>0: dx=min(max(2e-4/mx,1e-9),xmax/50, dx*2)
        else: dx=xmax/50
        Y2={s:max(0.0,Y[s]+dx*dY[s]) for s in SPL}
        ssum=sum(Y2.values());Y2={s:Y2[s]/ssum for s in SPL}
        T,p,rho,u=stateFromY(Y2,mflux,Pmom,Htot,rho)
        Y=Y2;x+=dx
        if len(hist)<4000: hist.append((x,T,dict(Y)))
        # stop when dT/dx tiny
        if x>1e-4 and abs(hist[-1][1]-hist[-2][1])/max(hist[-1][0]-hist[-2][0],1e-12)<1.0:
            break
    return dict(Tfrozen=Tf,pf=pf,hist=hist,Tend=T,Yend=Y,x=x)

# thesis condition
T1,p1,V=270.0,79.0,5000.0
for mode,Teq_expect in [('s3',7064),('s5',5360)]:
    r=relax(T1,p1,V,mode)
    print(f"{mode}: frozen T2={r['Tfrozen']:.0f} K -> relaxed T={r['Tend']:.0f} K at x={r['x']*1000:.2f} mm (equilibrium target ~{Teq_expect} K)")
    print(f"   Yend: "+", ".join(f"{s}={r['Yend'][s]:.4f}" for s in SPL))
    # relaxation length: 95% of temperature drop
    Tf,Te=r['Tfrozen'],r['Tend']
    x95=next((x for x,T,_ in r['hist'] if (Tf-T)>=0.95*(Tf-Te)),None)
    print(f"   95% relaxation length = {x95*1000:.3f} mm, points={len(r['hist'])}")
