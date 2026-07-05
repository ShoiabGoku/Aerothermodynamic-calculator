# -*- coding: utf-8 -*-
# Two-temperature (Park T-Tv) post-shock relaxation verification.
# Checks: (1) frozen 2-T jump equals the calorically-perfect (gamma=1.4) jump,
# (2) Tv rises from T1 and merges with T, (3) end state = 1-T equilibrium.
import math
Ru=8.314462618;kB=1.380649e-23;hPl=6.62607015e-34;NA=6.02214076e23
SP={'O2':dict(M=0.0319988,mol=True ,thv=2273.6,thr=2.079,sig=2,ge=3,hf=0.0),
    'N2':dict(M=0.0280134,mol=True ,thv=3393.5,thr=2.884,sig=2,ge=1,hf=0.0),
    'NO':dict(M=0.0300061,mol=True ,thv=2739.0,thr=2.452,sig=1,ge=4,hf=90290.0),
    'O' :dict(M=0.0159994,mol=False,ge=9,hf=249170.0),
    'N' :dict(M=0.0140067,mol=False,ge=4,hf=472680.0)}  # M kg/mol, hf J/mol
SPL=['O2','N2','NO','O','N'];Tref=298.15
def qVol(s,T):
    d=SP[s];m=d['M']/NA;qt=(2*math.pi*m*kB*T/(hPl*hPl))**1.5
    if not d['mol']:return qt*d['ge']
    return qt*(T/(d['sig']*d['thr']))*(1/(1-math.exp(-d['thv']/T)))*d['ge']
D_O2=2*SP['O']['hf'];D_N2=2*SP['N']['hf'];D_NO=SP['N']['hf']+SP['O']['hf']-SP['NO']['hf']
def Kp_O2(T):return qVol('O',T)**2/qVol('O2',T)*kB*T*math.exp(-D_O2/(NA*kB*T))
def Kp_N2(T):return qVol('N',T)**2/qVol('N2',T)*kB*T*math.exp(-D_N2/(NA*kB*T))
def Kp_NO(T):return qVol('N',T)*qVol('O',T)/qVol('NO',T)*kB*T*math.exp(-D_NO/(NA*kB*T))
# ---- energies (per kg of species) ----
def htr(s,T):   # trans+rot enthalpy + formation, J/kg (NO vibration)
    d=SP[s];cp=(3.5 if d['mol'] else 2.5)*Ru/d['M']
    return d['hf']/d['M']+cp*(T-Tref)
def ev(s,Tv):   # vibrational energy J/kg (ref: ev(Tref) subtracted for consistency)
    d=SP[s]
    if not d['mol']:return 0.0
    f=lambda t:Ru/d['M']*d['thv']/(math.exp(d['thv']/t)-1)
    return f(Tv)-f(Tref)
def cvv(s,Tv):
    d=SP[s]
    if not d['mol']:return 0.0
    x=d['thv']/Tv;e=math.exp(x)
    return Ru/d['M']*x*x*e/((e-1)**2)
def mixHtr(Y,T):return sum(Y[s]*htr(s,T) for s in SPL)
def mixEv(Y,Tv):return sum(Y[s]*ev(s,Tv) for s in SPL)
def mixMbar(Y):return 1.0/sum(Y[s]/SP[s]['M'] for s in SPL)
def solveT_tr(Y,target):  # T from trans-rot enthalpy
    lo,hi=200.0,40000.0
    for _ in range(64):
        m=0.5*(lo+hi)
        if mixHtr(Y,m)>target:hi=m
        else:lo=m
    return 0.5*(lo+hi)
# ---- Millikan-White + Park relaxation time (mixture) ----
def tauVib(s,Y,rho,p,T):
    d=SP[s];Mbar=mixMbar(Y)
    n_tot=rho/Mbar*NA               # number density 1/m3
    p_atm=p/101325.0
    inv=0.0;xs={t:Y[t]/SP[t]['M']*Mbar for t in SPL}  # mole fractions
    for t in SPL:
        if xs[t]<=0:continue
        mu=SP[s]['M']*SP[t]['M']/(SP[s]['M']+SP[t]['M'])*1000.0  # g/mol reduced
        A=1.16e-3*math.sqrt(mu)*d['thv']**(4.0/3.0)
        ptau=math.exp(A*(T**(-1.0/3.0)-0.015*mu**0.25)-18.42)    # atm*s
        inv+=xs[t]/(ptau/p_atm)
    tMW=1.0/inv if inv>0 else 1e6
    cbar=math.sqrt(8*Ru*T/(math.pi*SP[s]['M']))
    sig=1e-21*(50000.0/T)**2
    tP=1.0/(sig*cbar*n_tot)
    return tMW+tP
# ---- Park rates, dissociation at Ta=sqrt(T*Tv) ----
PARK=[dict(A=2.0e15,b=-1.5,E=4.95e5,t='dO2'),dict(A=7.0e15,b=-1.6,E=9.41e5,t='dN2'),
      dict(A=5.0e9 ,b=0.0 ,E=6.28e5,t='dNO'),dict(A=6.4e11,b=-1.0,E=3.19e5,t='xN2O'),
      dict(A=8.4e6 ,b=0.0 ,E=1.62e5,t='xNOO')]   # A m3/mol/s, E J/mol
def rates(Y,rho,T,Tv,mode,twoT):
    c={s:rho*Y[s]/SP[s]['M'] for s in SPL}
    w={s:0.0 for s in SPL}
    act={'s3':[0],'s4':[0,1],'s5':[0,1,2,3,4]}[mode]
    Ta=math.sqrt(T*max(Tv,50.0)) if twoT else T
    M3=lambda eff:sum(eff.get(s,1)*c[s] for s in SPL)
    for i in act:
        r=PARK[i]
        if r['t'] in('dO2','dN2','dNO'):
            kf=r['A']*Ta**r['b']*math.exp(-r['E']/(Ru*Ta))
            kfT=r['A']*T**r['b']*math.exp(-r['E']/(Ru*T))
            if r['t']=='dO2':
                kb=kfT/(Kp_O2(T)/(Ru*T));m3=M3({'O':5,'N':5})
                q=(kf*c['O2']-kb*c['O']**2)*m3;w['O2']-=q;w['O']+=2*q
            elif r['t']=='dN2':
                kb=kfT/(Kp_N2(T)/(Ru*T));m3=M3({'O':4.3,'N':4.3})
                q=(kf*c['N2']-kb*c['N']**2)*m3;w['N2']-=q;w['N']+=2*q
            else:
                kb=kfT/(Kp_NO(T)/(Ru*T));m3=M3({'NO':22,'O':22,'N':22})
                q=(kf*c['NO']-kb*c['N']*c['O'])*m3;w['NO']-=q;w['N']+=q;w['O']+=q
        else:
            kf=r['A']*T**r['b']*math.exp(-r['E']/(Ru*T))
            if r['t']=='xN2O':
                kb=kf/(Kp_N2(T)/Kp_NO(T))
                q=kf*c['N2']*c['O']-kb*c['NO']*c['N'];w['N2']-=q;w['O']-=q;w['NO']+=q;w['N']+=q
            else:
                kb=kf/(Kp_NO(T)/Kp_O2(T))
                q=kf*c['NO']*c['O']-kb*c['O2']*c['N'];w['NO']-=q;w['O']-=q;w['O2']+=q;w['N']+=q
    return w
def relax2T(T1,p1,V,mode,xmax=1.0):
    Y={'O2':0.233,'N2':0.767,'NO':0.0,'O':0.0,'N':0.0}
    Mb=mixMbar(Y);rho1=p1*Mb/(Ru*T1)
    mflux=rho1*V;Pmom=p1+rho1*V*V
    Htr_tot=mixHtr(Y,T1)+0.5*V*V           # trans-rot + KE conserved separately? no:
    Ev1=mixEv(Y,T1)                         # total h = htr+Ev+KE conserved
    Htot=mixHtr(Y,T1)+Ev1+0.5*V*V
    # frozen 2-T jump: vibration frozen at T1 (Ev const), chemistry frozen
    eps=0.1
    for _ in range(200):
        u2=eps*V;p2=p1+rho1*V*V*(1-eps)
        htr2=Htot-Ev1-0.5*u2*u2
        T2=solveT_tr(Y,htr2);rho2=p2*Mb/(Ru*T2)
        d=rho1/rho2-eps;eps+=0.6*d
        if abs(d)<1e-9:break
    T,p,rho,u=T2,p2,rho2,eps*V;Tv=T1;x=0.0;dx=1e-8
    hist=[(x,T,Tv,dict(Y))]
    def state(Yc,Evc,rg):
        r=rg
        for _ in range(50):
            uu=mflux/r;pp=Pmom-mflux*uu
            TT=solveT_tr(Yc,Htot-Evc-0.5*uu*uu)
            rn=pp*mixMbar(Yc)/(Ru*TT)
            if abs(rn-r)<1e-9*r:r=rn;break
            r=0.5*(r+rn)
        uu=mflux/r;pp=Pmom-mflux*uu
        return solveT_tr(Yc,Htot-Evc-0.5*uu*uu),pp,r,uu
    while x<xmax:
        w=rates(Y,rho,T,Tv,mode,True)
        dY={s:w[s]*SP[s]['M']/(rho*u) for s in SPL}
        # Tv ODE (non-preferential): sum Y cvv dTv/dx = sum Y (ev*(T)-ev(Tv))/(tau u)
        num=0.0;den=0.0
        for s in ('O2','N2','NO'):
            if Y[s]<=1e-8:continue
            tau=tauVib(s,Y,rho,p,T)
            num+=Y[s]*(ev(s,T)-ev(s,Tv))/tau
            den+=Y[s]*cvv(s,Tv)
        dTv= num/(den*u) if den>1e-12 else 0.0
        mx=max(max(abs(dY[s]) for s in SPL),abs(dTv)/max(Tv,1.0)*0.05)
        dx=min(max(2e-4/mx,1e-10),xmax/60,dx*2) if mx>0 else xmax/60
        Y={s:max(0.0,Y[s]+dx*dY[s]) for s in SPL}
        ss=sum(Y.values());Y={s:Y[s]/ss for s in SPL}
        Tv=Tv+dx*dTv
        molY=sum(Y[s] for s in('O2','N2','NO'))
        Evc=mixEv(Y,Tv) if molY>1e-6 else 0.0
        T,p,rho,u=state(Y,Evc,rho)
        if molY<=1e-6:Tv=T
        Tv=min(Tv,T+2000)  # guard
        x+=dx
        if len(hist)<5000:hist.append((x,T,Tv,dict(Y)))
        n=len(hist)
        if x>1e-4 and abs(hist[-1][1]-hist[-2][1])/max(hist[-1][0]-hist[-2][0],1e-15)<1.0 and abs(T-Tv)<0.02*T:
            break
    return dict(Tf=hist[0][1],Tvf=hist[0][2],Tend=T,Tvend=Tv,x=x,hist=hist,Y=Y)
T1,p1,V=270.0,79.0,5000.0
r=relax2T(T1,p1,V,'s5')
print(f"2-T s5: frozen T2={r['Tf']:.0f} K (expect ~ideal 12350), Tv0={r['Tvf']:.0f} K")
print(f"  end: T={r['Tend']:.0f} K, Tv={r['Tvend']:.0f} K (1-T equilibrium target ~5360 K), x={r['x']*1000:.1f} mm")
print(f"  Yend: "+", ".join(f"{s}={r['Y'][s]:.4f}" for s in SPL))
# where do T and Tv merge (within 5%)?
merge=next((x for x,T,Tv,_ in r['hist'] if abs(T-Tv)<0.05*T),None)
print(f"  T-Tv merge (5%) at x={merge*1000:.2f} mm; points={len(r['hist'])}")
# peak Tv overshoot?
Tvmax=max(t[2] for t in r['hist'])
print(f"  max Tv={Tvmax:.0f} K")
