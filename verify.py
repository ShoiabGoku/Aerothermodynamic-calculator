import math
Ru=8.314462618;kB=1.380649e-23;hPl=6.62607015e-34;NA=6.02214076e23
R_air=287.05;GAM=1.4;Tref=298.15
SP={'O2':dict(M=0.0319988,mol=True,thv=2273.6,thr=2.079,sig=2,ge=3,hf=0),
    'N2':dict(M=0.0280134,mol=True,thv=3393.5,thr=2.884,sig=2,ge=1,hf=0),
    'NO':dict(M=0.0300061,mol=True,thv=2739.0,thr=2.452,sig=1,ge=4,hf=90290),
    'O':dict(M=0.0159994,mol=False,ge=9,hf=249170),
    'N':dict(M=0.0140067,mol=False,ge=4,hf=472680)}
D_O2=2*SP['O']['hf']-SP['O2']['hf'];D_N2=2*SP['N']['hf']-SP['N2']['hf'];D_NO=SP['N']['hf']+SP['O']['hf']-SP['NO']['hf']
xO2=0.21;xN2=0.79;RATIO=(2*xN2)/(2*xO2)
def hS(s,T):
    d=SP[s]
    if not d['mol']: return d['hf']+2.5*Ru*(T-Tref)
    vib=lambda th:d['thv']/(math.exp(d['thv']/th)-1)
    return d['hf']+3.5*Ru*(T-Tref)+Ru*(vib(T)-vib(Tref))
def cpS(s,T):
    d=SP[s]
    if not d['mol']: return 2.5*Ru
    x=d['thv']/T;ex=math.exp(x);return 3.5*Ru+Ru*x*x*ex/((ex-1)**2)
def qVol(s,T):
    d=SP[s];m=d['M']/NA;qt=(2*math.pi*m*kB*T/(hPl*hPl))**1.5
    if not d['mol']: return qt*d['ge']
    return qt*(T/(d['sig']*d['thr']))*(1/(1-math.exp(-d['thv']/T)))*d['ge']
def Kp_O2(T): return qVol('O',T)**2/qVol('O2',T)*kB*T*math.exp(-D_O2/(NA*kB*T))
def Kp_N2(T): return qVol('N',T)**2/qVol('N2',T)*kB*T*math.exp(-D_N2/(NA*kB*T))
def Kp_NO(T): return qVol('N',T)*qVol('O',T)/qVol('NO',T)*kB*T*math.exp(-D_NO/(NA*kB*T))
def bisect(f,a,b,n):
    fa=f(a);fb=f(b);k=0
    while fa*fb>0 and k<60: b*=1.6;fb=f(b);k+=1
    for _ in range(n):
        m=0.5*(a+b);fm=f(m)
        if fa*fm<=0:b=m;fb=fm
        else:a=m;fa=fm
    return 0.5*(a+b)
def equil(T,p,mode):
    KO=Kp_O2(T);r=RATIO;pO=pN=pNO=0.0
    if mode=='s3':
        f=lambda po:(po*po/KO)+po+r*(2*po*po/KO+po)/2-p
        pO=bisect(f,0,p,80);pO2=pO*pO/KO;pN2=r*(2*pO2+pO)/2
    else:
        KN=Kp_N2(T);KNO=Kp_NO(T) if mode=='s5' else float('inf')
        def comp(po,pn):
            po2=po*po/KO;pn2=pn*pn/KN;pno=0 if KNO==float('inf') else po*pn/KNO
            return po2,pn2,pno
        def F(po,pn):
            po2,pn2,pno=comp(po,pn)
            dalton=po2+pn2+po+pn+pno-p
            Nat=2*pn2+pn+pno;Oat=2*po2+po+pno
            return dalton,Nat-r*Oat
        x=max(1e-3*p,1e-4);y=max(1e-3*p,1e-4)
        for _ in range(120):
            f0=F(x,y);dx=max(x*1e-6,1e-9);dy=max(y*1e-6,1e-9)
            fx=F(x+dx,y);fy=F(x,y+dy)
            J=[[(fx[0]-f0[0])/dx,(fy[0]-f0[0])/dy],[(fx[1]-f0[1])/dx,(fy[1]-f0[1])/dy]]
            det=J[0][0]*J[1][1]-J[0][1]*J[1][0]
            if abs(det)<1e-300:break
            sx=(-f0[0]*J[1][1]+f0[1]*J[0][1])/det;sy=(-f0[1]*J[0][0]+f0[0]*J[1][0])/det
            lam=1.0
            while lam>1e-4 and (x+lam*sx<=0 or y+lam*sy<=0):lam*=0.5
            x+=lam*sx;y+=lam*sy;x=max(x,1e-30);y=max(y,1e-30)
            if abs(sx*lam)<1e-6*max(x,1) and abs(sy*lam)<1e-6*max(y,1):break
        pO=x;pN=y;pO2,pN2,pNO=comp(pO,pN)
    parts={'O2':pO2,'N2':pN2,'O':pO,'N':pN,'NO':pNO};psum=sum(parts.values())
    xf={s:parts[s]/psum for s in parts};Mbar=sum(xf[s]*SP[s]['M'] for s in parts)
    Y={s:xf[s]*SP[s]['M']/Mbar for s in parts}
    rho=p*Mbar/(Ru*T);h=sum(Y[s]*hS(s,T)/SP[s]['M'] for s in parts)
    return dict(x=xf,Y=Y,Mbar=Mbar,rho=rho,h=h)
def frozenAir(T,p):
    Mbar=xO2*SP['O2']['M']+xN2*SP['N2']['M']
    YO2=xO2*SP['O2']['M']/Mbar;YN2=xN2*SP['N2']['M']/Mbar
    return dict(Mbar=Mbar,rho=p*Mbar/(Ru*T),h=YO2*hS('O2',T)/SP['O2']['M']+YN2*hS('N2',T)/SP['N2']['M'],
                Y=dict(O2=YO2,N2=YN2,O=0,N=0,NO=0))
def frozenCp(T):
    Mbar=xO2*SP['O2']['M']+xN2*SP['N2']['M']
    return (xO2*SP['O2']['M']/Mbar)*cpS('O2',T)/SP['O2']['M']+(xN2*SP['N2']['M']/Mbar)*cpS('N2',T)/SP['N2']['M']
def solveT(targetH,p,mode):
    hf=lambda T:(frozenAir(T,p)['h'] if mode in('real','frozen') else equil(T,p,mode)['h'])
    Tlo,Thi=200,30000;g=lambda T:hf(T)-targetH;glo=g(Tlo);ghi=g(Thi)
    if glo*ghi>0:return Thi if targetH>g(Thi) else Tlo
    for _ in range(90):
        Tm=0.5*(Tlo+Thi);gm=g(Tm)
        if glo*gm<=0:Thi=Tm;ghi=gm
        else:Tlo=Tm;glo=gm
    return 0.5*(Tlo+Thi)
def shockIdeal(M1,T1,p1):
    g=GAM;M2=math.sqrt((1+(g-1)/2*M1*M1)/(g*M1*M1-(g-1)/2))
    p2p1=1+2*g/(g+1)*(M1*M1-1);r2r1=(g+1)*M1*M1/((g-1)*M1*M1+2)
    p02p1=((g+1)**2*M1*M1/(4*g*M1*M1-2*(g-1)))**(g/(g-1))*((1-g+2*g*M1*M1)/(g+1))
    return dict(M2=M2,p2=p2p1*p1,T2=p2p1/r2r1*T1,rhoRatio=r2r1,p02=p02p1*p1)
def shockReacting(p1,T1,rho1,h1,u1,mode):
    eps=0.12;T2=4000;p2=rho2=0;comp=None
    for _ in range(200):
        p2=p1+rho1*u1*u1*(1-eps);h2=h1+0.5*u1*u1*(1-eps*eps);T2=solveT(h2,p2,mode)
        comp=frozenAir(T2,p2) if mode=='real' else equil(T2,p2,mode);rho2=comp['rho']
        epsNew=rho1/rho2;d=epsNew-eps;eps+=0.6*d
        if abs(d)<1e-7:break
    if mode=='real':cp2=frozenCp(T2)
    else:cp2=sum(comp['Y'][s]*cpS(s,T2)/SP[s]['M'] for s in comp['Y'] if comp['Y'][s]>0)
    Rmix=Ru/comp['Mbar'];g2=cp2/(cp2-Rmix);a2=math.sqrt(g2*Rmix*T2);u2=eps*u1
    return dict(T2=T2,p2=p2,rho2=rho2,comp=comp,M2=u2/a2,g2=g2,eps=eps)
def muSuth(T):return 1.716e-5*(T/273.15)**1.5*(273.15+110.4)/(T+110.4)
def suttonGraves(rho,V,Rn):return 1.7415e-4*math.sqrt(rho/Rn)*V**3
def fayRiddell(T0,p0,h0,rho0,fs_p,Rn,Tw):
    Pr=0.71;hw=frozenAir(Tw,p0)['h'];rhoe=rho0;mue=muSuth(T0)
    rhow=p0*(xO2*SP['O2']['M']+xN2*SP['N2']['M'])/(Ru*Tw);muw=muSuth(Tw)
    dudx=(1/Rn)*math.sqrt(2*max(p0-fs_p,1)/rhoe)
    return 0.763*Pr**-0.6*(rhoe*mue)**0.4*(rhow*muw)**0.1*math.sqrt(dudx)*(h0-hw)

T,p,V,Rn=270,79,5000,0.025
a=math.sqrt(GAM*R_air*T);M=V/a;rho=p/(R_air*T);h1=frozenAir(T,p)['h'];h0=h1+0.5*V*V
print(f"FREESTREAM M={M:.2f} a={a:.1f} rho={rho:.3e} h0={h0/1e6:.2f}MJ/kg Re/m={rho*V/muSuth(T):.2e}")
print(f"Sutton-Graves q={suttonGraves(rho,V,Rn)/1e6:.2f} MW/m2  (target ~4.4)\n")
idg=shockIdeal(M,T,p);T0id=T*(1+0.2*M*M)
qfr_id=fayRiddell(T0id,idg['p02'],h0,idg['p02']/(R_air*T0id),p,Rn,300)
print(f"IDEAL T2={idg['T2']:.0f} p2={idg['p2']/1000:.1f}kPa rho2/rho={idg['rhoRatio']:.2f} M2={idg['M2']:.3f} T0={T0id:.0f} p02={idg['p02']/1000:.1f}kPa qFR={qfr_id/1e6:.2f}MW/m2")
for mode in ['real','s3','s4','s5']:
    sh=shockReacting(p,T,rho,h1,V,mode);g2=sh['g2'];M2=sh['M2']
    p0=sh['p2']*(1+(g2-1)/2*M2*M2)**(g2/(g2-1));T0=solveT(h0,p0,mode if mode=='real' else mode)
    comp0=frozenAir(T0,p0) if mode=='real' else equil(T0,p0,mode)
    qfr=fayRiddell(T0,p0,h0,comp0['rho'],p,Rn,300)
    cc=None if mode=='real' else sh['comp']
    line=f"{mode.upper():5} T2={sh['T2']:.0f} p2={sh['p2']/1000:.1f}kPa rho2/rho={sh['rho2']/rho:.2f} M2={sh['M2']:.3f} T0={T0:.0f} p0={p0/1000:.1f}kPa qFR={qfr/1e6:.2f}MW/m2"
    if cc:line+=f"  x(O)={cc['x']['O']*100:.1f}% x(N)={cc['x'].get('N',0)*100:.2f}% x(NO)={cc['x'].get('NO',0)*100:.2f}%"
    print(line)
c5=equil(6000,1e5,'s5');Nat=2*c5['x']['N2']+c5['x']['N']+c5['x']['NO'];Oat=2*c5['x']['O2']+c5['x']['O']+c5['x']['NO']
print(f"\nElement check s5(6000K,1bar): N/O={Nat/Oat:.3f} (target 3.762) sumX={sum(c5['x'].values()):.4f}")
