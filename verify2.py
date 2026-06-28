# Verify the NEW perfect-gas physics modules before embedding in the web app.
import math

def isentropic(M, g=1.4):
    f = 1 + (g-1)/2*M*M
    T0T = f
    p0p = f**(g/(g-1))
    r0r = f**(1/(g-1))
    AAs = (1/M)*((2/(g+1))*f)**((g+1)/(2*(g-1))) if M>0 else float('inf')
    mu  = math.degrees(math.asin(1/M)) if M>=1 else None
    nu  = (math.degrees(math.sqrt((g+1)/(g-1))*math.atan(math.sqrt((g-1)/(g+1)*(M*M-1)))
           - math.atan(math.sqrt(M*M-1)))) if M>=1 else None
    return dict(T0T=T0T,p0p=p0p,r0r=r0r,AAs=AAs,mu=mu,nu=nu)

def M_from_AAs(AAs, g=1.4, branch='sup'):
    lo,hi = (1.0,50.0) if branch=='sup' else (1e-4,1.0)
    f = lambda M: isentropic(M,g)['AAs']-AAs
    for _ in range(200):
        mid=0.5*(lo+hi); v=f(mid)
        if branch=='sup':
            if v>0: hi=mid
            else: lo=mid
        else:
            if v>0: lo=mid
            else: hi=mid
    return 0.5*(lo+hi)

def M_from_nu(nu, g=1.4):
    f = lambda M: isentropic(M,g)['nu']-nu
    lo,hi=1.0,90.0
    for _ in range(200):
        mid=0.5*(lo+hi)
        if f(mid)>0: hi=mid
        else: lo=mid
    return 0.5*(lo+hi)

def normalshock(M1, g=1.4):
    M2=math.sqrt((1+(g-1)/2*M1*M1)/(g*M1*M1-(g-1)/2))
    p2p1=1+2*g/(g+1)*(M1*M1-1)
    r2r1=(g+1)*M1*M1/((g-1)*M1*M1+2)
    T2T1=p2p1/r2r1
    p02p01=(((g+1)*M1*M1/((g-1)*M1*M1+2))**(g/(g-1)))*((g+1)/(2*g*M1*M1-(g-1)))**(1/(g-1))
    return dict(M2=M2,p2p1=p2p1,r2r1=r2r1,T2T1=T2T1,p02p01=p02p01)

def theta_from_beta(M1, beta_deg, g=1.4):
    b=math.radians(beta_deg)
    Mn1=M1*math.sin(b)
    if Mn1<1: return None
    num=2/math.tan(b)*(M1*M1*math.sin(b)**2-1)
    den=M1*M1*(g+math.cos(2*b))+2
    return math.degrees(math.atan(num/den))

def oblique(M1, theta_deg, g=1.4, strong=False):
    mu=math.degrees(math.asin(1/M1))
    betas=[mu+1e-4+ i*(90-mu)/2000 for i in range(2001)]
    th=[theta_from_beta(M1,b,g) or -1 for b in betas]
    bmax=max(range(len(th)), key=lambda i:th[i]); theta_max=th[bmax]
    if theta_deg>theta_max+1e-9: return dict(detached=True, theta_max=theta_max)
    if not strong:
        lo,hi=mu, betas[bmax]
    else:
        lo,hi=betas[bmax],90
    f=lambda b:(theta_from_beta(M1,b,g) or -90)-theta_deg
    for _ in range(200):
        mid=0.5*(lo+hi)
        if not strong:
            if f(mid)>0: hi=mid
            else: lo=mid
        else:
            if f(mid)>0: lo=mid
            else: hi=mid
    beta=0.5*(lo+hi); b=math.radians(beta)
    Mn1=M1*math.sin(b); ns=normalshock(Mn1,g)
    M2=ns['M2']/math.sin(b-math.radians(theta_deg))
    return dict(beta=beta, theta_max=theta_max, M2=M2, **{k:ns[k] for k in('p2p1','r2r1','T2T1','p02p01')})

def Vbar_from_M(M,g=1.4):
    return 1/math.sqrt(1+2/((g-1)*M*M))
def M_from_Vbar(V,g=1.4):
    return math.sqrt(2*V*V/((g-1)*(1-V*V)))
def taylor_maccoll_sigma(M1, beta_deg, g=1.4):
    b=math.radians(beta_deg)
    delta=theta_from_beta(M1,beta_deg,g)
    if delta is None or delta<0: return None
    d=math.radians(delta)
    ns=normalshock(M1*math.sin(b),g)
    M2=ns['M2']/math.sin(b-d)
    V=Vbar_from_M(M2,g)
    Vr=V*math.cos(b-d); Vth=-V*math.sin(b-d)
    th=b; dth=-1e-4
    while th>1e-4:
        def deriv(theta,Vr,Vth):
            denom=(g-1)/2*(1-Vr*Vr-Vth*Vth)-Vth*Vth
            dVr=Vth
            dVth=( Vth*Vth*Vr - (g-1)/2*(1-Vr*Vr-Vth*Vth)*(2*Vr+Vth/math.tan(theta)) )/denom
            return dVr,dVth
        k1=deriv(th,Vr,Vth)
        k2=deriv(th+dth/2,Vr+dth/2*k1[0],Vth+dth/2*k1[1])
        k3=deriv(th+dth/2,Vr+dth/2*k2[0],Vth+dth/2*k2[1])
        k4=deriv(th+dth,Vr+dth*k3[0],Vth+dth*k3[1])
        Vr_n=Vr+dth/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        Vth_n=Vth+dth/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        if Vth_n>=0:
            frac=Vth/(Vth-Vth_n)
            th_s=th+frac*dth
            Vr_s=Vr+frac*(Vr_n-Vr)
            return dict(sigma=math.degrees(th_s), Msurf=M_from_Vbar(Vr_s,g), M2=M2)
        Vr,Vth,th=Vr_n,Vth_n,th+dth
    return None
def cone(M1, sigma_deg, g=1.4):
    lo=math.degrees(math.asin(1/M1))+0.05; hi=70
    def s(b):
        r=taylor_maccoll_sigma(M1,b,g); return (r['sigma'] if r else -90)
    for _ in range(80):
        mid=0.5*(lo+hi)
        if s(mid)>sigma_deg: hi=mid
        else: lo=mid
    beta=0.5*(lo+hi); r=taylor_maccoll_sigma(M1,beta,g)
    return dict(beta=beta, Msurf=r['Msurf'], M2behindshock=r['M2'])

print("== ISENTROPIC (g=1.4) ==")
for M in [0.5,2.0]:
    s=isentropic(M); print(f" M={M}: p0/p={s['p0p']:.4f} A/A*={s['AAs']:.4f} nu={s['nu']} mu={s['mu']}")
print(" inverse A/A*=1.6875 sup ->", round(M_from_AAs(1.6875,branch='sup'),4),"(expect 2.0)")
print(" inverse A/A*=1.6875 sub ->", round(M_from_AAs(1.6875,branch='sub'),4),"(expect ~0.3722)")
print(" inverse nu=26.3798 ->", round(M_from_nu(26.3798),4),"(expect 2.0)")
print("== NORMAL SHOCK M1=2 ==", {k:round(v,4) for k,v in normalshock(2.0).items()})
print("== OBLIQUE M1=2 theta=20 (weak) ==")
o=oblique(2.0,20.0); print(" ", {k:(round(v,4) if isinstance(v,float) else v) for k,v in o.items()},"(beta~53.4, M2~1.21)")
os_=oblique(2.0,20.0,strong=True); print(" strong beta=",round(os_['beta'],2),"M2=",round(os_['M2'],3))
print(" theta_max(M=2)=",round(o['theta_max'],2),"(expect ~22.97)")
print("== CONE M1=2 sigma=20 ==")
c=cone(2.0,20.0); print(" ", {k:round(v,3) for k,v in c.items()},"(beta~37-38, Msurf~1.6)")
c5=cone(5.0,15.0); print(" M1=5 sigma=15:", {k:round(v,3) for k,v in c5.items()})
