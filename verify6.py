# Verify inverse US Standard Atmosphere: pressure/density/temperature -> altitude.
import math
def atmosphere(zkm):
    g0=9.80665; Rs=287.0528; r0=6356.766
    h=r0*zkm/(r0+zkm)
    Hb=[0,11,20,32,47,51,71,84.852]
    Tb=[288.15,216.65,216.65,228.65,270.65,270.65,214.65,186.946]
    Lb=[-6.5,0,1,2.8,0,-2.8,-2,0]
    Pb=[101325,22632.06,5474.889,868.0187,110.9063,66.93887,3.956420,0.3733836]
    i=0
    while i<7 and not (h<Hb[i+1]): i+=1
    dH=h-Hb[i]; T=Tb[i]+Lb[i]*dH
    if abs(Lb[i])<1e-9: p=Pb[i]*math.exp(-g0*1000*dH/(Rs*Tb[i]))
    else: p=Pb[i]*((Tb[i]+Lb[i]*dH)/Tb[i])**(-g0*1000/(Rs*Lb[i]))
    return T,p
R_air=287.05
def prop(z,which):
    T,p=atmosphere(z)
    if which=='p': return p
    if which=='T': return T
    if which=='rho': return p/(R_air*T)
def invert_mono(target,which):     # p or rho: monotonic decreasing
    lo,hi=0.0,86.0
    if prop(lo,which)<target or prop(hi,which)>target: return None
    for _ in range(100):
        m=0.5*(lo+hi)
        if prop(m,which)>target: lo=m
        else: hi=m
    return 0.5*(lo+hi)
def invert_T(Ttar):
    sols=[]; plateaus=[]
    zs=[i*0.02 for i in range(0,4301)]
    prevT=atmosphere(zs[0])[0]; runstart=None
    for k in range(1,len(zs)):
        z=zs[k]; T=atmosphere(z)[0]; pz=zs[k-1]
        # crossing (strict sign change)
        if (prevT-Ttar)*(T-Ttar)<0:
            zc=pz+(Ttar-prevT)/(T-prevT)*(z-pz); sols.append(zc)
        prevT=T
    # plateau detection: contiguous band within 0.05 K
    band=None
    for z in zs:
        T=atmosphere(z)[0]
        if abs(T-Ttar)<0.06:
            if band is None: band=[z,z]
            else: band[1]=z
        else:
            if band and band[1]-band[0]>0.3: plateaus.append(tuple(band))
            band=None
    if band and band[1]-band[0]>0.3: plateaus.append(tuple(band))
    return sols,plateaus

print("== pressure inversion ==")
for zt in [0,11,30,50,70]:
    _,p=atmosphere(zt); zi=invert_mono(p,'p')
    print(f"  p at {zt} km = {p:.3f} Pa -> altitude {zi:.3f} km  (err {abs(zi-zt)*1000:.1f} m)")
print("== density inversion ==")
_,p=atmosphere(40); rho=p/(R_air*atmosphere(40)[0]); print(f"  rho at 40 km -> {invert_mono(rho,'rho'):.3f} km")
print("== temperature inversion (multi-valued) ==")
for Tt in [220,250,270.65,216.65,200,288.15]:
    sols,plat=invert_T(Tt)
    ss=", ".join(f"{s:.2f}" for s in sols) or "-"
    pp=", ".join(f"{a:.1f}-{b:.1f}" for a,b in plat) or "-"
    print(f"  T={Tt} K -> crossings [{ss}] km ; isothermal bands [{pp}] km")
