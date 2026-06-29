# Verify shock-tube / shock-tunnel physics before embedding.
import math

def shock_tube_p2p1(p4p1, g1, g4, a1, a4):
    # Solve the shock-tube equation for p2/p1 given diaphragm ratio p4/p1.
    def f(x):  # x = p2/p1
        term = 1 - (g4-1)*(a1/a4)*(x-1)/math.sqrt(2*g1*(2*g1+(g1+1)*(x-1)))
        if term<=0: return 1e9
        return x*term**(-2*g4/(g4-1)) - p4p1
    lo,hi=1.0000001,p4p1
    for _ in range(200):
        m=0.5*(lo+hi)
        if f(m)>0: hi=m
        else: lo=m
    return 0.5*(lo+hi)

def incident(p2p1,g1,a1,T1,rho1):
    Ms=math.sqrt((g1+1)/(2*g1)*(p2p1-1)+1)
    T2T1=(2*g1*Ms*Ms-(g1-1))*((g1-1)*Ms*Ms+2)/((g1+1)**2*Ms*Ms)
    rho2rho1=(g1+1)*Ms*Ms/((g1-1)*Ms*Ms+2)
    u2=(2*a1/(g1+1))*(Ms-1/Ms)             # mass/contact velocity (lab)
    Ws=Ms*a1                                # incident shock speed (lab)
    a2=a1*math.sqrt(T2T1)
    return dict(Ms=Ms,Ws=Ws,T2T1=T2T1,rho2rho1=rho2rho1,u2=u2,a2=a2,
                p2=p2p1*1, T2=T1*T2T1, rho2=rho1*rho2rho1)

def reflected(u2,a2,g1,T2,p2):
    K=u2*(g1+1)/(2*a2)
    Mr=(K+math.sqrt(K*K+4))/2               # closed form (u5=0)
    p5p2=1+2*g1/(g1+1)*(Mr*Mr-1)
    T5T2=(2*g1*Mr*Mr-(g1-1))*((g1-1)*Mr*Mr+2)/((g1+1)**2*Mr*Mr)
    Ur=Mr*a2-u2                              # reflected shock speed (toward driver)
    return dict(Mr=Mr,p5p2=p5p2,T5T2=T5T2,Ur=Ur,p5=p2*p5p2,T5=T2*T5T2)

def region3(p2,p4,T4,g4,a4):
    p3=p2; T3=T4*(p3/p4)**((g4-1)/g4); a3=a4*(p3/p4)**((g4-1)/(2*g4))
    u3=(2*a4/(g4-1))*(1-(p3/p4)**((g4-1)/(2*g4)))
    return dict(p3=p3,T3=T3,a3=a3,u3=u3)

# ---- Case A: air/air, p4/p1=10, T1=T4=300 ----
g=1.4; R=287.05; T1=T4=300.0; p1=1e4; p4=1e5  # p4/p1=10
a1=math.sqrt(g*R*T1); a4=math.sqrt(g*R*T4); rho1=p1/(R*T1)
p2p1=shock_tube_p2p1(p4/p1,g,g,a1,a4)
inc=incident(p2p1,g,a1,T1,rho1)
r3=region3(inc['T2']/T1*p1* (p2p1/(inc['T2']/T1)),0,0,0,0) if False else region3(p2p1*p1,p4,T4,g,a4)
ref=reflected(inc['u2'],inc['a2'],g,inc['T2'],p2p1*p1)
print("CASE A  air/air p4/p1=10, T=300K, a1=%.1f"%a1)
print(" p2/p1=%.4f  Ms=%.4f  Ws=%.1f m/s  u2(contact)=%.1f m/s"%(p2p1,inc['Ms'],inc['Ws'],inc['u2']))
print(" T2=%.1f K  region3: p3=p2? %.1f vs %.1f  u3=%.1f (should=u2=%.1f)  T3=%.1f"%(
      inc['T2'],r3['p3'],p2p1*p1,r3['u3'],inc['u2'],r3['T3']))
print(" reflected: Mr=%.4f  p5=%.1f Pa (p5/p1=%.2f)  T5=%.1f K (T5/T1=%.3f)  Ur=%.1f m/s"%(
      ref['Mr'],ref['p5'],ref['p5']/p1,ref['T5'],ref['T5']/T1,ref['Ur']))

# ---- Case B: cross-check incident Ms=2 -> T5/T1 ----
g=1.4;a1=347.2;T1=300.0;rho1=1.0
p2p1=1+2*g/(g+1)*(2.0**2-1)
inc=incident(p2p1,g,a1,T1,rho1); ref=reflected(inc['u2'],inc['a2'],g,inc['T2'],p2p1)
print("\nCASE B  Ms=2 check: p2/p1=%.3f T2/T1=%.4f u2=%.1f a2=%.1f Mr=%.4f T5/T1=%.3f (expect ~2.5)"%(
      p2p1,inc['T2T1'],inc['u2'],inc['a2'],ref['Mr'],ref['T5']/T1))

# ---- Case C: helium driver, strong shock ----
g1=1.4;R1=287.05;g4=1.667;R4=2077.0;T1=300;T4=300;p1=1e3;p4=2e6
a1=math.sqrt(g1*R1*T1);a4=math.sqrt(g4*R4*T4);rho1=p1/(R1*T1)
p2p1=shock_tube_p2p1(p4/p1,g1,g4,a1,a4); inc=incident(p2p1,g1,a1,T1,rho1)
ref=reflected(inc['u2'],inc['a2'],g1,inc['T2'],p2p1*p1)
print("\nCASE C  He driver -> air, p4/p1=2000, a4=%.0f (He)  -> p2/p1=%.2f Ms=%.3f T5=%.0fK"%(
      a4,p2p1,inc['Ms'],ref['T5']))
