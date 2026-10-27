import sympy as s

t,U,V,a,m,sig=s.symbols('t U V a m sig', nonzero=True)
E=1+t*t
W=U*t*t+a*t+V
H=m*t+sig*(t*t-1)/2
R=U*t*t-V
k=H/t
w=W/t
p=s.factor(s.diff(w,t)/s.diff(k,t))
pp=s.factor(s.diff(p,t)/s.diff(k,t))
g=(1-k*p/(2*w))**2 -p**2*(1/w+s.Rational(1,4))/4 +pp/2
Pi=4*E*(sig*E*W-H*R)**2-E*R**2*W*(W+4*t)+16*(U+V)*t**3*W**2
check=s.cancel(g*(4*sig**2*E**3*W**2)-Pi)
print('wprime:',p)
print('wsecond:',pp)
print('SYMBOLIC Pi FORMULA CHECK:',check)

k,w=s.symbols('k w')
z2,z3,z4,z5,z6=s.symbols('z2 z3 z4 z5 z6')
F=w**2+z2*k*w+z3*k**2+z4*w+z5*k+z6
A=s.diff(F,k); B=s.diff(F,w)
J=z5**2-z2*z4*z5+z3*z4**2+(z2**2-4*z3)*z6
invariant=s.expand(A**2-z2*A*B+z3*B**2-(4*z3-z2**2)*F-J)
print('CONIC INVARIANT CHECK:',invariant)

# Check earlier explicit P,Q formulas along the conic, using rational numeric coefficients
av,bv,mv,rv,sv=[s.Rational(x) for x in ('-0.0410','0.1331','0.3586','0.3060','0.4153')]
d=av-bv*rv*mv
zs={z2:-2*bv*rv,z3:bv*bv*(rv*rv-1),z4:-2*d,z5:2*d*bv*rv+2*bv*bv*mv,z6:d*d-bv*bv*(mv*mv+sv*sv)}
Fn=F.subs(zs)
An=A.subs(zs);Bn=B.subs(zs);Jn=J.subs(zs)
G=16*w**2*Bn**3+16*k*w*An*Bn**2+(4*k**2-4*w-w**2)*An**2*Bn-16*w**2*(An**2-zs[z2]*An*Bn+zs[z3]*Bn**2)
L=zs[z2]*k+zs[z4]
D=s.expand(L**2-4*(zs[z3]*k**2+zs[z5]*k+zs[z6]))
Pp=s.diff(D,k)/2
T=k*k+L/2-(D+L**2)/16
P=D*(-8*L*D+4*k*(zs[z2]*D+L*Pp)-2*zs[z2]*Pp*T+(L-4)*(zs[z2]**2*D+Pp**2)/8)-4*Jn*(D+L**2)
Q=4*D*(D+L**2-k*(Pp+zs[z2]*L))+8*Jn*L+T*(zs[z2]**2*D+Pp**2)-(L-4)*zs[z2]*Pp*D/4
rem=s.rem(G,Fn,w)
claimed=s.expand(P+Q*(2*w+L))
res=s.Poly(s.expand(rem-claimed),w)
print('PREVIOUS P/Q REMAINDER CHECK:',res.is_zero)
if not res.is_zero: print('leading discrepancy:', [(monom,str(c)[:140]) for monom,c in res.terms()][:2])