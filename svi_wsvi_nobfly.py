import argparse
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.optimize import brentq, minimize, minimize_scalar


def durrleman(k, w, wp, wpp):
    return (1 - k * wp / (2 * w))**2 - wp**2 * (1 / w + 0.25) / 4 + wpp / 2


def wsvi_curve(k, p, ell2):
    skew, d, theta, gamma = p
    k = np.asarray(k, dtype=float)
    M = k*k + ell2
    L = M * (skew*k - theta) + gamma*k*k
    delta = L*L - 4*M*M*d*k*k
    w = (-L + np.sqrt(delta)) / (2*M)
    A = 2*k*(w*w + skew*k*w - theta*w + d*k*k) + M*(skew*w + 2*d*k) + 2*gamma*k*w
    B = 2*M*w + L
    wp = -A / B
    Fkk = 2*(w*w + skew*k*w - theta*w + d*k*k) + 4*k*(skew*w + 2*d*k) + 2*M*d + 2*gamma*w
    Fkw = 2*k*(2*w + skew*k - theta) + M*skew + 2*gamma*k
    wpp = -(Fkk + 2*Fkw*wp + 2*M*wp*wp) / B
    return w, durrleman(k, w, wp, wpp)


def svi_curve(k, p):
    a, b, rho, m, sigma = p
    k = np.asarray(k, dtype=float)
    x = k-m
    h = np.hypot(x, sigma)
    w = a+b*(rho*x+h)
    wp = b*(rho+x/h)
    wpp = b*sigma*sigma/h**3
    return w, durrleman(k, w, wp, wpp)


def svi_to_conic(p):
    a, b, rho, m, sigma = p
    q = a-b*rho*m
    return np.array([1, -2*b*rho, b*b*(rho*rho-1), -2*q,
                     2*b*rho*q+2*b*b*m, q*q-b*b*(m*m+sigma*sigma)])


def conic_to_svi(y):
    z2, z3, z4, z5, z6 = y
    disc = z2*z2-4*z3
    J = z5*z5-z2*z4*z5+z3*z4*z4+disc*z6
    if disc <= 0 or J >= 0:
        return None
    b = np.sqrt(disc)/2
    rho = -z2/np.sqrt(disc)
    m = (2*z5-z2*z4)/disc
    sigma = 2*np.sqrt(-J)/disc
    a = -(z4+z2*m)/2
    return np.array([a, b, rho, m, sigma])


def polynomial_intervals(poly, eps=11):
    intervals = []
    for factor, _ in sp.factor_list(poly)[1]:
        if factor.degree() > 0:
            intervals.extend(bound for bound, _ in factor.intervals(eps=sp.Rational(1, 10**eps)))
    return sorted(set(intervals))


def rational_probes(intervals):
    if not intervals:
        return [sp.Integer(0)]
    probes = [intervals[0][0]-1]
    for left, right in zip(intervals[:-1], intervals[1:]):
        if left[1] >= right[0]:
            raise RuntimeError('Root-isolation intervals overlap')
        probes.append((left[1]+right[0])/2)
    probes.append(intervals[-1][1]+1)
    return probes


def radical_sign(a, b, delta):
    sa, sb = int(sp.sign(a)), int(sp.sign(b))
    if sa == 0:
        return sb
    if sb == 0 or sa == sb:
        return sa
    cmp = int(sp.sign(a*a-b*b*delta))
    return sa if cmp > 0 else (sb if cmp < 0 else 0)


def wsvi_oracle(p, ell2=sp.Rational(1, 140), digits=9):
    skew, d, theta, gamma = [sp.Rational(f'{float(x):.{digits}f}') for x in p]
    k, w = sp.symbols('k w')
    M = k*k+ell2
    F = sp.expand(M*(w*w+skew*k*w-theta*w+d*k*k)+gamma*k*k*w)
    A, B = sp.diff(F, k), sp.diff(F, w)
    G = sp.expand(16*w*w*B**3+16*k*w*A*B**2+(4*k*k-4*w-w*w)*A*A*B
                  -8*w*w*(sp.diff(F,k,2)*B*B-2*sp.diff(F,k,w)*A*B+sp.diff(F,w,2)*A*A))
    t0 = perf_counter()
    resultant = sp.Poly(sp.resultant(F,G,w),k)
    if resultant.is_zero:
        raise RuntimeError('Degenerate W-SVI resultant')
    intervals = polynomial_intervals(resultant)
    probes = rational_probes(intervals)
    signs = []
    for probe in probes:
        fq = sp.Poly(F.subs(k,probe),w)
        rq = sp.Poly(G.subs(k,probe),w).rem(fq)
        aa, bb, cc = fq.nth(2), fq.nth(1), fq.nth(0)
        delta = bb*bb-4*aa*cc
        if delta <= 0 or aa <= 0:
            raise RuntimeError('Upper W-SVI branch is not real')
        signs.append(radical_sign(2*aa*rq.nth(0)-bb*rq.nth(1),rq.nth(1),delta))
    root_term = float(sp.sqrt(skew*skew-4*d))
    left = (root_term+float(skew))/2
    right = (root_term-float(skew))/2
    valid = bool(theta > 0 and d < 0 and gamma >= 0 and 0 < left < 2 and 0 < right < 2)
    return {'certified':valid and all(s >= 0 for s in signs), 'degree':resultant.degree(),
            'roots':[float((a+b)/2) for a,b in intervals], 'signs':signs,
            'bad_strikes':[float(x) for x,s in zip(probes,signs) if s < 0],
            'rounded':np.array([float(skew),float(d),float(theta),float(gamma)]),
            'seconds':perf_counter()-t0}


def svi_oracle(p, digits=9):
    a, b, rho, m, sigma = [sp.Rational(f'{float(x):.{digits}f}') for x in p]
    t = sp.symbols('t')
    U, V = b*sigma*(1+rho)/2, b*sigma*(1-rho)/2
    E, W = 1+t*t, U*t*t+a*t+V
    H, R = m*t+sigma*(t*t-1)/2, U*t*t-V
    numerator = 4*E*(sigma*E*W-H*R)**2-E*R*R*W*(W+4*t)+16*(U+V)*t**3*W*W
    t0 = perf_counter()
    poly = sp.Poly(sp.expand(numerator),t)
    intervals = [v for v in polynomial_intervals(poly) if v[1] > 0]
    if any(lo <= 0 and hi > 0 for lo,hi in intervals):
        raise RuntimeError('Positive root isolation crosses zero')
    probes = []
    if intervals:
        probes.append(intervals[0][0]/2)
        for left,right in zip(intervals[:-1],intervals[1:]):
            probes.append((left[1]+right[0])/2)
        probes.append(intervals[-1][1]+1)
    else:
        probes = [sp.Integer(1)]
    signs = [int(sp.sign(poly.eval(x))) for x in probes]
    left, right = float(b*(1-rho)), float(b*(1+rho))
    wmin = float(a+b*sigma*sp.sqrt(1-rho*rho)) if abs(rho)<1 else -np.inf
    valid = bool(b>0 and sigma>0 and abs(rho)<1 and wmin>0 and 0<left<2 and 0<right<2)
    bad = [float(m+sigma*(x-1/x)/2) for x,s in zip(probes,signs) if s<0]
    return {'certified':valid and all(s >= 0 for s in signs), 'degree':poly.degree(),
            'roots':[float((lo+hi)/2) for lo,hi in intervals], 'signs':signs,
            'bad_strikes':bad, 'rounded':np.array([float(a),float(b),float(rho),float(m),float(sigma)]),
            'seconds':perf_counter()-t0}


def weighted_objective(X, weights):
    C, target = X[:,1:],X[:,0]
    factor = 1/np.sum(weights)
    def loss(y):
        r = C@y+target
        return float((weights@(r*r))*factor)
    def grad(y):
        r = C@y+target
        return 2*factor*C.T@(weights*r)
    return loss, grad


def fit_wsvi(kobs, wobs, weights, ell2, margin=0.002, rounds=12):
    M = kobs*kobs+ell2
    X = np.column_stack((M*wobs*wobs, M*kobs*wobs, M*kobs*kobs, -M*wobs, kobs*kobs*wobs))
    X /= (2*M*wobs)[:,None]
    loss, grad = weighted_objective(X,weights)
    y = np.array([0., -0.15/140, 1/140, 0.8/140])
    active = [0.]
    t0 = perf_counter()
    for j in range(rounds):
        ks = np.asarray(active)
        def structure(q):
            skew, d, theta, gamma = q
            curvature = 2*(-d/theta-gamma/ell2)
            return np.array([theta-1e-6, -d-1e-7, gamma, 4+d+2*skew-1e-4,
                             4+d-2*skew-1e-4, -curvature-1e-5])
        cons = [{'type':'ineq','fun':structure},
                {'type':'ineq','fun':lambda q, ks=ks:wsvi_curve(ks,q,ell2)[1]-margin}]
        result = minimize(loss,y,jac=grad,method='SLSQP',constraints=cons,
                          bounds=[(-0.08,0.08),(-0.01,-1e-8),(0.001,0.03),(0,0.03)],
                          options={'ftol':1e-15,'maxiter':400})
        if not result.success:
            raise RuntimeError(f'W-SVI optimizer failed: {result.message}')
        y = result.x
        cert = wsvi_oracle(y)
        if cert['certified']:
            return cert['rounded'], {'rounds':j+1,'elapsed':perf_counter()-t0,
                                     'objective':loss(y),'certificate':cert}
        missing = [v for v in cert['bad_strikes'] if not np.any(np.isclose(v,active,atol=1e-8,rtol=0))]
        if not missing:
            raise RuntimeError('W-SVI failed global certificate without new violating strikes')
        active.extend(missing)
    raise RuntimeError('W-SVI refinement limit exceeded')


def approximate_svi_min(p):
    a,b,rho,m,sigma = p
    x = np.linspace(-15,15,321)
    g = svi_curve(m+sigma*np.sinh(x),p)[1]
    ids = np.where((g[1:-1]<g[:-2])&(g[1:-1]<=g[2:]))[0]+1
    values = [(g[i],m+sigma*np.sinh(x[i])) for i in np.r_[0,ids,len(x)-1]]
    for i in ids:
        res = minimize_scalar(lambda u:svi_curve(m+sigma*np.sinh(u),p)[1],
                              bounds=(x[i-1],x[i+1]),method='bounded')
        values.append((res.fun,m+sigma*np.sinh(res.x)))
    return min(values,key=lambda pair:pair[0])


def fit_svi(kobs, wobs, weights, target, margin=0.001, rounds=12):
    D = np.column_stack((wobs*wobs,kobs*wobs,kobs*kobs,wobs,kobs,np.ones_like(kobs)))
    zfit = np.r_[1.,np.linalg.lstsq(D[:,1:],-D[:,0],rcond=None)[0]]
    B0 = 2*wobs+zfit[1]*kobs+zfit[3]
    if np.any(B0<=0):
        raise RuntimeError('Initial SVI conic denominator is not positive')
    loss,grad = weighted_objective(D/B0[:,None],weights)
    level = brentq(lambda a:approximate_svi_min(np.r_[a,target[1:]])[0]-1e-5,target[0],0)
    y = svi_to_conic(np.r_[level,target[1:]])[1:]
    active = [approximate_svi_min(target)[1]]
    t0 = perf_counter()
    def structure(q):
        z2,z3,z4,z5,z6 = q
        disc = z2*z2-4*z3
        J = z5*z5-z2*z4*z5+z3*z4*z4+disc*z6
        raw = conic_to_svi(q)
        if raw is None or abs(raw[2])>=1:
            return -np.ones(13)
        a,b,rho,m,sigma = raw
        wmin = a+b*sigma*np.sqrt(1-rho*rho)
        return np.array([-z3-1e-7, -J-1e-10, wmin-1e-5,
                         z3+4+2*z2-1e-6, z3+4-2*z2-1e-6,
                         sigma-0.18,0.85-sigma,b-0.07,0.23-b,m+0.25,0.85-m,rho+0.4,0.8-rho])
    for j in range(rounds):
        ks = np.asarray(active)
        def pointwise(q):
            raw = conic_to_svi(q)
            if raw is None or abs(raw[2])>=1 or np.any(svi_curve(ks,raw)[0]<=0):
                return -1000*np.ones(len(ks))
            return svi_curve(ks,raw)[1]-margin
        result = minimize(loss,y,jac=grad,method='SLSQP',
                          constraints=[{'type':'ineq','fun':structure},{'type':'ineq','fun':pointwise}],
                          bounds=[(-.3,.3),(-.25,-1e-7),(-.6,.6),(-.3,.3),(-.4,.4)],
                          options={'ftol':1e-11,'maxiter':180})
        if not result.success:
            raise RuntimeError(f'SVI optimizer failed: {result.message}')
        y = result.x
        raw = conic_to_svi(y)
        cert = svi_oracle(raw)
        if cert['certified']:
            return cert['rounded'],{'rounds':j+1,'elapsed':perf_counter()-t0,
                                    'objective':loss(y),'certificate':cert}
        missing = [v for v in cert['bad_strikes'] if not np.any(np.isclose(v,active,atol=1e-8,rtol=0))]
        if not missing:
            raise RuntimeError('SVI failed global certificate without new violating strikes')
        active.extend(missing)
    raise RuntimeError('SVI refinement limit exceeded')


def plot_comparison(model, kobs, wobs, target, fits, curve, outdir):
    x = np.linspace(-.3,.3,801) if model == 'wsvi' else np.linspace(-.8,1.9,901)
    original_w, original_g = curve(x,target)
    values = [(name,)+curve(x,params) for name,params,info in fits]
    fig,axes = plt.subplots(2,1,figsize=(8.5,8.5),sharex=True)
    axes[0].plot(kobs,wobs,'.',ms=3,label='Arbitrage-violating target')
    axes[1].plot(x,original_g,label='Target g(k)')
    for name,w,g in values:
        axes[0].plot(x,w,lw=1.8,label=name)
        axes[1].plot(x,g,lw=1.8,label=name)
    axes[1].axhline(0,color='black',lw=1)
    axes[1].fill_between(x,original_g,0,where=original_g<0,alpha=.15)
    axes[0].set_ylabel('Total implied variance w(k)')
    axes[1].set_ylabel('Durrleman g(k)')
    axes[1].set_xlabel('Forward log-moneyness k')
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend()
    fig.tight_layout()
    combined=outdir/f'{model}_comparison.png'
    fig.savefig(combined,dpi=180)
    plt.close(fig)
    for idx,kind in enumerate(('smile','durrleman')):
        fig,ax=plt.subplots(figsize=(8.5,4.25))
        if idx==0:
            ax.plot(kobs,wobs,'.',ms=3,label='Arbitrage-violating target')
        else:
            ax.plot(x,original_g,label='Target g(k)')
        for name,w,g in values:
            ax.plot(x,w if idx==0 else g,lw=1.8,label=name)
        if idx==1:
            ax.axhline(0,color='black',lw=1)
            ax.fill_between(x,original_g,0,where=original_g<0,alpha=.15)
        ax.set_xlabel('Forward log-moneyness k')
        ax.set_ylabel('Total implied variance w(k)' if idx==0 else 'Durrleman g(k)')
        ax.legend()
        ax.grid(alpha=.2)
        fig.tight_layout()
        fig.savefig(outdir/f'{model}_{kind}.png',dpi=180)
        plt.close(fig)
    return combined


def main():
    ap=argparse.ArgumentParser(description='Direct-fit SVI or symmetric quartic W-SVI with global butterfly-arbitrage testing')
    ap.add_argument('--model',choices=['wsvi','svi'],default='wsvi')
    ap.add_argument('--boost',type=float,default=9.)
    ap.add_argument('--width',type=float)
    ap.add_argument('--outdir',type=Path,default=Path('.'))
    args=ap.parse_args()
    width=args.width if args.width is not None else (0.065 if args.model=='wsvi' else 0.30)
    if args.boost<0 or width<=0:
        ap.error('Require --boost >= 0 and --width > 0')
    args.outdir.mkdir(parents=True,exist_ok=True)
    if args.model=='wsvi':
        ell2=1/140
        target=np.array([0.,-0.15/140,1/140,1.2/140])
        kobs=np.linspace(-.20,.20,101)
        curve=lambda k,p:wsvi_curve(k,p,ell2)
        fitter=lambda w,weights:fit_wsvi(kobs,w,weights,ell2)
        atm=0.05
    else:
        target=np.array([-0.041,0.1331,0.306,0.3586,0.4153])
        kobs=np.linspace(-.65,1.65,71)
        curve=svi_curve
        fitter=lambda w,weights:fit_svi(kobs,w,weights,target)
        atm=0.25
    wobs=curve(kobs,target)[0]
    weights=1+args.boost*np.exp(-0.5*(kobs/width)**2)
    fits=[]
    for name,ww in [('Unweighted',np.ones_like(kobs)),('ATM-weighted',weights)]:
        params,info=fitter(wobs,ww)
        err=curve(kobs,params)[0]-wobs
        mask=np.abs(kobs)<=atm
        print(f'{name}: parameters={np.array2string(params,precision=9)}')
        print(f'  variance RMSE={np.sqrt(np.mean(err**2)):.10g}; ATM RMSE={np.sqrt(np.mean(err[mask]**2)):.10g}')
        print(f"  rounds={info['rounds']}; seconds={info['elapsed']:.4f}; global no-butterfly certificate={info['certificate']['certified']}")
        print(f"  certificate polynomial degree={info['certificate']['degree']}; real-root count={len(info['certificate']['roots'])}")
        fits.append((name,params,info))
    output=plot_comparison(args.model,kobs,wobs,target,fits,curve,args.outdir)
    print(f'Plots: {output}, {args.outdir / (args.model+"_smile.png")}, {args.outdir / (args.model+"_durrleman.png")}')


if __name__=='__main__':
    main()