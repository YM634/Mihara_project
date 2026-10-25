#!/usr/bin/env python3
"""Generate/run independent SC_model jobs; A is the applied scalar field."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import time

DEFAULT_EXE = '/Users/yutamurakami/Documents/Softwares/Lib_TDDM_ver2/Execute_sc/sc_model.ex'

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe', default=DEFAULT_EXE)
    p.add_argument('--output', type=Path, default=Path(__file__).resolve().parent / 'jobs')
    p.add_argument('--L', type=int, default=50)
    p.add_argument('--dt', type=float, default=0.02)
    p.add_argument('--tmin', type=float, default=-700)
    p.add_argument('--tmax', type=float, default=700)
    p.add_argument('--F0', type=float, nargs='+', default=[0.001, 0.01])
    p.add_argument('--epsilon', type=float, default=0.001)
    p.add_argument('--omega', type=float, nargs='+', default=[1.0])
    p.add_argument('--tau-f', type=float, default=0.1)
    p.add_argument('--tau-g', type=float, default=math.sqrt(2)/0.01)
    p.add_argument('--polarization', type=float, nargs=2, default=[1, 0], metavar=('PX', 'PY'))
    p.add_argument('--run', action='store_true')
    a = p.parse_args()
    if a.L < 2 or a.dt <= 0 or not a.tmin < 0 < a.tmax or min(a.F0) <= 0 or a.epsilon <= 0 or min(a.tau_f,a.tau_g) <= 0:
        p.error('invalid grid, amplitudes or widths')
    nt = round((a.tmax-a.tmin)/a.dt)
    if nt < 4 or not math.isclose(nt*a.dt, a.tmax-a.tmin, abs_tol=1e-8):
        p.error('time interval must be an integer number of dt (at least 4)')
    polnorm = math.hypot(*a.polarization)
    if polnorm == 0:
        p.error('polarization must be nonzero')
    pol = [v/polnorm for v in a.polarization]
    exe = str(Path(a.exe).expanduser().resolve())
    a.output.mkdir(parents=True, exist_ok=True)
    # A(0)=0 is assumed by equilibrium initialization. Do not silently offset A.
    if max(math.exp(-a.tmin**2/(2*tau*tau)) for tau in [a.tau_f,a.tau_g]) > 1e-5:
        p.error('tmin is too close to the pulse: require envelope(tmin) <= 1e-5')
    specs = [(0.,0.,'zero')]
    for f0 in a.F0:
        specs.append((f0,0.,'f'))
        for omega in a.omega:
            specs.extend((f0,omega,kind) for kind in ['cos_plus','cos_minus','sin_plus','sin_minus'])
    for idx,(f0,omega,kind) in enumerate(specs):
        d = a.output / f'{idx:04d}_{kind}_F{f0:.8g}_w{omega:.8g}'
        d.mkdir(exist_ok=False)
        handles = [open(d/name,'w') for name in ['Ax_t.txt','Ay_t.txt','Ex_t.txt','Ey_t.txt']]
        try:
            for it in range(nt+1):
                t = a.tmin+it*a.dt
                f = f0/(2*math.pi)*math.exp(-t*t/(2*a.tau_f*a.tau_f))
                df = -t/(a.tau_f*a.tau_f)*f
                g = dg = 0.
                if kind.startswith(('cos','sin')):
                    env = f0/(2*math.pi)*math.exp(-t*t/(2*a.tau_g*a.tau_g))
                    denv = -t/(a.tau_g*a.tau_g)*env
                    if kind.startswith('cos'):
                        g = env*math.cos(omega*t)
                        dg = denv*math.cos(omega*t)-omega*env*math.sin(omega*t)
                    else:
                        g = env*math.sin(omega*t)
                        dg = denv*math.sin(omega*t)+omega*env*math.cos(omega*t)
                sign = -1 if kind.endswith('minus') else 1
                aval = f+sign*a.epsilon*g
                eval_ = -df-sign*a.epsilon*dg
                for h,v in zip(handles,[pol[0]*aval,pol[1]*aval,pol[0]*eval_,pol[1]*eval_]):
                    h.write(f'{v:.17e}\n')
        finally:
            for h in handles: h.close()
        params = dict(T=0.02,mu=0.0,U_el=-4.0,g_elph=0.0,w0_ph=0.3,SC_symmetry='S',elph_type='S',L=a.L,Nt=nt,dt=a.dt,Ax_t='Ax_t.txt',Ay_t='Ay_t.txt',Ex_t='Ex_t.txt',Ey_t='Ey_t.txt',it_obs=0)
        (d/'in_params').write_text(''.join(f'__{k}= {v}\n' for k,v in params.items()))
        meta = dict(params,kind=kind,F0=f0,omega=omega,epsilon=a.epsilon,tmin=a.tmin,tmax=a.tmax,tau_f=a.tau_f,tau_g=a.tau_g,polarization=pol,exe=exe,field='vector potential',gaussian='exp(-t^2/(2 tau^2))')
        (d/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
        # Portable local job, also usable as the payload of a scheduler job.
        import shlex
        (d/'job.sh').write_text('#!/bin/sh\nset -eu\ncd "$(dirname "$0")"\nexec '+shlex.quote(exe)+' in_params > log 2> stderr\n')
        if a.run:
            start = time.perf_counter()
            with open(d/'log','w') as log, open(d/'stderr','w') as err:
                result = subprocess.run([exe,'in_params'],cwd=d,stdout=log,stderr=err)
            timing = dict(seconds=time.perf_counter()-start,returncode=result.returncode)
            (d/'timing.json').write_text(json.dumps(timing,indent=2)+'\n')
            if result.returncode or not (d/'obs_t.dat').exists():
                raise RuntimeError(f'job failed: {d}; see log/stderr')
            print(d.name, timing, flush=True)
        else:
            print(d)

if __name__ == '__main__': main()
