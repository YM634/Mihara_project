"""Export the omega1=1.40 set as five-column text for Mihara."""
from pathlib import Path
import hashlib
import numpy as np
try: import tomllib
except ModuleNotFoundError: import tomli as tomllib
root=Path(__file__).resolve().parents[2]
src=root/'Attractive_Hubbard/results/scan_x_float64'
out=root/'For_Mihara'
out.mkdir(exist_ok=True)
if not (out/"README.md").is_file():
 raise RuntimeError("For_Mihara/README.md (shared definitions) is required before export")
checks=[]
for amp in ['0.001','0.01']:
 for kind in ['cos_plus','cos_minus','sin_plus','sin_minus']:
  original=src/f'omega_028_F{amp}_{kind}'
  c=tomllib.loads((original/'config.toml').read_text())
  assert float(c['omega'])==1.4
  dest=out/f'omega1_1.40_F_{amp}_{kind}';dest.mkdir(exist_ok=True)
  datafile=dest/'observables.dat';params=dest/'Params.md'
  if datafile.exists() or params.exists():raise RuntimeError(f'Refusing overwrite: {dest}')
  n=140001
  raw=np.memmap(original/'observables.f64',dtype='<f8',mode='r',shape=(n,10),order='F')
  data=np.asarray(raw[:,[0,3,4,5,6]])
  np.savetxt(datafile,data,fmt='%.16e',header='1: time  2: Phi_re  3: Phi_im  4: J_x  5: J_y')
  loaded=np.loadtxt(datafile)
  if not np.array_equal(loaded,data):raise RuntimeError(f'Export roundtrip mismatch: {dest}')
  trig='cos' if kind.startswith('cos') else 'sin'
  sign='+' if kind.endswith('plus') else '-'
  f=float(amp);gcoef=f*.001/(2*np.pi)
  text=f'''# 個別計算条件

共通の外場定義、模型、ハミルトニアン、時間発展、データ列は[README.md](../README.md)を参照してください。

条件名: {dest.name}

F = {amp}

kind = {kind}

パルスfの係数 F/(2*pi) = {f/(2*np.pi):.16e}

テスト波epsilon*gの係数 epsilon*F/(2*pi) = {gcoef:.16e}

この条件の符号は {sign}。

元データ: Attractive_Hubbard/results/scan_x_float64/{original.name}/observables.f64
'''
  params.write_text(text)
  checks.append(f'{hashlib.sha256(datafile.read_bytes()).hexdigest()}  {dest.name}/observables.dat')
  print(dest.name,'verified',loaded.shape,flush=True)
(out/'SHA256SUMS').write_text('\n'.join(checks)+'\n')
print('total MB',sum(p.stat().st_size for p in out.rglob('*') if p.is_file())/1e6)
