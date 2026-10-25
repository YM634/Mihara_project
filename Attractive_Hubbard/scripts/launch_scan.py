"""Launch two persistent single-threaded workers; completed cases are restartable."""
import os, sys, subprocess, pathlib, json, datetime, shutil, argparse
root = pathlib.Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('output', nargs='?', type=pathlib.Path, default=root/'results/scan_x_float64')
parser.add_argument('--config', type=pathlib.Path, default=root/'configs/scan_x.toml')
args = parser.parse_args()
config = args.config.resolve()
out = args.output.resolve()
out.mkdir(parents=True, exist_ok=True)
manifest = out / 'launch.json'
if manifest.exists():
    for pid in json.loads(manifest.read_text())['pids']:
        try: os.kill(pid, 0)
        except ProcessLookupError: pass
        else: raise SystemExit(f'Previous worker PID {pid} still exists; refusing duplicate launch')
saved = out / 'scan_config.toml'
if saved.exists() and saved.read_bytes()!=config.read_bytes():
    raise SystemExit('Config differs from existing scan')
if config != saved:
    shutil.copyfile(config, saved)
env = dict(os.environ, JULIA_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', JULIA_DEPOT_PATH=str(root.parent/'Ono論文再現/.julia_depot'))
pids=[]
for worker in (1,2):
    with (out/f'worker{worker}.log').open('ab') as log:
        p=subprocess.Popen(['/opt/local/bin/julia','--startup-file=no',f'--project={root}',str(root/'scripts/scan_worker.jl'),str(saved),str(out),str(worker),'2'],env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    pids.append(p.pid)
manifest.write_text(json.dumps(dict(pids=pids,started_at=datetime.datetime.now().isoformat(),workers=2,threads_per_worker=1,expected_cases=810),indent=2)+'\n')
print(manifest.read_text())
