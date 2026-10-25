"""Build an instrumented C++ reference; leave all original files untouched."""
from pathlib import Path
import argparse,hashlib,json,subprocess
p=argparse.ArgumentParser(description=__doc__)
p.add_argument("library_root",nargs="?",default="/Users/yutamurakami/Documents/Softwares/Lib_TDDM_ver2")
p.add_argument("--fixed",action="store_true")
args=p.parse_args()
ROOT=Path(__file__).resolve().parents[1]
LIB=Path(args.library_root)
folder=LIB/'Programs_MF/SC_phonon'
src=folder/('sc_dynamics_main_fixMF.cpp' if args.fixed else 'SC_dynamics_main.cpp')
header=src.stem+'.hpp'
s=src.read_text().replace(f'#include "{header}"',f'#include "{header}"\n#include <chrono>')
observable='evaluate_observales_fixMF' if args.fixed else 'evaluate_observales'
s=s.replace('setprecision(12)','setprecision(17)').replace('setprecision(10)','setprecision(17)')
s=s.replace('model.set_equilibrium();','auto eq_start=std::chrono::steady_clock::now();\n        model.set_equilibrium();\n        auto eq_end=std::chrono::steady_clock::now();')
s=s.replace('for(int it=0 ; it<Nt ; it++){','''{
                std::ofstream initial("initial_state.tsv"); initial<<std::setprecision(17);
                for(int ik=0;ik<model.Nsys;ik++){
                    auto &r=model.rho_el_k_tstp[ik];
                    initial<<2*real(r(0,1))<<" "<<-2*imag(r(0,1))<<" "<<real(r(0,0)-r(1,1))<<"\\n";
                }
            }
            auto evolution_start=std::chrono::steady_clock::now();
            for(int it=0 ; it<Nt ; it++){''',1)
s=s.replace(f'model.{observable}(Nt);','''model.OBSERVABLE(Nt);
            auto evolution_end=std::chrono::steady_clock::now();
            std::ofstream timing("kernel_timing.tsv");timing<<std::setprecision(17);
            timing<<"equilibrium_seconds "<<std::chrono::duration<double>(eq_end-eq_start).count()<<"\\n";
            timing<<"evolution_seconds "<<std::chrono::duration<double>(evolution_end-evolution_start).count()<<"\\n";''',1)
s=s.replace('model.OBSERVABLE(Nt);',f'model.{observable}(Nt);')
ref=ROOT/'reference';ref.mkdir(exist_ok=True)
main=ref/('SC_dynamics_main_fixed_instrumented.cpp' if args.fixed else 'SC_dynamics_main_instrumented.cpp');main.write_text(s)
files=[main,folder/'SC_model_basic.cpp',LIB/'Lib_utils/basic_const.cpp',LIB/'Lib_utils/Interpolation.cpp',LIB/'Lib_QM/basic_functions.cpp',LIB/'Lib_utils/input_param.cpp']
exe=ref/('sc_model_fixed_reference.ex' if args.fixed else 'sc_model_reference.ex')
cmd=['clang++','-O3','-std=c++11','-Wno-deprecated-declarations','-I/opt/local/include/eigen3',f'-I{folder}',*[str(f) for f in files],'-o',str(exe)]
with open(ref/('build_fixed.log' if args.fixed else 'build.log'),'w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True)
metadata={'fixed_mean_field':args.fixed,'command':cmd,'compiler':subprocess.check_output(['clang++','--version'],text=True),'changes':'main output precision17, equilibrium state dump and steady_clock timing only','sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [src,*files,folder/'SC_model_basic.hpp']}}
(ref/('build_fixed.json' if args.fixed else 'build.json')).write_text(json.dumps(metadata,indent=2)+'\n')
print(exe)
