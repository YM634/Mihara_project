using TOML, Dates
include(joinpath(@__DIR__, "..", "src", "AttractiveHubbard.jl"))
using .AttractiveHubbard
length(ARGS)==4 || error("scan_worker.jl config output worker_id worker_count")
raw=TOML.parsefile(ARGS[1]);out=ARGS[2];worker=parse(Int,ARGS[3]);workers=parse(Int,ARGS[4])
Threads.nthreads()==1 || error("Launch each worker with JULIA_NUM_THREADS=1")
model=raw["model"];scan=raw["scan"]
config(extra)=Config(Float64;(Symbol(k)=>v for (k,v) in merge(model,extra))...)
# Compile before timing; equilibrium is identical for all fields.
simulate(Config(Float64;L=4,tmin="-0.02",tmax="0.02",parallel_k=false))
c0=config(Dict("kind"=>"zero"));eq=equilibrium(c0)
mkpath(out)
function atomic_toml(path,data)
    open(path*".tmp","w") do io;TOML.print(io,data);end
    mv(path*".tmp",path;force=true)
end
function run_case(label,extra)
    dir=joinpath(out,label);mkpath(dir)
    isfile(joinpath(dir,"complete.toml")) && return
    c=config(extra);started=time();r=simulate(c;initial=eq.state)
    all(isfinite,r.observables) && all(isfinite,r.state) || error("Non-finite output: $label")
    for (name,data) in [("observables.f64",r.observables),("initial_state.f64",eq.state),("final_state.f64",r.state)]
        open(joinpath(dir,name*".tmp"),"w") do io;write(io,data);end
        mv(joinpath(dir,name*".tmp"),joinpath(dir,name);force=true)
    end
    atomic_toml(joinpath(dir,"config.toml"),merge(model,extra))
    atomic_toml(joinpath(dir,"complete.toml"),Dict("rows"=>size(r.observables,1),"columns"=>10,"state_shape"=>[3,c.L^2],"dtype"=>"Float64","endianness"=>"little","layout"=>"column-major","column_names"=>["time","n_up","n_down","Phi_re","Phi_im","Jx","Jy","Ekin","Eint","Etot"],"evolution_seconds"=>r.evolution_seconds,"total_seconds"=>time()-started,"finished_at"=>string(now()),"worker"=>worker,"threads"=>Threads.nthreads()))
    println(now()," completed ",label," evolution=",r.evolution_seconds," total=",time()-started);flush(stdout)
end
if worker==1
    for amp in scan["amplitudes"]
        run_case("common_F$(amp)_f",Dict("F0"=>amp,"kind"=>"f","omega"=>"0"))
    end
end
omegas=range(scan["omega_min"],scan["omega_max"];length=scan["omega_count"])
for i in worker:workers:length(omegas)
    for amp in scan["amplitudes"],kind in ["cos_plus","cos_minus","sin_plus","sin_minus"]
        run_case("omega_$(lpad(i-1,3,'0'))_F$(amp)_$(kind)",Dict("F0"=>amp,"kind"=>kind,"omega"=>string(omegas[i])))
    end
end
println(now()," WORKER COMPLETE ",worker);flush(stdout)
