using TOML,DelimitedFiles,Statistics,Dates
include(joinpath(@__DIR__,"..","src","AttractiveHubbard.jl"))
using .AttractiveHubbard
ROOT=normpath(joinpath(@__DIR__,".."))
fixed="--fixed" in ARGS
angle_arg=findfirst(a->startswith(a,"--angle="),ARGS)
angle=angle_arg===nothing ? nothing : parse(Float64,split(ARGS[angle_arg],"=";limit=2)[2])
angle===nothing || isfinite(angle) || error("Angle must be finite")
angle_tag=angle===nothing ? "" : (isinteger(angle) ? string(Int(angle)) : replace(string(angle),"."=>"p","-"=>"m"))
folder=angle===nothing ? (fixed ? "cpp_fixed_comparison" : "cpp_comparison") : "angle$(angle_tag)_"*(fixed ? "fixed" : "self_consistent")
OUT=joinpath(ROOT,"results",folder)
mkpath(OUT)
reference=joinpath(ROOT,"reference",fixed ? "sc_model_fixed_reference.ex" : "sc_model_reference.ex")
isfile(reference) || error("Run scripts/build_reference.py first")
base=TOML.parsefile(joinpath(ROOT,"configs","comparison.toml"))
base["update_mean_field"]=!fixed
if angle!==nothing;base["polarization"]=[cosd(angle),sind(angle)];end
cases=[("zero",Dict("kind"=>"zero")),("gaussian_x",Dict("kind"=>"f")),("ono_diagonal",Dict("kind"=>"cos_plus","polarization"=>[1,1])),("doped_sin",Dict("kind"=>"sin_plus","mu"=>"0.5","polarization"=>[1,.5])),("L20_gaussian",Dict("L"=>20,"kind"=>"f"))]
if angle!==nothing
    cases=[("gaussian",Dict("kind"=>"f")),("cos_plus",Dict("kind"=>"cos_plus")),("sin_plus",Dict("kind"=>"sin_plus"))]
end
warmup=@elapsed simulate(Config(Float64;L=4,tmin="-0.02",tmax="0.02",kind="zero",update_mean_field=!fixed))
function make_input(d,c,f)
    mkpath(d)
    for (name,v) in [("Ax_t.txt",f.ax),("Ay_t.txt",f.ay),("Ex_t.txt",f.ex),("Ey_t.txt",f.ey)]
        open(joinpath(d,name),"w") do io;for x in v;println(io,x);end;end
    end
    pairs=["T"=>c.temperature,"mu"=>c.mu,"U_el"=>c.U,"g_elph"=>0,"w0_ph"=>.3,"SC_symmetry"=>"S","elph_type"=>"S","L"=>c.L,"Nt"=>length(f.time)-1,"dt"=>c.dt,"Ax_t"=>"Ax_t.txt","Ay_t"=>"Ay_t.txt","Ex_t"=>"Ex_t.txt","Ey_t"=>"Ey_t.txt","it_obs"=>0]
    open(joinpath(d,"in_params"),"w") do io;for (k,v) in pairs;println(io,"__",k,"= ",v);end;end
end
function run_cpp(exe,d;origin=0.)
    cached=isfile(joinpath(d,"obs_t.dat"))
    seconds=cached ? (isfile(joinpath(d,"wall_seconds.txt")) ? parse(Float64,read(joinpath(d,"wall_seconds.txt"),String)) : NaN) : NaN
    if !cached
    seconds=@elapsed open(joinpath(d,"log"),"w") do log
        open(joinpath(d,"stderr"),"w") do err
            run(pipeline(Cmd(Cmd([exe,"in_params"]);dir=d),stdout=log,stderr=err))
        end
    end
    write(joinpath(d,"wall_seconds.txt"),string(seconds))
    end
    a=readdlm(joinpath(d,"obs_t.dat"),Float64;skipstart=1)
    j=readdlm(joinpath(d,"Jcurr_t.dat"),Float64;skipstart=1)
    hcat(a[:,1].+origin,a[:,2:5],j[:,2:3],a[:,8:9],a[:,7]),seconds
end
summaries=isfile(joinpath(OUT,"summary.toml")) ? Any[TOML.parsefile(joinpath(OUT,"summary.toml"))["cases"]...] : Any[]
for (name,override) in cases
    raw=merge(base,override)
    saved=findfirst(r->r["name"]==name,summaries)
    if saved!==nothing
        merge(Dict("update_mean_field"=>true),summaries[saved]["config"])==raw || error("Completed case has different parameters: $name")
        continue
    end
    c=Config(Float64;(Symbol(k)=>v for (k,v) in raw)...)
    d=joinpath(OUT,name);mkpath(d)
    if isfile(joinpath(d,"config.toml")); merge(Dict("update_mean_field"=>true),TOML.parsefile(joinpath(d,"config.toml")))==raw || error("Existing case has different parameters");end
    open(joinpath(d,"config.toml"),"w") do io;TOML.print(io,raw);end
    f=sample_fields(c)
    clean(v)=[abs(x)<1e-300 ? 0.0 : x for x in v]
    f=merge(f,(;ax=clean(f.ax),ay=clean(f.ay),ex=clean(f.ex),ey=clean(f.ey)))
    prepared=prepare_fields(c,f)
    jruns=Any[];cruns=Any[]
    result=nothing;reference_obs=nothing;cppdir=""
    for repeat in 1:((fixed || angle!==nothing) ? 1 : 3)
        GC.gc()
        jw=@elapsed result=simulate(c;fields=prepared)
        push!(jruns,Dict("setup_seconds"=>result.setup_seconds,"equilibrium_seconds"=>result.equilibrium_seconds,"evolution_seconds"=>result.evolution_seconds,"warm_total_seconds"=>jw))
        cppdir=joinpath(d,"cpp_$repeat");make_input(cppdir,c,f)
        reference_obs,wall=run_cpp(reference,cppdir;origin=c.tmin)
        kt=Dict(String(row[1])=>Float64(row[2]) for row in eachrow(readdlm(joinpath(cppdir,"kernel_timing.tsv"))))
        kt["wall_seconds"]=wall;push!(cruns,kt)
    end
    energy_shift=0.0
    if fixed
        writedlm(joinpath(d,"cpp_observables_raw.tsv"),reference_obs)
        energy_shift=-c.U*abs2(complex(reference_obs[1,4],reference_obs[1,5]))
        reference_obs[:,9].+=energy_shift
        reference_obs[:,10].+=energy_shift
    end
    write_trace(joinpath(d,"julia.tsv"),result)
    writedlm(joinpath(d,"cpp_observables.tsv"),reference_obs)
    initial=permutedims(readdlm(joinpath(cppdir,"initial_state.tsv"),Float64))
    controlled=simulate(c;initial=initial,fields=prepared)
    write_trace(joinpath(d,"julia_cpp_initial.tsv"),controlled)
    a=result.observables;b=reference_obs
    all(isfinite,a) && all(isfinite,b) || error("Nonfinite output")
    size(a)==size(b) || error("Shape mismatch")
    err=vec(maximum(abs.(a-b);dims=1))
    controlled_err=vec(maximum(abs.(controlled.observables-b);dims=1))
    # Positive real equilibrium phase is used by both solvers, no fitted phase alignment.
    phierr=maximum(abs.(complex.(a[:,4],a[:,5]).-complex.(b[:,4],b[:,5])))
    row=Dict("energy_reference_shift"=>energy_shift,"same_initial_max_Phi_complex_error"=>maximum(abs.(complex.(controlled.observables[:,4],controlled.observables[:,5]).-complex.(b[:,4],b[:,5]))),"name"=>name,"config"=>raw,"julia_runs"=>jruns,"cpp_runs"=>cruns,"julia_equilibrium_phi_real"=>real(result.equilibrium.phi),"cpp_equilibrium_phi_real"=>b[1,4],"julia_equilibrium_residual"=>result.equilibrium.residual,"julia_equilibrium_iterations"=>result.equilibrium.iterations,"max_errors"=>err,"same_initial_max_errors"=>controlled_err,"max_Phi_complex_error"=>phierr,"max_Delta_complex_error"=>abs(c.U)*phierr,"cpp_max_Phi_change"=>maximum(abs.(complex.(b[:,4].-b[1,4],b[:,5].-b[1,5]))),"speedup_dynamics"=>median(r["evolution_seconds"] for r in cruns)/median(r["evolution_seconds"] for r in jruns))
    if name=="gaussian_x" && !fixed
        original="/Users/yutamurakami/Documents/Softwares/Lib_TDDM_ver2/Execute_sc/sc_model.ex"
        od=joinpath(d,"original_executable");make_input(od,c,f)
        orig,wall=run_cpp(original,od;origin=c.tmin)
        row["original_binary_wall_seconds"]=wall
        row["original_vs_instrumented_max_errors"]=vec(maximum(abs.(orig-b);dims=1))
        analytic_config=Config(Float64;(Symbol(k)=>v for (k,v) in merge(raw,Dict("field_mode"=>"analytic")))...)
        analytic=simulate(analytic_config)
        write_trace(joinpath(d,"julia_analytic_field.tsv"),analytic)
        row["analytic_vs_table4_max_Phi_error"]=maximum(abs.(complex.(analytic.observables[:,4],analytic.observables[:,5]).-complex.(a[:,4],a[:,5])))
    end
    # Shared-accuracy acceptance limits; these compare implementations, not mesh convergence.
    phierr<1e-9 || error("OP mismatch: $name $phierr")
    maximum(err[6:10])<1e-8 || error("Current/energy mismatch: $name")
    push!(summaries,row)
    open(joinpath(OUT,"summary.toml"),"w") do io
        TOML.print(io,Dict("polarization_angle_degrees"=>(angle===nothing ? "unspecified" : string(angle)),"fixed_mean_field"=>fixed,"julia_version"=>string(VERSION),"threads"=>Threads.nthreads(),"cpu"=>Sys.CPU_NAME,"warmup_seconds"=>warmup,"field_export_cutoff"=>1e-300,"columns"=>["time","n_up","n_down","Phi_re","Phi_im","Jx","Jy","Ekin","Eint","Etot"],"cases"=>summaries))
    end
    println(name," max Phi difference=",phierr," speedup=",row["speedup_dynamics"]);flush(stdout)
end
println("Saved ",OUT)
