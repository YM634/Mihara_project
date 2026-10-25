using TOML
include(joinpath(@__DIR__,"..","src","AttractiveHubbard.jl"))
using .AttractiveHubbard
length(ARGS)>=2 || error("Usage: run.jl config.toml output-directory [Float64|Double64|BigFloat]")
raw=TOML.parsefile(ARGS[1])
mode=length(ARGS)>2 ? ARGS[3] : get(raw,"precision","Float64")
model=copy(raw);pop!(model,"precision",nothing)
if mode=="Float64"
    T=Float64
elseif mode=="Double64"
    @eval using DoubleFloats
    T=DoubleFloats.Double64
elseif mode=="BigFloat"
    setprecision(BigFloat,106);T=BigFloat
else
    error("Unknown precision")
end
out=ARGS[2];ispath(out) && error("Output already exists")
c=Config(T;(Symbol(k)=>v for (k,v) in model)...)
# JIT warmup: same types and code paths, tiny physical grid.
warmup=@elapsed simulate(Config(T;L=4,tmin="-0.02",tmax="0.02",tau_g="1",update_mean_field=c.update_mean_field,parallel_k=c.parallel_k))
mkpath(out)
r=simulate(c)
write_seconds=@elapsed write_trace(joinpath(out,"observables.tsv"),r)
open(joinpath(out,"timing.toml"),"w") do io
    TOML.print(io,Dict("precision"=>mode,"update_mean_field"=>c.update_mean_field,"hamiltonian_OP"=>(c.update_mean_field ? "self_consistent" : string(r.equilibrium.phi)),"julia_version"=>string(VERSION),"threads"=>Threads.nthreads(),"parallel_k"=>c.parallel_k,"active_k_threads"=>(c.parallel_k ? Threads.nthreads() : 1),"warmup_seconds"=>warmup,"setup_seconds"=>r.setup_seconds,"equilibrium_seconds"=>r.equilibrium_seconds,"evolution_seconds"=>r.evolution_seconds,"write_seconds"=>write_seconds,"equilibrium_iterations"=>r.equilibrium.iterations,"equilibrium_residual"=>string(r.equilibrium.residual),"Phi_equilibrium"=>string(r.equilibrium.phi),"input"=>TOML.parsefile(ARGS[1])))
end
println("Phi=",r.equilibrium.phi," evolution=",r.evolution_seconds," s")
