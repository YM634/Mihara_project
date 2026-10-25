# One fixed test frequency. Deliberately no 2D scan or spectral calculation.
using TOML
include(joinpath(@__DIR__,"..","src","AttractiveHubbard.jl"))
using .AttractiveHubbard
length(ARGS)>=2 || error("Usage: run_set.jl one_frequency.toml output-directory [Float64|Double64]")
raw=TOML.parsefile(ARGS[1]);mode=length(ARGS)>2 ? ARGS[3] : get(raw,"precision","Float64")
if mode=="Float64"
    T=Float64
elseif mode=="Double64"
    @eval using DoubleFloats
    T=DoubleFloats.Double64
else
    error("Unsupported precision")
end
out=ARGS[2];ispath(out) && error("Output exists");mkpath(out)
simulate(Config(T;L=4,tmin="-0.02",tmax="0.02",update_mean_field=get(raw["model"],"update_mean_field",true),parallel_k=get(raw["model"],"parallel_k",true)))
timings=Any[]
for amp in raw["set"]["amplitudes"],kind in raw["set"]["kinds"]
    input=merge(raw["model"],Dict("F0"=>amp,"kind"=>kind))
    c=Config(T;(Symbol(k)=>v for (k,v) in input)...)
    r=simulate(c)
    d=joinpath(out,"F$(amp)_$(kind)");mkpath(d)
    write_trace(joinpath(d,"observables.tsv"),r)
    open(joinpath(d,"config.toml"),"w") do io;TOML.print(io,input);end
    push!(timings,Dict("F0"=>amp,"kind"=>kind,"equilibrium_seconds"=>r.equilibrium_seconds,"evolution_seconds"=>r.evolution_seconds))
    println(amp," ",kind," ",r.evolution_seconds," s");flush(stdout)
end
open(joinpath(out,"timing.toml"),"w") do io;TOML.print(io,Dict("precision"=>mode,"threads"=>Threads.nthreads(),"parallel_k"=>get(raw["model"],"parallel_k",true),"runs"=>timings));end
