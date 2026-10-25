# Run once in a fresh Julia process for each JULIA_NUM_THREADS setting.
using TOML,Serialization,Statistics
using DoubleFloats
include(joinpath(@__DIR__,"..","src","AttractiveHubbard.jl"))
using .AttractiveHubbard
ROOT=normpath(joinpath(@__DIR__,".."))
OUT=joinpath(ROOT,"results","thread_scaling","threads$(Threads.nthreads())")
ispath(OUT) && error("Output exists: $OUT")
mkpath(OUT)
cases=[("Float64_L50",Float64,50),("Float64_L200",Float64,200),("Double64_L50",Double64,50)]
rows=Any[]
for (label,T,L) in cases, update in (true,false)
    name=label*(update ? "_update" : "_fixed")
    c=Config(T;L=L,dt="0.01",tmin="-7",tmax="7",F0="0.3",epsilon="0.01",omega="1",tau_f="0.1",tau_g="1",polarization=[cosd(22),sind(22)],kind="cos_plus",update_mean_field=update)
    # Compile same precision and MF mode on a short grid before timing.
    simulate(Config(T;L=8,tmin="-0.02",tmax="0.02",update_mean_field=update))
    f=prepare_fields(c);g=grid(c);eq=equilibrium(c,g)
    seconds=Float64[];eqseconds=Float64[];r=nothing
    for trial in 1:3
        GC.gc();r=simulate(c;fields=f)
        push!(seconds,r.evolution_seconds);push!(eqseconds,r.equilibrium_seconds)
    end
    # Serialization preserves Double64, including its low component.
    serialize(joinpath(OUT,name*".jls"),(observables=r.observables,state=r.state))
    row=Dict("case"=>name,"precision"=>label,"L"=>L,"update_mean_field"=>update,"evolution_seconds"=>seconds,"median_seconds"=>median(seconds),"equilibrium_seconds"=>eqseconds)
    push!(rows,row)
    println("threads=",Threads.nthreads()," ",name," median=",median(seconds)," s");flush(stdout)
end
open(joinpath(OUT,"timing.toml"),"w") do io
    TOML.print(io,Dict("julia_threads"=>Threads.nthreads(),"hardware_threads"=>Sys.CPU_THREADS,"cpu"=>Sys.CPU_NAME,"julia_version"=>string(VERSION),"steps"=>1400,"angle_degrees"=>22,"trials"=>3,"cases"=>rows))
end
