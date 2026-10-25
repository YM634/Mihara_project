using TOML, Dates
include(joinpath(@__DIR__,"..","src","RiceMele.jl"))
using .RiceMele
length(ARGS)>=2 || error("Usage: benchmark.jl config.toml Float64|Double64 [output-directory]")
raw=TOML.parsefile(ARGS[1])
mode=ARGS[2]
if mode=="Double64"
    @eval using DoubleFloats
    T=DoubleFloats.Double64
elseif mode=="Float64"
    T=Float64
elseif mode=="BigFloat"
    setprecision(BigFloat,106);T=BigFloat
else
    error("Unknown precision")
end
out=length(ARGS)>2 ? ARGS[3] : joinpath(@__DIR__,"..","results",mode)
ispath(out) && error("Output exists: $out")
mkpath(out)
c=Config(T;(Symbol(k)=>v for (k,v) in raw if k!="amplitudes")...)
amplitudes=parse.(T,raw["amplitudes"])
length(amplitudes)==2 || error("Expected two amplitudes")
# Compile the same typed kernels using a tiny run, exclude this from steady timing.
w=Config(T;nk=8,tmin="-0.02",tmax="0.02",dt="0.01",taug="3")
warm=@elapsed simulate(w,amplitudes[1],:cos)
rows=Any[]
traces=Dict()
for (label,amp) in zip(["low","high"],amplitudes),kind in [:f,:cos,:sin]
    GC.gc()
    measurement=@timed simulate(c,amp,kind)
    trace=measurement.value
    all(isfinite,trace.current) || error("Non-finite current")
    traces[(label,kind)]=trace
    path=joinpath(out,"$(label)_$(kind).tsv")
    write_seconds=@elapsed open(path,"w") do io
        println(io,"# time\tA\tcurrent")
        for i in eachindex(trace.time)
            # string(Double64) retains its high precision, never cast to Float64.
            println(io,trace.time[i],'\t',trace.A[i],'\t',trace.current[i])
        end
    end
    row=Dict("amplitude"=>string(amp),"field"=>string(kind),"label"=>label,"evolution_seconds"=>measurement.time,"write_seconds"=>write_seconds,"allocated_bytes"=>measurement.bytes,"gc_seconds"=>measurement.gctime)
    push!(rows,row)
    println(label," ",kind,": ",round(measurement.time,digits=3)," s (write ",round(write_seconds,digits=3)," s)");flush(stdout)
end
ratio=amplitudes[2]/amplitudes[1]
dc=second_variation(traces[("low",:f)].current,traces[("low",:cos)].current,traces[("high",:f)].current,traces[("high",:cos)].current,ratio)
ds=second_variation(traces[("low",:f)].current,traces[("low",:sin)].current,traces[("high",:f)].current,traces[("high",:sin)].current,ratio)
open(joinpath(out,"variations.tsv"),"w") do io
    println(io,"# time\tdelta_cos_quadratic\tdelta_sin_quadratic")
    for i in eachindex(dc)
        println(io,traces[("high",:f)].time[i],'\t',dc[i],'\t',ds[i])
    end
end
n=round(Int,(c.tmax-c.tmin)/c.dt)
summary=Dict("precision"=>mode,"julia_version"=>string(VERSION),"threads"=>Threads.nthreads(),"cpu"=>Sys.CPU_NAME,"timestamp"=>string(now()),"config"=>raw,"steps"=>n,"runs"=>rows,"warmup_seconds"=>warm,"evolution_total_seconds"=>sum(r["evolution_seconds"] for r in rows),"write_total_seconds"=>sum(r["write_seconds"] for r in rows),"max_abs_delta_cos"=>string(maximum(abs,dc)),"max_abs_delta_sin"=>string(maximum(abs,ds)))
open(joinpath(out,"timing.toml"),"w") do io;TOML.print(io,summary);end
println("Six-trajectory set: ",summary["evolution_total_seconds"]," s evolution")
