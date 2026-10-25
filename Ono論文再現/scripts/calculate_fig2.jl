using TOML, Dates
include(joinpath(@__DIR__,"..","src","RiceMele.jl"))
using .RiceMele
out=joinpath(@__DIR__,"..","results","fig2_float64_stronger")
mkpath(out)
# Both the 2D map and diagonal/inset use actual sampled rows, no frequency interpolation of dynamics.
mapfreq=collect(0.:.04:4.)
diagfreq=sort(unique(round.(vcat(collect(0.:.02:3.),collect(.48:.002:.54),collect(.96:.005:1.08),collect(1.98:.005:2.04));digits=6)))
ws=sort(unique(vcat(mapfreq,diagfreq)))
if !isempty(ARGS); ws=parse.(Float64,ARGS);end
base=Dict("nk"=>800,"dt"=>"0.01","tmin"=>"-1000","tmax"=>"1000","hx"=>"1","hy"=>"0.5","hz"=>"0.1","gamma"=>"0.01","tauf"=>"0.1","epsilon"=>"1e-3")
open(joinpath(out,"parameters.toml"),"w") do io
    TOML.print(io,Dict("config"=>base,"amplitudes"=>["1e-3","1e-2"],"precision"=>"Float64","difference"=>"central epsilon; two-amplitude linear cancellation","map_positive_frequencies"=>mapfreq,"diagonal_frequencies"=>diagfreq,"julia_threads"=>Threads.nthreads(),"julia_version"=>string(VERSION),"cpu"=>Sys.CPU_NAME))
end
makeconfig(w)=Config(Float64;(Symbol(k)=>v for (k,v) in base)...,omega1=w)
# Warm the kernel before measuring.
simulate(Config(Float64;nk=8,tmin="-0.02",tmax="0.02",taug="3"),.01,:cos)
if !isfile(joinpath(out,"baseline.tsv"))
    c=makeconfig(.5)
    lo=simulate(c,.001,:f);hi=simulate(c,.01,:f);cs=simulate(c,.01,:cos)
    open(joinpath(out,"baseline.tsv"),"w") do io
        println(io,"# time A_high J_low_f J_high_f J_high_cos")
        for i in 1:10:length(hi.time)
            println(io,hi.time[i],' ',hi.A[i],' ',lo.current[i],' ',hi.current[i],' ',cs.current[i])
        end
    end
end
lockobj=ReentrantLock()
start=time()
Threads.@threads :dynamic for j in eachindex(ws)
    w=ws[j];tag=lpad(string(round(Int,w*1e6)),7,'0')
    path=joinpath(out,"row_$(tag).bin")
    isfile(path) && continue
    c=makeconfig(w);tick=time()
    # Central variation includes epsilon, exactly removing the epsilon^2 term at quadratic order.
    function variation(amp,plus,minus)
        p=simulate(c,amp,plus).current
        m=simulate(c,amp,minus).current
        (p.-m)./2
    end
    lc=variation(.001,:cos,:cos_minus);hc=variation(.01,:cos,:cos_minus)
    ls=variation(.001,:sin,:sin_minus);hs=variation(.01,:sin,:sin_minus)
    dc=(hc.-10 .*lc)./.9
    ds=(hs.-10 .*ls)./.9
    open(path*".tmp","w") do io
        write(io,dc);write(io,ds)
    end
    mv(path*".tmp",path)
    seconds=time()-tick
    lock(lockobj) do
        open(joinpath(out,"timing.tsv"),"a") do io;println(io,w,'\t',seconds);end
        println("row ",j,"/",length(ws)," omega1=",w," seconds=",round(seconds,digits=2)," elapsed=",round(time()-start,digits=1));flush(stdout)
    end
end
println("Completed ",length(ws)," rows in ",time()-start," seconds")
