include(joinpath(@__DIR__,"..","src","RiceMele.jl"))
using .RiceMele
out=joinpath(@__DIR__,"..","results","fig2_float64_stronger")
w=.5
for (tag,lo,hi,eps) in [("half_amplitude",.0005,.005,.001),("half_epsilon",.001,.01,.0005)]
    path=joinpath(out,"check_$(tag).bin")
    isfile(path) && continue
    c=Config(Float64;nk=800,dt="0.01",tmin="-1000",tmax="1000",epsilon=eps,omega1=w)
    variation(amp,p,m)=(simulate(c,amp,p).current.-simulate(c,amp,m).current)./2
    lc=variation(lo,:cos,:cos_minus);hc=variation(hi,:cos,:cos_minus)
    ls=variation(lo,:sin,:sin_minus);hs=variation(hi,:sin,:sin_minus)
    dc=(hc.-10 .*lc)./.9;ds=(hs.-10 .*ls)./.9
    open(path,"w") do io;write(io,dc);write(io,ds);end
    println(tag," done");flush(stdout)
end
