using TOML,Serialization,DoubleFloats,DelimitedFiles
ROOT=normpath(joinpath(@__DIR__,".."))
OUT=joinpath(ROOT,"results","thread_scaling")
base=TOML.parsefile(joinpath(OUT,"threads1","timing.toml"))
rows=Any[]
for nt in (1,2,4,8)
    meta=TOML.parsefile(joinpath(OUT,"threads$nt","timing.toml"))
    for item in meta["cases"]
        name=item["case"]
        a=deserialize(joinpath(OUT,"threads1",name*".jls"));b=deserialize(joinpath(OUT,"threads$nt",name*".jls"))
        phi(x)=complex.(x.observables[:,4],x.observables[:,5])
        e=vec(maximum(abs.(a.observables-b.observables);dims=1))
        stateerror=maximum(abs.(a.state-b.state));phierror=maximum(abs.(phi(a)-phi(b)))
        exact=isequal(a.observables,b.observables) && isequal(a.state,b.state)
        limit=item["precision"]=="Double64_L50" ? parse(Double64,"1e-27") : 1e-12
        maximum(e)<limit && stateerror<limit || error("Thread mismatch $nt $name")
        reference=only(filter(c->c["case"]==name,base["cases"]))
        push!(rows,merge(item,Dict("threads"=>nt,"speedup"=>reference["median_seconds"]/item["median_seconds"],"bitwise_equal"=>exact,"max_OP_error"=>string(phierror),"max_current_error"=>string(maximum(e[6:7])),"max_energy_error"=>string(e[10]),"max_state_error"=>string(stateerror))))
    end
end
open(joinpath(OUT,"comparison.toml"),"w") do io
    TOML.print(io,Dict("cpu"=>base["cpu"],"hardware_threads"=>base["hardware_threads"],"trials"=>3,"cases"=>rows))
end
println("All thread counts agree. Bitwise identical: ",all(r["bitwise_equal"] for r in rows))

# Confirm that the one-thread reference still reproduces the pre-threading 22-degree traces.
for (suffix,oldfolder) in [("update","angle22_self_consistent"),("fixed","angle22_fixed")]
    old=readdlm(joinpath(ROOT,"results",oldfolder,"cos_plus","julia.tsv"),Float64;skipstart=1)
    current=deserialize(joinpath(OUT,"threads1","Float64_L50_$(suffix).jls")).observables
    @assert isequal(old,current)
end
println("Previous single-thread 22-degree traces are bitwise unchanged.")
