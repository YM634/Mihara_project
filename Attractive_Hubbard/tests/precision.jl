using Test,DoubleFloats
include(joinpath(@__DIR__,"..","src","AttractiveHubbard.jl"))
using .AttractiveHubbard
@testset "Double64 and 106-bit BigFloat" begin
    setprecision(BigFloat,106)
    config(T,update)=Config(T;update_mean_field=update,L=4,mu="0.3",tmin="-0.1",tmax="0.1",tolerance="1e-24",tau_g="0.05",tau_f="0.05",F0="0.5")
    for update in (true,false)
        d=simulate(config(Double64,update));b=simulate(config(BigFloat,update))
        @test eltype(d.observables)===Double64
        @test maximum(abs.(BigFloat.(d.state).-b.state))<big"1e-26"
        @test abs(BigFloat(real(d.equilibrium.phi))-real(b.equilibrium.phi))<big"1e-26"
        @test d.equilibrium.residual<parse(Double64,"1e-24")
    end
end
