using Test, DoubleFloats
include(joinpath(@__DIR__,"..","src","RiceMele.jl"))
using .RiceMele
@testset "Double64 against 106-bit BigFloat" begin
    setprecision(BigFloat,106)
    config(T)=Config(T;nk=8,dt="0.01",tmin="-0.1",tmax="0.1",taug="0.02",tauf="0.02")
    dd=simulate(config(Double64),parse(Double64,"0.1"),:cos)
    big=simulate(config(BigFloat),parse(BigFloat,"0.1"),:cos)
    @test maximum(abs.(BigFloat.(dd.current).-big.current)) < big"1e-28"
    z=simulate(config(Double64),zero(Double64),:zero)
    @test maximum(abs,z.current)<parse(Double64,"1e-29")
end
