using Test, LinearAlgebra
include(joinpath(@__DIR__,"..","src","RiceMele.jl"))
using .RiceMele
const sx=ComplexF64[0 1;1 0]
const sy=ComplexF64[0 -im;im 0]
const sz=ComplexF64[1 0;0 -1]
H(k,A,c)=-c.hx*cos((k-A)/2)*sx-c.hy*sin((k-A)/2)*sy-c.hz*sz
function reference_current(c,amp,kind)
    ks=[-pi+2pi*i/c.nk for i in 0:c.nk-1]
    rho=[(Matrix{ComplexF64}(I,2,2)-H(k,0,c)/sqrt(real(tr(H(k,0,c)^2))/2))/2 for k in ks]
    function dr(r,k,t)
        h=H(k,field(t,amp,kind,c),c)
        eq=(Matrix{ComplexF64}(I,2,2)-h/sqrt(real(tr(h*h))/2))/2
        -im*(h*r-r*h)-c.gamma*(r-eq)
    end
    for it in 0:round(Int,(c.tmax-c.tmin)/c.dt)-1
        t=c.tmin+it*c.dt;d=c.dt
        for i in eachindex(ks)
            r=rho[i];k=ks[i]
            a=dr(r,k,t);b=dr(r+d*a/2,k,t+d/2);e=dr(r+d*b/2,k,t+d/2);f=dr(r+d*e,k,t+d)
            rho[i]=r+d*(a+2*b+2*e+f)/6
        end
    end
    A=field(c.tmax,amp,kind,c)
    sum(real(tr(r*((c.hx/2)*sin((k-A)/2)*sx-(c.hy/2)*cos((k-A)/2)*sy))) for (r,k) in zip(rho,ks))/c.nk
end
@testset "Dynamics, units and order separation" begin
    c=Config(Float64;nk=16,dt="0.01",tmin="-1",tmax="1",taug="0.1",epsilon="0.2")
    z=simulate(c,0.,:zero)
    @test maximum(abs,z.current)<1e-14
    @test length(z.time)==201
    @test field(0.,0.2,:f,c)≈0.2/(2pi)
    @test field(0.,0.2,:sin,c)≈field(0.,0.2,:f,c)
    @test field(0.,0.2,:cos,c)≈1.2*field(0.,0.2,:f,c)
    for kind in [:f,:cos,:sin]
        s=simulate(c,0.2,kind)
        @test s.current[end]≈reference_current(c,0.2,kind) atol=5e-15
    end
    hx,hy,hz=-.7,-.3,-.1
    r=[.2,.1,-.3];rho=(I+r[1]*sx+r[2]*sy+r[3]*sz)/2
    h=hx*sx+hy*sy+hz*sz; eq=(I-h/sqrt(hx^2+hy^2+hz^2))/2
    dr=-im*(h*rho-rho*h)-c.gamma*(rho-eq)
    @test collect(rhs(r...,hx,hy,hz,c.gamma))≈real.([tr(dr*sx),tr(dr*sy),tr(dr*sz)])
    # Synthetic response includes linear and quadratic pieces, no dynamics assumptions.
    J(F,g)=3F*(1+g)+7F^2*(1+g)^2
    lo=.01;hi=.1;e=.001
    result=second_variation([J(lo,0)],[J(lo,e)],[J(hi,0)],[J(hi,e)],hi/lo)
    @test result[1]≈7hi^2*((1+e)^2-1) rtol=1e-10
    # RK4 convergence against finer solution in an appreciably driven transient.
    cf(dt)=Config(Float64;nk=16,dt=dt,tmin="-1",tmax="1",taug="0.1")
    js=[simulate(cf(d),1.,:f).current[end] for d in ["0.04","0.02","0.01","0.005"]]
    @test abs(js[1]-js[4])>10abs(js[2]-js[4])
    @test abs(js[2]-js[4])>10abs(js[3]-js[4])
end
