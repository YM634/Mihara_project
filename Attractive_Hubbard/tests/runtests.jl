using Test,LinearAlgebra
include(joinpath(@__DIR__,"..","src","AttractiveHubbard.jl"))
using .AttractiveHubbard
const SX=ComplexF64[0 1;1 0];const SY=ComplexF64[0 -im;im 0];const SZ=ComplexF64[1 0;0 -1]
function matrix_evolve(c,initial)
    g=grid(c);f=prepare_fields(c);n=c.L^2
    rho=[(Matrix{ComplexF64}(I,2,2)+initial[1,i]*SX+initial[2,i]*SY+initial[3,i]*SZ)/2 for i in 1:n]
    equilibrium_phi=sum(m[1,2] for m in rho)/n
    function rhs(r,ax,ay)
        phi=c.update_mean_field ? sum(m[1,2] for m in r)/n : equilibrium_phi
        [begin
            epsminus=-2*c.hopping*(g.cx[i]*cos(ax)+g.sx[i]*sin(ax)+g.cy[i]*cos(ay)+g.sy[i]*sin(ay))-c.mu
            epsplus=-2*c.hopping*(g.cx[i]*cos(ax)-g.sx[i]*sin(ax)+g.cy[i]*cos(ay)-g.sy[i]*sin(ay))-c.mu
            h=ComplexF64[epsminus c.U*phi;c.U*conj(phi) -epsplus]
            -im*(h*r[i]-r[i]*h)
        end for i in 1:n]
    end
    dt=c.dt
    for it in 1:length(f.samples.time)-1
        a=rhs(rho,f.ax[2*it-1],f.ay[2*it-1])
        b=rhs([rho[i]+dt*a[i]/2 for i in 1:n],f.ax[2*it],f.ay[2*it])
        d=rhs([rho[i]+dt*b[i]/2 for i in 1:n],f.ax[2*it],f.ay[2*it])
        e=rhs([rho[i]+dt*d[i] for i in 1:n],f.ax[2*it+1],f.ay[2*it+1])
        rho=[rho[i]+dt*(a[i]+2*b[i]+2*d[i]+e[i])/6 for i in 1:n]
    end
    reduce(hcat,[[real(tr(r*SX)),real(tr(r*SY)),real(tr(r*SZ))] for r in rho])
end
@testset "Attractive Hubbard dynamics" begin
    c=Config(Float64;L=6,mu="0.4",tmin="-0.2",tmax="0.2",F0="0.8",tau_f="0.1",tau_g="0.1",polarization=[1,1],kind="sin_plus")
    eq=equilibrium(c)
    @test eq.residual<c.tolerance
    @test real(eq.phi)>0
    r=simulate(c;initial=eq.state)
    @test maximum(abs.(r.state-matrix_evolve(c,eq.state)))<2e-14
    # Scalar-identity part of the Nambu Hamiltonian is retained in the matrix check.
    @test size(r.observables)==(41,10)
    for j in 0:4
        v=[Float64(i)^j for i in 0:10]
        @test AttractiveHubbard.midpoint(v,3)≈3.5^j atol=1e-12
    end
    ca=Config(Float64;L=4,tmin="-0.1",tmax="0.1",field_mode="analytic")
    fa=prepare_fields(ca)
    @test fa.ax[2]≈AttractiveHubbard.scalar_field(ca.tmin+ca.dt/2,ca)[1]
    q=Config(Float64;L=50,kind="zero",tmin="-0.1",tmax="0.1")
    z=simulate(q)
    @test -q.U*real(z.equilibrium.phi)≈1.38130780884 atol=1e-10
    @test maximum(abs.(z.observables[:,4].-z.observables[1,4]))<1e-11
    @test maximum(abs,z.observables[:,6:7])<1e-13
    @test maximum(abs.(z.observables[:,2:3].-.5))<1e-13
    @test_throws Exception Config(Float64;U="4")
    @test_throws Exception Config(Float64;field_mode="unknown")
end

@testset "Frozen equilibrium mean field" begin
    c=Config(Float64;L=6,mu="0.4",tmin="-0.2",tmax="0.2",F0="0.8",tau_f="0.1",tau_g="0.1",polarization=[1,1],kind="sin_plus",update_mean_field=false)
    eq=equilibrium(c);r=simulate(c;initial=eq.state)
    @test !c.update_mean_field
    @test maximum(abs.(r.state-matrix_evolve(c,eq.state)))<2e-14
    # Only the coupling in H is fixed; the measured anomalous expectation evolves.
    @test maximum(abs.(r.observables[:,4].-r.observables[1,4]))>1e-6
    sc=Config(Float64;L=6,mu="0.4",tmin="-0.2",tmax="0.2",F0="0.8",tau_f="0.1",tau_g="0.1",polarization=[1,1],kind="sin_plus",update_mean_field=true)
    @test maximum(abs.(r.state-simulate(sc;initial=eq.state).state))>1e-6
    # Energy is <H_frozen> plus its equilibrium double-counting constant.
    g=grid(c);f=prepare_fields(c);ax=f.ax[end];ay=f.ay[end];E=0.
    for i in 1:c.L^2
        rho=(Matrix{ComplexF64}(I,2,2)+r.state[1,i]*SX+r.state[2,i]*SY+r.state[3,i]*SZ)/2
        em=-2*c.hopping*(g.cx[i]*cos(ax)+g.sx[i]*sin(ax)+g.cy[i]*cos(ay)+g.sy[i]*sin(ay))-c.mu
        ep=-2*c.hopping*(g.cx[i]*cos(ax)-g.sx[i]*sin(ax)+g.cy[i]*cos(ay)-g.sy[i]*sin(ay))-c.mu
        h=ComplexF64[em c.U*eq.phi;c.U*conj(eq.phi) -ep]
        E+=real(tr(h*rho))/c.L^2
    end
    E-=c.U*abs2(eq.phi)
    @test r.observables[end,10]≈E atol=1e-14
    z=simulate(Config(Float64;L=6,kind="zero",tmin="-0.1",tmax="0.1",update_mean_field=false))
    @test maximum(abs.(z.observables[:,4].-z.observables[1,4]))<1e-11
    @test_throws Exception Config(Float64;update_mean_field="false")
end

@testset "Threaded and explicit serial paths" begin
    for update in (true,false)
        cfg(parallel)=Config(Float64;L=12,mu="0.3",tmin="-0.2",tmax="0.2",tau_g="0.1",polarization=[1,0.4],update_mean_field=update,parallel_k=parallel)
        a=simulate(cfg(true));b=simulate(cfg(false))
        @test isequal(a.state,b.state)
        @test isequal(a.observables,b.observables)
    end
    @test_throws Exception Config(Float64;parallel_k="true")
end
