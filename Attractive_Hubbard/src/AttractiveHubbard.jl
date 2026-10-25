module AttractiveHubbard
export Config, grid, equilibrium, sample_fields, prepare_fields, simulate, meanphi, rhs!, write_trace

struct Config{T}
    L::Int
    temperature::T
    U::T
    mu::T
    hopping::T
    dt::T
    tmin::T
    tmax::T
    tolerance::T
    max_iterations::Int
    F0::T
    epsilon::T
    omega::T
    tau_f::T
    tau_g::T
    polarization::NTuple{2,T}
    kind::Symbol
    field_mode::Symbol
    update_mean_field::Bool
    parallel_k::Bool
end
function Config(::Type{T}; L=50,temperature="0.02",U="-4",mu="0",hopping="1",dt="0.01",tmin="-7",tmax="7",tolerance="1e-12",max_iterations=100000,F0="0.3",epsilon="0.01",omega="1",tau_f="0.1",tau_g="1",polarization=[1,0],kind="cos_plus",field_mode="table4",update_mean_field=true,parallel_k=true) where T
    parallel_k isa Bool || error("parallel_k must be a Boolean")
    update_mean_field isa Bool || error("update_mean_field must be a Boolean")
    cv(x)=x isa AbstractString ? parse(T,x) : T(x)
    p=cv.(polarization);length(p)==2 || error("polarization needs 2 components")
    norm=sqrt(p[1]^2+p[2]^2);norm>0 || error("Zero polarization")
    c=Config{T}(L,cv(temperature),cv(U),cv(mu),cv(hopping),cv(dt),cv(tmin),cv(tmax),cv(tolerance),max_iterations,cv(F0),cv(epsilon),cv(omega),cv(tau_f),cv(tau_g),(p[1]/norm,p[2]/norm),Symbol(kind),Symbol(field_mode),update_mean_field,parallel_k)
    L>=2 && c.temperature>0 && c.U<0 && c.hopping>0 && c.dt>0 && c.tmax>c.tmin && c.tolerance>0 && max_iterations>0 && c.tau_f>0 && c.tau_g>0 || error("Invalid model parameters")
    c.kind in (:zero,:f,:cos_plus,:cos_minus,:sin_plus,:sin_minus) || error("Unknown field")
    c.field_mode in (:table4,:analytic) || error("Unknown field mode")
    nt=round(Int,(c.tmax-c.tmin)/c.dt)
    nt>=4 && abs(nt*c.dt-(c.tmax-c.tmin))<c.dt/100000 || error("Invalid time grid")
    c
end
function grid(c::Config{T}) where T
    angles=[T(i)*T(2)*T(pi)/T(c.L) for i in 0:c.L-1]
    # Same ordering as C++: ik=ikx+L*iky.
    cx=[cos(angles[ix]) for iy in 1:c.L for ix in 1:c.L]
    sx=[sin(angles[ix]) for iy in 1:c.L for ix in 1:c.L]
    cy=[cos(angles[iy]) for iy in 1:c.L for ix in 1:c.L]
    sy=[sin(angles[iy]) for iy in 1:c.L for ix in 1:c.L]
    (;cx,sx,cy,sy)
end
function meanphi(r::AbstractMatrix{T}) where T
    x=zero(T);y=zero(T)
    @inbounds for i in axes(r,2);x+=r[1,i];y+=r[2,i];end
    complex(x,-y)/(T(2)*T(size(r,2)))
end
function equilibrium(c::Config{T},g=grid(c)) where T
    n=c.L^2;r=zeros(T,3,n);r[1,:].=T(7)/T(10)
    next=similar(r);residual=T(Inf)
    for iteration in 1:c.max_iterations
        phi=meanphi(r);hx=c.U*real(phi);hy=-c.U*imag(phi);residual=zero(T)
        @inbounds for i in 1:n
            hz=-T(2)*c.hopping*(g.cx[i]+g.cy[i])-c.mu
            E=sqrt(hx^2+hy^2+hz^2)
            # C++ uses Fermi functions clipped at beta*E>50.
            weight=E==0 ? inv(T(2)*c.temperature) : (E/c.temperature>50 ? one(T) : tanh(E/(T(2)*c.temperature)))/E
            a=-hx*weight;b=-hy*weight;d=-hz*weight
            next[1,i]=a;next[2,i]=b;next[3,i]=d
            residual+=sqrt(((a-r[1,i])^2+(b-r[2,i])^2+(d-r[3,i])^2)/T(2))
        end
        residual/=T(n);r,next=next,r
        residual<c.tolerance && return (;state=r,phi=meanphi(r),iterations=iteration,residual)
    end
    error("Equilibrium did not converge: residual=$residual")
end
function scalar_field(t,c::Config{T}) where T
    c.kind===:zero && return (zero(T),zero(T))
    f=c.F0/(T(2)*T(pi))*exp(-t*t/(T(2)*c.tau_f^2))
    df=-t/c.tau_f^2*f
    c.kind===:f && return (f,-df)
    env=c.F0/(T(2)*T(pi))*exp(-t*t/(T(2)*c.tau_g^2))
    denv=-t/c.tau_g^2*env
    if c.kind in (:cos_plus,:cos_minus)
        g=env*cos(c.omega*t);dg=denv*cos(c.omega*t)-c.omega*env*sin(c.omega*t)
    else
        g=env*sin(c.omega*t);dg=denv*sin(c.omega*t)+c.omega*env*cos(c.omega*t)
    end
    s=c.kind in (:cos_minus,:sin_minus) ? -one(T) : one(T)
    (f+s*c.epsilon*g,-df-s*c.epsilon*dg)
end
function sample_fields(c::Config{T}) where T
    nt=round(Int,(c.tmax-c.tmin)/c.dt)
    times=[c.tmin+T(i)*c.dt for i in 0:nt]
    values=[scalar_field(t,c) for t in times]
    (;time=times,ax=[c.polarization[1]*v[1] for v in values],ay=[c.polarization[2]*v[1] for v in values],ex=[c.polarization[1]*v[2] for v in values],ey=[c.polarization[2]*v[2] for v in values])
end
# Five-node fourth-degree interpolation of the same table supplied to C++.
function midpoint(values::Vector{T},i::Int) where T
    n=length(values)-1;start=clamp(i-2,0,n-4);z=T(i-start)+one(T)/T(2)
    result=zero(T)
    for j in 0:4
        weight=one(T)
        for k in 0:4
            j==k && continue
            weight*=(z-T(k))/T(j-k)
        end
        result+=weight*values[start+j+1]
    end
    result
end
function prepare_fields(c::Config{T},f=sample_fields(c)) where T
    n=length(f.time)-1;ax=zeros(T,2*n+1);ay=similar(ax)
    for i in 0:n
        ax[2*i+1]=f.ax[i+1];ay[2*i+1]=f.ay[i+1]
        i==n && continue
        if c.field_mode===:table4
            ax[2*i+2]=midpoint(f.ax,i);ay[2*i+2]=midpoint(f.ay,i)
        else
            a,_=scalar_field(f.time[i+1]+c.dt/T(2),c)
            ax[2*i+2]=c.polarization[1]*a;ay[2*i+2]=c.polarization[2]*a
        end
    end
    (;ax,ay,samples=f)
end
function rhs!(dr,r,c,g,ax,ay,fixed_phi=nothing)
    phi=fixed_phi===nothing ? meanphi(r) : fixed_phi;hx=c.U*real(phi);hy=-c.U*imag(phi)
    ca=cos(ax);cb=cos(ay)
    @inbounds for i in axes(r,2)
        hz=-2*c.hopping*(g.cx[i]*ca+g.cy[i]*cb)-c.mu
        x,y,z=r[1,i],r[2,i],r[3,i]
        dr[1,i]=2*(hy*z-hz*y)
        dr[2,i]=2*(hz*x-hx*z)
        dr[3,i]=2*(hx*y-hy*x)
    end
    dr
end
function observe(r,c,g,ax,ay,fixed_phi=nothing)
    T=eltype(r);sx,cx=sincos(ax);sy,cy=sincos(ay)
    phi=meanphi(r);rz=zero(T);jx=zero(T);jy=zero(T);ekin=zero(T)
    @inbounds for i in axes(r,2)
        rz+=r[3,i]
        jx+=2*c.hopping*(g.sx[i]*cx-r[3,i]*g.cx[i]*sx)
        jy+=2*c.hopping*(g.sy[i]*cy-r[3,i]*g.cy[i]*sy)
        d0=-2*c.hopping*(g.sx[i]*sx+g.sy[i]*sy)
        hz=-2*c.hopping*(g.cx[i]*cx+g.cy[i]*cy)-c.mu
        ekin+=d0+hz*r[3,i]
    end
    n=T(size(r,2));density=(1+rz/n)/2;ekin/=n
    # Frozen Hamiltonian energy, including the constant that matches equilibrium.
    eint=fixed_phi===nothing ? c.U*abs2(phi) : 2*c.U*real(conj(fixed_phi)*phi)-c.U*abs2(fixed_phi)
    (density,density,real(phi),imag(phi),jx/n,jy/n,ekin,eint,ekin+eint)
end
# Keep the global OP reduction in the original summation order. Only disjoint
# k-point writes run in parallel, with a barrier between every RK stage.
function threaded_stage!(k,tmp,rstage,rbase,c,g,ax,ay,coefficient,fixed_phi)
    phi=fixed_phi===nothing ? meanphi(rstage) : fixed_phi
    hx=c.U*real(phi);hy=-c.U*imag(phi);ca=cos(ax);cb=cos(ay)
    n=size(rstage,2);chunks=min(n,Threads.nthreads())
    Threads.@threads :dynamic for chunk in 1:chunks
        first=fld((chunk-1)*n,chunks)+1;last=fld(chunk*n,chunks)
        @inbounds for i in first:last
            hz=-2*c.hopping*(g.cx[i]*ca+g.cy[i]*cb)-c.mu
            x,y,z=rstage[1,i],rstage[2,i],rstage[3,i]
            a=2*(hy*z-hz*y);b=2*(hz*x-hx*z);d=2*(hx*y-hy*x)
            k[1,i]=a;k[2,i]=b;k[3,i]=d
            tmp[1,i]=rbase[1,i]+coefficient*a
            tmp[2,i]=rbase[2,i]+coefficient*b
            tmp[3,i]=rbase[3,i]+coefficient*d
        end
    end
    nothing
end
function threaded_final!(r,tmp,k1,k2,k3,c,g,ax,ay,sixth,fixed_phi)
    phi=fixed_phi===nothing ? meanphi(tmp) : fixed_phi
    hx=c.U*real(phi);hy=-c.U*imag(phi);ca=cos(ax);cb=cos(ay)
    n=size(r,2);chunks=min(n,Threads.nthreads())
    Threads.@threads :dynamic for chunk in 1:chunks
        first=fld((chunk-1)*n,chunks)+1;last=fld(chunk*n,chunks)
        @inbounds for i in first:last
            hz=-2*c.hopping*(g.cx[i]*ca+g.cy[i]*cb)-c.mu
            x,y,z=tmp[1,i],tmp[2,i],tmp[3,i]
            a=2*(hy*z-hz*y);b=2*(hz*x-hx*z);d=2*(hx*y-hy*x)
            r[1,i]=r[1,i]+sixth*(k1[1,i]+2*k2[1,i]+2*k3[1,i]+a)
            r[2,i]=r[2,i]+sixth*(k1[2,i]+2*k2[2,i]+2*k3[2,i]+b)
            r[3,i]=r[3,i]+sixth*(k1[3,i]+2*k2[3,i]+2*k3[3,i]+d)
        end
    end
    nothing
end

function simulate(c::Config{T}; initial=nothing, fields=nothing) where T
    init_start=time_ns();g=grid(c);f=fields===nothing ? prepare_fields(c) : fields
    setup_seconds=(time_ns()-init_start)/1e9
    eq_start=time_ns()
    eq=initial===nothing ? equilibrium(c,g) : (;state=copy(initial),phi=meanphi(initial),iterations=0,residual=zero(T))
    equilibrium_seconds=(time_ns()-eq_start)/1e9
    fixed_phi=c.update_mean_field ? nothing : eq.phi
    r=copy(eq.state);size(r)==(3,c.L^2) || error("Bad initial state shape")
    k1=similar(r);k2=similar(r);k3=similar(r);k4=similar(r);tmp=similar(r)
    n=length(f.samples.time);obs=Matrix{T}(undef,n,10)
    dt=c.dt;half=dt/T(2);sixth=dt/T(6)
    evolution_start=time_ns()
    for it in 1:n
        ax=f.ax[2*it-1];ay=f.ay[2*it-1]
        o=observe(r,c,g,ax,ay,fixed_phi);obs[it,1]=f.samples.time[it]
        @inbounds for j in 1:9;obs[it,j+1]=o[j];end
        it==n && break
        if c.parallel_k && Threads.nthreads()>1
            threaded_stage!(k1,tmp,r,r,c,g,ax,ay,half,fixed_phi)
            threaded_stage!(k2,tmp,tmp,r,c,g,f.ax[2*it],f.ay[2*it],half,fixed_phi)
            threaded_stage!(k3,tmp,tmp,r,c,g,f.ax[2*it],f.ay[2*it],dt,fixed_phi)
            threaded_final!(r,tmp,k1,k2,k3,c,g,f.ax[2*it+1],f.ay[2*it+1],sixth,fixed_phi)
        else
        rhs!(k1,r,c,g,ax,ay,fixed_phi)
        @. tmp=r+half*k1
        rhs!(k2,tmp,c,g,f.ax[2*it],f.ay[2*it],fixed_phi)
        @. tmp=r+half*k2
        rhs!(k3,tmp,c,g,f.ax[2*it],f.ay[2*it],fixed_phi)
        @. tmp=r+dt*k3
        rhs!(k4,tmp,c,g,f.ax[2*it+1],f.ay[2*it+1],fixed_phi)
        @. r=r+sixth*(k1+2*k2+2*k3+k4)
        end
    end
    evolution_seconds=(time_ns()-evolution_start)/1e9
    (;observables=obs,state=r,equilibrium=eq,fields=f.samples,setup_seconds,equilibrium_seconds,evolution_seconds)
end
function write_trace(path,result)
    open(path,"w") do io
        println(io,"# time n_up n_down Phi_re Phi_im Jx Jy Ekin Eint Etot")
        for i in axes(result.observables,1)
            println(io,join(string.(result.observables[i,:]),'\t'))
        end
    end
end
end
