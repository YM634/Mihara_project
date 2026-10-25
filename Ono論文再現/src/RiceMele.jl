module RiceMele
export Config, simulate, field, rhs, second_variation, fourier, chi2

struct Config{T}
    nk::Int
    dt::T
    tmin::T
    tmax::T
    hx::T
    hy::T
    hz::T
    gamma::T
    tauf::T
    taug::T
    epsilon::T
    omega1::T
end
function Config(::Type{T}; nk=800, dt="0.01", tmin="-700", tmax="700", hx="1", hy="0.5", hz="0.1", gamma="0.01", tauf="0.1", taug=nothing, epsilon="1e-6", omega1="0.5") where T
    cv(x)=x isa AbstractString ? parse(T,x) : T(x)
    g=cv(gamma)
    c=Config{T}(nk,cv(dt),cv(tmin),cv(tmax),cv(hx),cv(hy),cv(hz),g,cv(tauf),taug===nothing ? sqrt(T(2))/g : cv(taug),cv(epsilon),cv(omega1))
    nk>1 && c.dt>0 && c.tmax>c.tmin && c.tauf>0 && c.taug>0 && c.gamma>=0 || error("Invalid parameters")
    n=round(Int,(c.tmax-c.tmin)/c.dt)
    abs(T(n)*c.dt-(c.tmax-c.tmin)) < c.dt/T(100000) || error("Noninteger step count")
    c
end

# A(t) is the vector potential, not the electric field. E=-dA/dt.
function field(t, amplitude, kind::Symbol, c::Config{T}) where T
    a=amplitude/(T(2)*T(pi))
    f=a*exp(-t*t/(T(2)*c.tauf*c.tauf))
    kind===:f && return f
    kind===:zero && return zero(T)
    g=a*exp(-t*t/(T(2)*c.taug*c.taug))
    kind===:cos && return f+c.epsilon*g*cos(c.omega1*t)
    kind===:cos_minus && return f-c.epsilon*g*cos(c.omega1*t)
    kind===:sin && return f+c.epsilon*g*sin(c.omega1*t)
    kind===:sin_minus && return f-c.epsilon*g*sin(c.omega1*t)
    error("Unknown field kind: $kind")
end

# rho=(I+r.sigma)/2, dr/dt=2h cross r-Gamma(r-r0), r0=-h/|h|.
@inline function rhs(x,y,z,hx,hy,hz,gamma)
    invE=inv(sqrt(hx*hx+hy*hy+hz*hz))
    (2*(hy*z-hz*y)-gamma*(x+hx*invE),
     2*(hz*x-hx*z)-gamma*(y+hy*invE),
     2*(hx*y-hy*x)-gamma*(z+hz*invE))
end
@inline function hfield(sk,ck,sa,ca,c)
    # q=(k-A)/2. Angle addition avoids Nk transcendental calls per RK stage.
    sq=sk*ca-ck*sa; cq=ck*ca+sk*sa
    (-c.hx*cq,-c.hy*sq,-c.hz,sq,cq)
end
function simulate(c::Config{T}, amplitude::T, kind::Symbol; save_stride=1) where T
    save_stride>=1 || error("save_stride must be positive")
    n=round(Int,(c.tmax-c.tmin)/c.dt)
    # [-pi, pi), no endpoint duplication. No physical spin degeneracy factor.
    ks=[-T(pi)+T(2)*T(pi)*T(i)/T(c.nk) for i in 0:c.nk-1]
    sk=sin.(ks./T(2)); ck=cos.(ks./T(2))
    x=Vector{T}(undef,c.nk);y=similar(x);z=similar(x)
    @inbounds for i in eachindex(x)
        hx=-c.hx*ck[i];hy=-c.hy*sk[i];hz=-c.hz
        E=sqrt(hx*hx+hy*hy+hz*hz)
        x[i]=-hx/E;y[i]=-hy/E;z[i]=-hz/E
    end
    ts=T[]; currents=T[]; fields=T[]
    dt=c.dt; half=dt/T(2); sixth=dt/T(6)
    sa,ca=sincos(field(c.tmin,amplitude,kind,c)/T(2))
    for step in 0:n
        t=c.tmin+T(step)*dt
        if step%save_stride==0 || step==n
            j=zero(T)
            @inbounds for i in eachindex(x)
                _,_,_,sq,cq=hfield(sk[i],ck[i],sa,ca,c)
                j+=x[i]*(c.hx/T(2))*sq-y[i]*(c.hy/T(2))*cq
            end
            push!(ts,t); push!(currents,j/T(c.nk));push!(fields,field(t,amplitude,kind,c))
        end
        step==n && break
        sm,cm=sincos(field(t+half,amplitude,kind,c)/T(2))
        se,ce=sincos(field(t+dt,amplitude,kind,c)/T(2))
        @inbounds for i in eachindex(x)
            h1=hfield(sk[i],ck[i],sa,ca,c)
            hm=hfield(sk[i],ck[i],sm,cm,c)
            he=hfield(sk[i],ck[i],se,ce,c)
            a,b,d=x[i],y[i],z[i]
            k1=rhs(a,b,d,h1[1],h1[2],h1[3],c.gamma)
            k2=rhs(a+half*k1[1],b+half*k1[2],d+half*k1[3],hm[1],hm[2],hm[3],c.gamma)
            k3=rhs(a+half*k2[1],b+half*k2[2],d+half*k2[3],hm[1],hm[2],hm[3],c.gamma)
            k4=rhs(a+dt*k3[1],b+dt*k3[2],d+dt*k3[3],he[1],he[2],he[3],c.gamma)
            x[i]=a+sixth*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
            y[i]=b+sixth*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
            z[i]=d+sixth*(k1[3]+2*k2[3]+2*k3[3]+k4[3])
        end
        sa,ca=se,ce
    end
    (;time=ts,current=currents,A=fields)
end

# At fixed shape, J(F)=aF+bF^2+... . This cancels the linear term and
# corrects the finite-low-amplitude contribution to the quadratic term.
function second_variation(low_f,low_test,high_f,high_test,ratio)
    ((high_test.-high_f).-ratio.*(low_test.-low_f))./(1-inv(ratio))
end
function fourier(t,q,omega)
    length(t)==length(q) && length(t)>1 || error("Invalid time series")
    s=zero(complex(q[1]))
    for i in 1:length(t)-1
        s+=(t[i+1]-t[i])/2*(q[i]*cis(omega*t[i])+q[i+1]*cis(omega*t[i+1]))
    end
    s
end
# Main Eq.(14): returns symmetrized chi_bar, not chi_bar/2.
function chi2(t,dc,ds,omega2,Fhigh,c)
    fhat=Fhigh*c.tauf/sqrt(typeof(Fhigh)(2)*typeof(Fhigh)(pi))*exp(-omega2^2*c.tauf^2/2)
    typeof(Fhigh)(2)*typeof(Fhigh)(pi)*(fourier(t,dc,c.omega1+omega2)-im*fourier(t,ds,c.omega1+omega2))/(c.epsilon*Fhigh*fhat)
end
end
