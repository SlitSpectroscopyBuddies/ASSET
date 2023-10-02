using EasyFITS
using PyPlot
using LinearAlgebra
using LazyAlgebra
using DelimitedFiles
using InterpolationKernels
using LinearInterpolators
using Statistics
using PointSpreadFunctions
using OptimPackNextGen
import OptimPackNextGen.Powell.Bobyqa
import OptimPackNextGen.Brent
using SparseArrays
using Dates

struct jwstPSF <: AbstractPSF{1}
    α::Float64
    β::Float64
    jwstPSF(p::AbstractArray{Float64,1}) = new(p[1], p[2])
end
    

function (P::jwstPSF)(r::T,λ::T) where {T<:AbstractFloat}
    σ² = P.α*λ*λ + P.β
    return exp(-1/2 *r*r/σ²)/sqrt(2*pi*σ²)
end

struct jwstpolyPSF <: AbstractPSF{1}
    a::Float64
    b::Float64
    c::Float64
    jwstpolyPSF(p::AbstractArray{Float64,1}) = new(p[1], p[2], p[3])
end
    

function (P::jwstpolyPSF)(r::T,λ::T) where {T<:AbstractFloat}
    b = max(P.b, - (P.a*λ*λ + P.c)/λ+1e-8)
    σ² = P.a*λ*λ + b*λ + P.c
    return exp(-1/2 *r*r/σ²)/sqrt(2*pi*σ²)
end

function psf_map!(map::AbstractArray{T,N},
                 h::AbstractPSF,
                 ρ::AbstractArray{T,N},
                 λ::AbstractArray{T,N}) where {N, T<:AbstractFloat}
    sz=size(map)
    @assert sz == size(ρ)
    @assert sz == size(λ)
    
    if ndims(map) == 3
        indices_list = CartesianIndices(map[1:end,1:end,1:end])
    else
        indices_list = CartesianIndices(map[1:end,1:end])
    end
    for ind in indices_list
        map[ind] = h(ρ[ind], λ[ind])
    end
end
      
function psf_map(h::AbstractPSF,
                 ρ::AbstractArray{T,N},
                 λ::AbstractArray{T,N}) where {N, T<:AbstractFloat}
    sz=size(ρ)
    @assert sz == size(λ)
 
    map = zeros(sz)    
    psf_map!(map, h, ρ, λ)
    return map
end    

function spectrum_interpolator(ker::Kernel,   
                               camera_λ_grid::AbstractArray{T,N},                    
                               spectrum_λ_grid::AbstractRange) where {N, T<:AbstractFloat}
    return SparseInterpolator(ker, camera_λ_grid, spectrum_λ_grid)
end

function fit_spectrum!(z::AbstractVector{T},
                       d::AbstractArray{T,N},
                       w::AbstractArray{T,N},
                       H::AbstractArray{T,N},
                       F::SparseInterpolator,
                       μ::T;
                       plot_matrix=false) where {N, T<:AbstractFloat}
    sz=size(d)
    n=length(z)
    @assert sz == size(w)
    @assert sz == size(H)
    @assert μ >=0
    
    A=spzeros(n,n)
    b=zeros(n)
    
    A = LinearInterpolators.SparseInterpolators.AtWA(F,H.*w.*H)
       
    b = F'*(H.*w.*d)
    
    if μ > 0
        A += spdiagm( 0 => 2*μ*ones(n), 1 => -μ*ones(n-1), -1 => -μ*ones(n-1))
    end
    if plot_matrix
        disp_born=max(abs(minimum(A)), maximum(A))
        f1=figure()
        ax=f1.add_subplot(1,1,1)
        I=ax.imshow(A, interpolation="none", cmap=get_cmap("RdGy"), vmin=-disp_born, vmax=disp_born)
        colorbar(I)
    end
    
    z.=A\b
end    

function fit_background!(b::AbstractArray{T,2},
                       d::AbstractArray{T,N},
                       w::AbstractArray{T,N},
                       M::AbstractArray{T,N},
                       μ::T;
                       ratio=1.,
                       maxiter=1000) where {N, T<:AbstractFloat}
    
    function fg_bg!(x::AbstractArray{T,2}, g::AbstractArray{T,2}) where {T <: AbstractFloat}
        fill!(g,0.0);
        r = (x .+ M - d)
        wr = w .*r
        vcopy!(g, sum(wr, dims=3)[:,:])
        f =  vdot(r,wr)/2
        if μ>0
            f += edge_preserving_smoothing(x,g,μ,ratio)
        end
        return f 
    end
    b0=copy(b)
    b .= vmlmb!(fg_bg!, b0;  mem=3, lower=0., maxeval=maxiter, maxiter=maxiter, verb=true)                            
end


function cost(z::AbstractVector{T},
              d::AbstractArray{T,N},
              w::AbstractArray{T,N},
              H::AbstractArray{T,N},
              F::SparseInterpolator) where {N, T<:AbstractFloat}
    sz=size(d)
    @assert sz == size(w)
    @assert sz == size(H)
    
    r = H.*(F*z)
    wr = w.*r
 
    return vdot(r,wr)
end    

 
function cost(z::AbstractVector{T},
              d::AbstractArray{T,N},
              w::AbstractArray{T,N},
              ρ::AbstractArray{T,N},#centré
              λ::AbstractArray{T,N}, 
              h::AbstractPSF,
              F::SparseInterpolator) where {N, T<:AbstractFloat}
    sz=size(d)
    @assert sz == size(w)
    @assert sz == size(ρ)
    @assert sz == size(λ)
    
    H = psf_map(h, ρ, λ)
 
    return cost(z,d,w,H,F)
end    

function cost_tikhonov(z::AbstractVector{T},
              d::AbstractArray{T,N},
              w::AbstractArray{T,N},
              H::AbstractArray{T,N}, 
              F::SparseInterpolator,
              μ::T) where {N, T<:AbstractFloat}
              
    dz = z - cat(z[2:end], 0, dims=1)
    return cost(z, d, w, H, F) + (μ/2)*vdot(dz, dz)
    
end    

function fit_psf_parameters!(z::AbstractVector{T},
                       d::AbstractArray{T,N},
                       w::AbstractArray{T,N},
                       ρ::AbstractArray{T,N},
                       λ::AbstractArray{T,N}, 
                       par::AbstractVector{T},
                       paru::AbstractVector{T},
                       parl::AbstractVector{T},
                       psf_model,
                       F::SparseInterpolator;
                       kwds...) where {N, T<:AbstractFloat}
   
    
    return Bobyqa.minimize!(x->cost(z, d, w, ρ, λ, psf_model(x), F), 
                            par, paru, parl,2.5e-1, 1e-5; kwds...)[2];
end         


function fit_psf_center!(z::AbstractVector{T},
              d::AbstractArray{T,N},
              w::AbstractArray{T,N},
              ρ::AbstractArray{T,N},
              λ::AbstractArray{T,N},
              ρ0::AbstractVector{T},
              ρ0u::AbstractVector{T},
              ρ0l::AbstractVector{T}, 
              h::AbstractPSF,
              F::SparseInterpolator;
              kwds...) where {N, T<:AbstractFloat}
    
    nt = length(ρ0)         
    #@assert nt == size(d)[3]
     return Bobyqa.minimize!(x->cost(z, d, w, ρ .- reshape(x,(1,1,nt)), λ, h, F), 
                            ρ0, ρ0u, ρ0l,5e-1, 1e-5; kwds...)[2];
end         
 

function fit_psf_center!(z::AbstractVector{T},
              d::AbstractArray{T,2},
              w::AbstractArray{T,2},
              ρ::AbstractArray{T,2},
              λ::AbstractArray{T,2},
              ρ0::T,
              h::AbstractPSF,
              F::SparseInterpolator;
              kwds...) where {T<:AbstractFloat}
           
    return Brent.fmin(x->cost(z, d, w, ρ .- x, λ, h, F), 
                            ρ0-1., ρ0+1.; kwds...)[1];
end         


function extract_spectrum!(z::AbstractVector{T},
                       d::AbstractArray{T,N},
                       w::AbstractArray{T,N},
                       ρ::AbstractArray{T,N},
                       λ::AbstractArray{T,N}, 
                       μ::T,
                       psf_params::AbstractVector{T},
                       ρ0::AbstractVector{T},
                       psf_model,
                       F::SparseInterpolator,
                       psf_params_l::AbstractVector{T},
                       psf_params_u::AbstractVector{T},
                       ρ0_dev::AbstractVector{T};
                       ztol=1e-6,
                       maxiter=100,
                       save=false,
                       kwds...) where {N, T<:AbstractFloat}
                       
    f=0.
    fsave=[]
    iter=0
    zres = 1.
    ressave=[]
    nd=ndims(d)
    ztemp=copy(z)
    
    if save
        ndir="temp_results_reconst_"*Dates.format(now(),"yyyy-mm-dd")
        mkdir(ndir)
    end
    
    for iter = 1:maxiter    
        if nd == 3
            ρc = ρ .- reshape(ρ0,(1,1,length(ρ0)));
        elseif nd == 2
            ρc = ρ .- ρ0
        end
        h = psf_model(psf_params)
        H = psf_map(h, ρc, λ)
        z=fit_spectrum!(z, d , w, H, F, μ) 
            
        zres = sum(abs.(z - ztemp))/sum(abs.(z))
          
        if save 
            f=cost_tikhonov(z, d , w, H, F ,μ)  
            push!(fsave,f)
            push!(ressave,zres)
            writedlm(ndir*"/"*"spectum_at_iter_$iter"*".txt", z)
            writedlm(ndir*"/"*"psf_parameters_at_iter_$iter"*".txt", psf_params)
            writedlm(ndir*"/"*"psf_centers_at_iter_$iter"*".txt", ρ0)
            writedlm(ndir*"/"*"cost_function"*".txt", fsave)
            writedlm(ndir*"/"*"relative_absolute_error"*".txt", ressave)
        end
        
        if  zres < ztol
            break
        else
            ztemp=copy(z)
                psf_params = fit_psf_parameters!(z, d  , w, ρc,  λ,  psf_params, 
                                                 psf_params_l, psf_params_u, psf_model, F)
            if nd == 3
                ρ0=fit_psf_center!(z, d , w, ρ, λ,  ρ0, ρ0 .- ρ0_dev, ρ0 .+ ρ0_dev, psf_model(psf_params), F)
            elseif nd == 2
                ρ0[1] =fit_psf_center!(z,  d , w, ρ, λ, ρ0[1], psf_model(psf_params), F)
            end
            
        end    
    end
    return z
end


function extract_spectrum_and_bg!(z::AbstractVector{T},
                       b::AbstractArray{T,2},
                       d::AbstractArray{T,N},
                       w::AbstractArray{T,N},
                       ρ::AbstractArray{T,N},
                       λ::AbstractArray{T,N}, 
                       μ::T,
                       psf_params::AbstractVector{T},
                       ρ0::AbstractVector{T},
                       psf_model,
                       F::SparseInterpolator,
                       psf_params_l::AbstractVector{T},
                       psf_params_u::AbstractVector{T},
                       ρ0_dev::AbstractVector{T};
                       ztol=1e-6,
                       maxiter=100,
                       save=false,
                       kwds...) where {N, T<:AbstractFloat}
                       
    f=0.
    fsave=[]
    iter=0
    zres = 1.
    ressave=[]
    nd=ndims(d)
    ztemp=copy(z)
    
    if save
        ndir="temp_results_reconst_"*Dates.format(now(),"yyyy-mm-dd")
        mkdir(ndir)
    end
    
    for iter = 1:maxiter    
        if nd == 3
            ρc = ρ .- reshape(ρ0,(1,1,length(ρ0)));
        elseif nd == 2
            ρc = ρ .- ρ0
        end
        h = psf_model(psf_params)
        H = psf_map(h, ρc, λ)
        z=fit_spectrum!(z, d .-b , w, H, F, μ) 
        M = H .* (F*z)
        fit_background!(b, d, w, M, μ; ratio=μ*1000, maxiter=10000)
            
        zres = sum(abs.(z - ztemp))/sum(abs.(z))
          
        if save 
            f=cost_tikhonov(z, d .-b , w, H, F ,μ)  
            push!(fsave,f)
            push!(ressave,zres)
            writedlm(ndir*"/"*"spectum_at_iter_$iter"*".txt", z)
            writedlm(ndir*"/"*"psf_parameters_at_iter_$iter"*".txt", psf_params)
            writedlm(ndir*"/"*"psf_centers_at_iter_$iter"*".txt", ρ0)
            writedlm(ndir*"/"*"cost_function"*".txt", fsave)
            writedlm(ndir*"/"*"relative_absolute_error"*".txt", ressave)
        end
        
        if  zres < ztol
            break
        else
            ztemp=copy(z)
                psf_params = fit_psf_parameters!(z, d .-b , w, ρc,  λ,  psf_params, 
                                                 psf_params_l, psf_params_u, psf_model, F)
            if nd == 3
                ρ0=fit_psf_center!(z, d.-b  , w, ρ, λ,  ρ0, ρ0 .- ρ0_dev, ρ0 .+ ρ0_dev, psf_model(psf_params), F)
            elseif nd == 2
                ρ0[1] =fit_psf_center!(z,  d.-b  , w, ρ, λ, ρ0[1], psf_model(psf_params), F)
            end
            
        end    
    end
    return z
end



function extract_spectrum_then_bg!(z::AbstractVector{T},
                       b::AbstractArray{T,2},
                       d::AbstractArray{T,N},
                       w::AbstractArray{T,N},
                       ρ::AbstractArray{T,N},
                       λ::AbstractArray{T,N}, 
                       μ::T,
                       psf_params::AbstractVector{T},
                       ρ0::AbstractVector{T},
                       psf_model,
                       F::SparseInterpolator,
                       psf_params_l::AbstractVector{T},
                       psf_params_u::AbstractVector{T},
                       ρ0_dev::AbstractVector{T};
                       ztol=1e-6,
                       maxiter=100,
                       save=false,
                       kwds...) where {N, T<:AbstractFloat}
                       
    f=0.
    fsave=[]
    iter=0
    zres = 1.
    ressave=[]
    nd=ndims(d)
    ztemp=copy(z)
    
    if save
        ndir="temp_results_reconst_"*Dates.format(now(),"yyyy-mm-dd")
        mkdir(ndir)
    end
    
    for iter = 1:maxiter    
        if nd == 3
            ρc = ρ .- reshape(ρ0,(1,1,length(ρ0)));
        elseif nd == 2
            ρc = ρ .- ρ0
        end
        h = psf_model(psf_params)
        H = psf_map(h, ρc, λ)
        z=fit_spectrum!(z, d , w, H, F, μ) 
            
        zres = sum(abs.(z - ztemp))/sum(abs.(z))
          
        if save 
            f=cost_tikhonov(z, d , w, H, F ,μ)  
            push!(fsave,f)
            push!(ressave,zres)
            writedlm(ndir*"/"*"spectum_at_iter_$iter"*".txt", z)
            writedlm(ndir*"/"*"psf_parameters_at_iter_$iter"*".txt", psf_params)
            writedlm(ndir*"/"*"psf_centers_at_iter_$iter"*".txt", ρ0)
            writedlm(ndir*"/"*"cost_function"*".txt", fsave)
            writedlm(ndir*"/"*"relative_absolute_error"*".txt", ressave)
        end
        
        if  zres < ztol
            break
        else
            ztemp=copy(z)
                psf_params = fit_psf_parameters!(z, d , w, ρc,  λ,  psf_params, 
                                                 psf_params_l, psf_params_u, psf_model, F)
            if nd == 3
                ρ0=fit_psf_center!(z, d , w, ρ, λ,  ρ0, ρ0 .- ρ0_dev, ρ0 .+ ρ0_dev, psf_model(psf_params), F)
            elseif nd == 2
                ρ0[1] =fit_psf_center!(z,  d , w, ρ, λ, ρ0[1], psf_model(psf_params), F)
            end
            
        end    
    end
    return z
end


