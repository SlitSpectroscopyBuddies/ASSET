"""
        oneDimensionalPSF(x) -> h
       
""" oneDimensionalPSF

struct oneDimensionalPSF{V<:AbstractVector} <: NonParametricPSF
    h::V
end

function (P::oneDimensionalPSF)(ρ::AbstractArray{T,N},
                                λ::AbstractArray{T,N}) where {T<:AbstractFloat}
    λref=maxium(λ)
    γ[λ .!=0] .= λref./λ[λ .!=0]
    X = γ.*ρ
    xmin = minimum(X)
    xmax = maximum(X)
    x = range(minimum(X), stop=maximum(X), length=length(P.h))
        return SparseInterpolator(convert(Kernel{T},ker),
                                  convert_eltype(T,X),
                                  convert_eltype(T,x))
end

@inline parameters(P::oneDimensionalPSF) = getfield(P, :h)

function getfwhm(P::chromGaussianPSF, ρ::T,λ::T) where {T<:AbstractFloat}
    @error "Not implented yet"
end

"""
    psf_map!(map, h, ρ, λ)

store in the `AbstractArray` `map` the result of applying the psf function
stored in the `NonParametricPSF` `h`, to each pixel `i` of the spatial and spectral
maps `ρ` and `λ`.

See also [`psf_map`](@ref)

"""
function psf_map!(map::AbstractArray{T,N},
                 h::NonParametricPSF,
                 ρ::AbstractArray{T,N},
                 λ::AbstractArray{T,N}) where {N, T<:AbstractFloat}
    
    @assert axes(map) == axes(ρ) == axes(λ)

        map = h(ρ, λ)*h[:]
end
  
function psf_map(h::NonParametricPSF,
                 ρ::AbstractArray{T,N},
                 λ::AbstractArray{T,N}) where {N, T<:AbstractFloat}

    @assert axes(ρ) == axes(λ)
    
    map = similar(ρ)
    psf_map!(map, h, ρ, λ)
    return map
end    




"""
    fit_spectrum_and_psf!(z, psf, psf_center, F, D, Reg; kwds...)

yields the object spectrum `z`, the off-axis PSF `psf` and its center along the
spectral axis `psf_center`, extracted from the `CalibratedData` `D` via an a
posteriori likelihood minimization with regularization `Reg`. The model of the
object is defined by `Diag(H)*F*z` where `H` is found via the `psf_center` and `psf`
arguments, using the method `psf_map!`. The optimization problem is solved by
the `vmlmb` method defined in the `OptimPackNextGen` package by calling the
method `fit_spectrum!`.

An auto-calibration step can be done to better estimate the parameters, in`psf`, and
center, `psf_center`, of the PSF. The Bobyqa method of Powell is used to
estimate these quantities.

# Keywords
 - `auto_calib` : (`Val(true)` by default) precise if an auto-calibration step
   of the PSF must be done after extracting the spectrum.
 - `psf_params_bnds` : (a vector of zero-values Tuple by default) defines the
   boundaries of the parameters of the PSF. If `auto_calib=true`, the user must
   specify them.
 - `psf_center_bnds` : (a vector of zero-values Tuple by default) defines the
   boundaries of the center of the PSF along the spectral axis. If
   `auto_calib=true`, the user must specify them.
 - `max_iter` : (`1000` by default) defines the maximum number of iterations
   that can do the method (useful when `auto_calib=true`).
 - `loss_tol` : (`(0,1e-6)` by default) defines the absolute and relative
   tolerance between two consecutive iteration of the loss function as a stop
   criterion.
 - `z_tol` : (`(0,1e-6)` by default) defines the absolute and relative
   tolerance between two consecutive iteration of the estimate `z` as a stop
   criterion.
 - Other keywords can be given which are forwarded to the Bobyqa nethod.

See also: [`OptimPackNextGen.vmlbm`](@ref),
[`OptimPackNextGen.Powell.Bobyqa`](@ref), [`psf_map!`](@ref),
[`fit_spectrum!`](@ref)
"""
function fit_spectrum_and_psf!(z::AbstractVector{T},
                      psf::NonParametricPSF,
                      psf_center::AbstractVector{T},
                      F::SparseInterpolator{T},
                      D::CalibratedData{T},
                      Reg::Regularization;
                      auto_calib::Val = Val(true),
                      psf_params_bnds::AbstractVector{Tuple{T,T}} = [(0.,0.) for i in 1:length(psf[:])],
                      psf_center_bnds::AbstractVector{Tuple{T,T}} = [(0.,0.) for i in 1:length(psf_center)],
                      max_iter::Int = 1000,
                      loss_tol::Tuple{Real,Real} = (0.0, 1e-6),
                      z_tol::Tuple{Real,Real} = (0.0, 1e-6),
                      kwds...) where {T}
    
    # Initialization
    iter = 0
    H = similar(D.d)
    ρ_map_centered = D.ρ_map .- reshape(psf_center, 1, 1, length(psf_center))
    psf_map!(H, psf, ρ_map_centered, D.λ_map)
    z_last = copy(z)
    loss_last = loss(D, H, F, z, Reg)
    while true
        # Extract spectrum and PSF FIXME
        
        f=BilinearProblem(D.d,D.w,F, Fz, Reg, Reg);#FIXME regul sur z et regul sur h
info, hopt, zopt = AMORS.solve!(f, h0, z0,μ=μ, ν=ν, first=Val(:x),xtol=1e-3,ytol=1e-3,maxiter=1000)

        # Stop criterions
        if (auto_calib != Val(true)) 
            break
        end
        (auto_calib == Val(true)) && (loss_temp = loss(D, H, F, z, Reg)) 
        if (iter > 1) && ((iter >= max_iter) ||
           test_tol(loss_temp, loss_last, loss_tol) || 
           test_tol(z, z_last, z_tol))
            break
        end
        # Auto-calibration step
        if auto_calib == Val(true)
            check_bnds(psf_params_bnds)
            check_bnds(psf_center_bnds)
            
            psf = fit_psf_params(psf, psf_center, z, F, D; 
                                 psf_params_bnds=psf_params_bnds)
            fit_psf_center!(psf_center, psf, z, F, D;
                            psf_center_bnds=psf_center_bnds)
            # re-center the spatial map with the new centers of the psf
            copyto!(ρ_map_centered, D.ρ_map .- reshape(psf_center, 1, 1, length(psf_center)))
            psf_map!(H, psf, ρ_map_centered, D.λ_map)
        end
        iter +=1
        display(loss_temp)
        loss_last = loss_temp
        copyto!(z_last, z)
    end
    return z, psf, psf_center
end

