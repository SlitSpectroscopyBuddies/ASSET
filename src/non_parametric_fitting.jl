#
# non_parametric_fitting.jl
#
# ------------------------------------------------
#
# This file is part of ASSET


"""
    psf_map!(map, h, ρ, λ)

store in the `AbstractArray` `map` the result of applying the psf function
stored in the `NonParametricPSF` `h`, to each pixel `i` of the spatial and spectral
maps `ρ` and `λ`.

See also [`psf_map`](@ref)

"""
function psf_map!(map::AbstractArray{T,N},
                  P::NonParametricPSF,
                  ρ::AbstractArray{T,N},
                  λ::AbstractArray{T,N};
                  λref=maximum(λ)) where {N,T<:AbstractFloat}
    
    @assert axes(map) == axes(ρ) == axes(λ)
    
    map .= P(ρ, λ;λref)*P[:]
end

"""
psf_map(h, ρ, λ)

yield the result of `psf_map!` and store it in a new `AbstractArray`.

See also [`psf_map!`](@ref)

"""
function psf_map(P::NonParametricPSF,
                 ρ::AbstractArray{T,N},
                 λ::AbstractArray{T,N};
                 λref=maximum(λ)) where {N,T<:AbstractFloat}

    @assert axes(ρ) == axes(λ)
    
    map = similar(ρ)
    psf_map!(map, P, ρ, λ;λref)
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
#FIXME: update doc
"""
function fit_spectrum_and_psf!(z::AbstractVector{T},
    psf::NonParametricPSF,
    psf_center::AbstractVector{T},
    F::SparseInterpolator{T},
    D::CalibratedData{T},
    Reg::Regularization;
    auto_calib::Val = Val(true),
    shift_bnds::Tuple{T,T} = (0.,0.),
    psf_center_bnds::AbstractVector{Tuple{T,T}} = [(0.,0.) for i in 1:length(psf_center)],
    max_iter::Int = 1000,
    loss_tol::Tuple{Real,Real} = (0.0, 1e-6),
    z_tol::Tuple{Real,Real} = (0.0, 1e-6),
    h_tol::Tuple{Real,Real} = (0.0, 1e-6),
    kwds...) where {T}
    
    # Initialization
    iter = 0
    #H = similar(D.d)
    ρ_map_centered = D.ρ_map .- reshape(psf_center, 1, 1, length(psf_center))
    #psf_map!(H, psf, ρ_map_centered, D.λ_map)
    z_last = copy(z)
    h_last=copy(psf[:])
    loss_last = loss(CalibratedData(D.d, D.w, ρ_map_centered, D.λ_map), psf, F, z, Reg)
    Fh = psf(ρ_map_centered, D.λ_map)
    while true
        f = BilinearProblem(D.d, D.w, Fh,  F, psf.R, Reg);
        # Extract spectrum and PSF
        AMORS.solve!(f, psf[:], z, first=Val(:x),
                                        xtol=1e-3,ytol=1e-3,maxiter=1000)

        # Stop criterions
        if (auto_calib != Val(true)) 
            break
        end
        (auto_calib == Val(true)) && (loss_temp = loss(CalibratedData(D.d, D.w, ρ_map_centered, D.λ_map), psf, F, z, Reg)) 
        if (iter > 1) && ((iter >= max_iter) ||
           test_tol(loss_temp, loss_last, loss_tol) || 
           test_tol(z, z_last, z_tol) ||
           test_tol(psf[:], h_last, h_tol))
            break
        end
        # Auto-calibration step
        if auto_calib == Val(true)
            #FIXME: seems not to be working
            check_bnds(psf_center_bnds)            
            fit_psf_center!(psf_center, psf, z, F, D;
                            psf_center_bnds=psf_center_bnds)
            # re-center the spatial map with the new centers of the psf
            copyto!(ρ_map_centered, D.ρ_map .- reshape(psf_center, 1, 1, length(psf_center)))
            #psf_map!(H, psf, ρ_map_centered, D.λ_map)
            Fh = psf(ρ_map_centered, D.λ_map)
        end
        iter +=1
        loss_last = loss_temp
        copyto!(z_last, z)
        copyto!(h_last, psf[:])
    end
    return z, psf, psf_center
end
#=
function fit_psf_shift!(psf::NonParametricPSF,
                        psf_center::AbstractVector{T},
                        z::AbstractVector{T},
                        F::SparseInterpolator{T},
                        D::CalibratedData{T};
                        shift_bnds::Tuple{T,T} = (0.,0.),
                        rho_tol::Union{UndefInitializer,Tuple{Real,Real}} = undef,
                        kwds...) where {T}
                        
    # create a vector containing the parameters of the PSF
    # define the size of the trust region for bobyqa
     # define boundaries to use for bobyqa
    bnd_min = shift_bnds[1]
    bnd_max = shift_bnds[2]

    # define the likelihood to be minimized
    d, w, ρ_map, λ_map = D.d, D.w, D.ρ_map, D.λ_map
    ρ_map_centered = ρ_map .- reshape(psf_center, 1, 1, length(psf_center))
    H_p = zeros(T, size(d))
    function likelihood!(a)
        psf_p = typeof(psf)(psf[:],a, psf.ker,psf.R)
        psf_map!(H_p, psf_p, ρ_map_centered, λ_map)
        L = Lkl(Diag(H_p) * F, d, w)
        return L(z)
    end
    shift_param = shift(psf)  
    shift_param=fmin(likelihood!,bnd_min,bnd_max;kwds...)[1]
    end
    return typeof(psf)(psf_param)
end
=#
