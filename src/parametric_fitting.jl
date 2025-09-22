#
# parametric_fitting.jl
#
# ------------------------------------------------
#
# This file is part of ASSET


"""
    psf_map!(map, h, ρ, λ)

store in the `AbstractArray` `map` the result of applying the psf function
stored in the `ParametricPSF` `h`, to each pixel `i` of the spatial and spectral
maps `ρ` and `λ`.

See also [`psf_map`](@ref)

"""
function psf_map!(map::AbstractArray{T,N},
                 h::ParametricPSF,
                 ρ::AbstractArray{T,N},
                 λ::AbstractArray{T,N}) where {N, T<:AbstractFloat}
    
    @assert axes(map) == axes(ρ) == axes(λ)

    @inbounds for i in eachindex(map, ρ, λ)
        map[i] = h(ρ[i], λ[i])
    end
end

"""
    psf_map(h, ρ, λ)

yield the result of `psf_map!` and store it in a new `AbstractArray`.

See also [`psf_map!`](@ref)

"""
function psf_map(h::ParametricPSF,
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
                      psf::ParametricPSF,
                      #psf_center::AbstractVector{T},
                      F::SparseInterpolator{T},
                      D::CalibratedData{T},
                      Reg::Regularization;
                      auto_calib::Val = Val(true),
                      psf_params_bnds::AbstractVector{Tuple{T,T}} = [(0.,0.) for i in 1:length(psf[:])],
                      #psf_center_bnds::AbstractVector{Tuple{T,T}} = [(0.,0.) for i in 1:length(psf_center)],
                      psf_center_bnd::T=0.,
                      max_iter::Int = 1000,
                      loss_tol::Tuple{Real,Real} = (0.0, 1e-6),
                      z_tol::Tuple{Real,Real} = (0.0, 1e-6),
                      kwds...) where {T}
    
    # Initialization
    iter = 0
    H = similar(D.d)
    #ρ_map_centered = D.ρ_map .- reshape(psf_center, 1, 1, length(psf_center))
    #psf_map!(H, psf, ρ_map_centered, D.λ_map)
    psf_map!(H, psf, D.ρ_map, D.λ_map)
    z_last = copy(z)
    #loss_last = loss(CalibratedData(D.d, D.w, ρ_map_centered, D.λ_map), psf, F, z, Reg)
    loss_last = loss(CalibratedData(D.d, D.w, D.ρ_map, D.λ_map), psf, F, z, Reg)
    axs_D = size(D.ρ_map)
    psf_center = ( length(axs_D) == 3 ? zeros(axs_D[3]) : zeros(1) )
        psf_center_bnds = [(-psf_center_bnd, psf_center_bnd) for k=1:length(psf_center)]
            
        #display(psf)
    while true
        # Extract spectrum
        fit_spectrum!(z, F, H, D, Reg; kwds...)
        # Stop criterions
        if (auto_calib != Val(true)) 
            break
        end
        #(auto_calib == Val(true)) && (loss_temp = loss(CalibratedData(D.d, D.w, ρ_map_centered, D.λ_map), psf, F, z, Reg)) 
        (auto_calib == Val(true)) && (loss_temp = loss(CalibratedData(D.d, D.w, D.ρ_map, D.λ_map), psf, F, z, Reg)) 
        if (iter > 0) && ((iter >= max_iter) ||
           test_tol(loss_temp, loss_last, loss_tol) || 
           test_tol(z, z_last, z_tol))
            break
        end
        # Auto-calibration step
        if auto_calib == Val(true)
            check_bnds(psf_params_bnds)
            check_bnds(psf_center_bnds)
            fill!(psf_center,0.)
            fit_psf_center!(psf_center, psf, z, F, D;
                            psf_center_bnds=psf_center_bnds)
            #display(psf_center)
            psf = fit_psf_params(psf, psf_center, z, F, D; 
                                 psf_params_bnds=psf_params_bnds)
            #display(psf)
            # re-center the spatial map with the new centers of the psf
            #copyto!(ρ_map_centered, D.ρ_map .- reshape(psf_center, 1, 1, length(psf_center)))
            #psf_map!(H, psf, ρ_map_centered, D.λ_map)
            copyto!(D.ρ_map, D.ρ_map .- reshape(psf_center, 1, 1, length(psf_center)))
            psf_map!(H, psf, D.ρ_map, D.λ_map)
        
        end
        iter +=1
        loss_last = loss_temp
        copyto!(z_last, z)
    end
        #display(psf)
    return z, psf, psf_center
end




"""
    fit_psf_params(psf, psf_center, z, F, D; kwds...)

yields a new PSF structure of same type than `psf`, where its parameters are
estimated by Powell's Bobyqa method as defined in `OptimPackNextGen`. The center
of the PSF is defined by the vector `psf_center`, while the spectrum of the
object of interest can be retireve using the vector `z`, the `SparseInterpolator
`F`, and the `CalibratedData` `D`.

#Keywords
 - `psf_params_bnds` : (a vector of zero-values Tuple by default) defines the
   boundaries of the parameters of the PSF. The user must specify them.
 - `rho_tol` : (`undef` by default) defines the size of the trust-region used to
   estimate the parameters of the PSF.
 - Other keywords can be given which are forwarded to the Bobyqa nethod.

See also [`OptimPackNextGen.Powell.Bobyqa`](@ref)

"""
function fit_psf_params(psf::ParametricPSF,
                        psf_center::AbstractVector{T},
                        z::AbstractVector{T},
                        F::SparseInterpolator{T},
                        D::CalibratedData{T};
                        psf_params_bnds::AbstractVector{Tuple{T,T}} = [(0.,0.) for i in 1:length(psf[:])],
                        rho_tol::Union{UndefInitializer,Tuple{Real,Real}} = undef,
                        kwds...) where {T}
    
    @assert length(psf[:]) == length(psf_params_bnds)
    
    # create a vector containing the parameters of the PSF
    # define the size of the trust region for bobyqa
    if rho_tol == undef               
        rhobeg=(psf_params_bnds[1][2] - psf_params_bnds[1][1])/2                                     
        for i = 2:length(psf_params_bnds)
            rhobeg = min(rhobeg, (psf_params_bnds[i][2] - psf_params_bnds[i][1])/2)
        end
        rho_tol = (0.99*rhobeg, 1e-3*rhobeg)
    end
    # define boundaries to use for bobyqa
    bnd_min = zeros(length(psf_params_bnds))
    bnd_max = zeros(length(psf_params_bnds))
    for p in eachindex(psf_params_bnds)
        bnds_param = psf_params_bnds[p]
        bnd_min[p] = bnds_param[1]
        bnd_max[p] = bnds_param[2]
    end

    # define the likelihood to be minimized
    d, w, ρ_map, λ_map = D.d, D.w, D.ρ_map, D.λ_map
    ρ_map_centered = ρ_map .- reshape(psf_center, 1, 1, length(psf_center))
    H_p = zeros(T, size(d))
    function likelihood!(p)
        psf_p = typeof(psf)(p)
        psf_map!(H_p, psf_p, ρ_map_centered, λ_map)
        L = Lkl(Diag(H_p) * F, d, w)
        return L(z)
    end
    psf_param = collect(psf[:])  
    if length(psf_param)<2
        psf_param=fmin(likelihood!,bnd_min[1],bnd_max[1];kwds...)[1]
    else
    #    status, psf_param, fp = bobyqa!(likelihood!, psf_param, bnd_min, bnd_max, rho_tol[1], rho_tol[2]; kwds...)
    psf_param = bobyqa(likelihood!, psf_param, xl=bnd_min, xu=bnd_max, rhobeg=rho_tol[1], rhoend=rho_tol[2]; kwds...)[1]
    end
    return typeof(psf)(psf_param)
end


