

"""
    extract_spectrum!(z, F, psf_center, psf, D, Reg [, Bkg]; kwds...)

estimates the spectrum of an object observed with long-slit spectroscopy,
when the data can be corrupted by a background component and noise. The direct model
of such data can be writen:
`
                    d = Bkg + Diag(H)*F*z + n
`
where `*` the matrix product. `z` is the spectrum of the object of interest and
`F` its associated interpolation operator, `H` the operator modeling the Point
Spread Function (PSF) of the instrument, while `Bkg` is the background component
to disentangle from the object of interest and `n` accounts for noises.

The data and there associated weights are given by the `CalibratedData`
structure `D`, while the array `H` is produced using the method `psf_map` with
the arguments `psf` and `psf_center`.

The estimation of the different unknowns relies on the minimization of the a
posteriori likelihood, where `Reg` is the regularization function of `z` (of
type `InverseProblem.Regularization`) and where `Bkg` is a structure of type
`AbstractBkg`, which contains a way to produce the background component and its
regularization function. If `Bkg` is not given, the algorithm will suppose that
there is no background in the data and will only call `object_step!`. If there
is a `Bkg` given, the algorithm will alternate the estimation of the background
using `fit_bkg!` with the extraction of the spectrum (done in `object_step!`).

Finally, an auto-calibration step can be activated which will refine the
parameters and center of the PSF.

# Keywords
 - `auto_calib` : (`Val(true)` by default) precise if an auto-calibration step
   of the PSF must be done after extracting the spectrum. It can also take the
   value `Val(:delay)` which will run the algorithm without auto_calib until it
   converges, before activating the `auto_calib` and re-running the algorithm.
 - `mask_width` : (`3` by default) defines the number of fwhm of the psf will be
   used to hide the object on the first iteration of the background estimtion
   (if thee is one).
 - `max_iter` : (`1000` by default) defines the maximum number of iterations
   that can do the method (useful when `auto_calib=true`).
 - `loss_tol` : (`(0,1e-6)` by default) defines the absolute and relative
   tolerance between two consecutive iteration of the loss function as a stop
   criterion.
 - `z_tol` : (`(0,1e-6)` by default) defines the absolute and relative
   tolerance between two consecutive iteration of the estimate `z` as a stop
   criterion.
 - `extract_kwds` : (`(verb=true,)` by default) where other keywords can be
   given which are forwarded to the `object_step` method. 

See also [`AbstractBkg`](@ref), [`AbstractPSF`](@ref)

#FIXME: size(psf_center) = (1, 1, size(ρ_map,3)) to specify in the doc?
#FIXME: keywords for psf_params_bnds and psf_center_bnds to detail also here to
make sure the user specifies them if auto_calib?
#FIXME: change object_step! name in doc?
"""
function extract_spectrum!(z::AbstractVector{T},
    F::SparseInterpolator{T},
    psf_center::AbstractArray{T},
    psf::AbstractPSF,
    D::CalibratedData{T},
    Reg::Regularization,
    Bkg::Union{<:AbstractBkg,UndefInitializer} = undef;
    auto_calib::Val = Val(true),
    mask_width::Real = 3.0,
    max_iter::Int = 1000,
    loss_tol::Tuple{Real,Real} = (0.0, 1e-6),
    z_tol::Tuple{Real,Real} = (0.0, 1e-6),
    extract_kwds::K = (verb=true,)) where {T,K<:NamedTuple}
    
    @assert size(D.d) == LinearInterpolators.output_size(F)
    
    # Initialization
    iter = 0
    Res = CalibratedData(copy(D.d), copy(D.w), copy(D.ρ_map), copy(D.λ_map))
    mask_z = similar(Res.w)
    H = similar(Res.d)
    ρ_map_centered = Res.ρ_map .- reshape(psf_center, 1, 1, length(psf_center))
    # psf_param = collect(psf[:])#FIXME: doesn't seem to be used, remove?
    psf_map!(H, psf, ρ_map_centered, Res.λ_map)
    z_last = copy(z)
    loss_last = loss(D, H, F, z, Reg, Bkg)
    while true
        if Bkg != undef
            # Estimate the background component in the data
            if iter == 0
                # Mask object for first iteration
                copyto!(mask_z, mask_object(ρ_map_centered, Res.λ_map, psf; 
                                            mask_width=mask_width))
                
                res_bkg_masked_object = CalibratedData(D.d, D.w.*mask_z, D.ρ_map, D.λ_map) 
                fit_bkg!(Bkg, res_bkg_masked_object)                                            
            else
                fit_bkg!(Bkg, Res)
            end
            # copy the background subtracted to data to Res for object spectrum extraction
            copyto!(Res.d, D.d - Bkg)
        end

        # Estimate object's spectrum FIXME: change name?
        #FIXME: check psf_center is updated after this step?
        psf = object_step!(z, psf, psf_center, F, Res, Reg; 
                           auto_calib=auto_calib, extract_kwds...)[2]
        # re-center the spatial map with the new centers of the psf
        copyto!(ρ_map_centered, Res.ρ_map .- reshape(psf_center, 1, 1, length(psf_center)))
        psf_map!(H, psf, ρ_map_centered, Res.λ_map)

        # Stop criterions
        loss_temp = loss(D, H, F, z, Reg, Bkg)
        if (iter >= max_iter) || test_tol(loss_temp, loss_last, loss_tol) || 
                                 test_tol(z, z_last, z_tol)
            if auto_calib == Val(:delay)
                auto_calib = Val(true)
            else
                break
            end
        end
        copyto!(Res.d, D.d - H .* (F*z))
        iter += 1
        loss_last = loss_temp
        copyto!(z_last, z)
    end

    return z, psf, psf_center
end




"""
    mask_object(ρ_map, λ_map, psf [; mask_width=3])

yields a mask array of same size than `ρ_map` and `λ_map`, where all the pixels
of distance less than the fwhm of `psf` times the `mask_width` are flagged as
zeros, while the rest are at unitary level.

To do so, the user needs to make sure that the `ρ_map` has its origin centered
on the object of interest.

See also [`AbstractPSF`](@ref), [`get_fwhm`](@ref)

"""
function mask_object(ρ_map::AbstractArray{T,N},#must be centered
                     λ_map::AbstractArray{T,N},
                     psf::AbstractPSF;
                     mask_width::Union{T,AbstractVector{T}} = 3.0) where {T,N}
    
    @assert axes(ρ_map) == axes(λ_map)
    
    mask_z = ones(T, size(ρ_map))
    for i in eachindex(ρ_map)
        fwmh=getfwhm(psf, ρ_map[i], λ_map[i])
        if psf(ρ_map[i], λ_map[i]) <= mask_width * fwmh
            mask_z[i] = T(0)
        end
    end

    return mask_z
end




"""
    object_step!(z, psf, psf_center, F, D, Reg; kwds...)

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

FIXME: change name into object_spectrum_extraction?
"""
function object_step!(z::AbstractVector{T},
                      psf::AbstractPSF,
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
                      kwds...) where {T,N}
    
    # Initialization
    iter = 0
    H = similar(D.d)
    ρ_map_centered = D.ρ_map .- reshape(psf_center, 1, 1, length(psf_center))
    psf_map!(H, psf, ρ_map_centered, D.λ_map)
    z_last = copy(z)
    loss_last = loss(D, H, F, z, Reg)
    while true
        # Extract spectrum
        fit_spectrum!(z, F, H, D, Reg; kwds...)
        # Stop criterions
        (auto_calib == Val(true)) && (loss_temp = loss(D, H, F, z, Reg)) 
        if (auto_calib != Val(true)) || (iter >= max_iter) || 
                                        test_tol(loss_temp, loss_last, loss_tol) || 
                                        test_tol(z, z_last, z_tol)
            break
        end
        # Auto-calibration step
        if auto_calib == Val(true)
            psf = fit_psf_params(psf, psf_center, z, F, D; 
                                 psf_params_bnds=psf_params_bnds)
            fit_psf_center!(psf_center, psf, z, F, D;
                            psf_center_bnds=psf_center_bnds)
            # re-center the spatial map with the new centers of the psf
            copyto!(ρ_map_centered, D.ρ_map .- reshape(psf_center, 1, 1, length(psf_center)))
            psf_map!(H, psf, ρ_map_centered, D.λ_map)
        end
        iter +=1
        loss_last = loss_temp
        copyto!(z_last, z)
    end

    return z, psf, psf_center
end




"""
    psf_map!(map, h, ρ, λ)

store in the `AbstractArray` `map` the result of applying the psf function
stored in the `AbstractPSF` `h`, to each pixel `i` of the spatial and spectral
maps `ρ` and `λ`.

See also [`psf_map`](@ref)

"""
function psf_map!(map::AbstractArray{T,N},
                 h::AbstractPSF,
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
function psf_map(h::AbstractPSF,
                 ρ::AbstractArray{T,N},
                 λ::AbstractArray{T,N}) where {N, T<:AbstractFloat}

    @assert axes(ρ) == axes(λ)
    
    map = similar(ρ)
    psf_map!(map, h, ρ, λ)
    return map
end    




"""
    fit_spectrum!(z, F, H, D, Reg ; kwds...)

check if the `Regularization` `R` can be used in a direct inversion framework
before calling the solve function (either `solve_analytic` or `solve_vmlmb`).
The keywords `kwds` are forwarded to the method that solves the problem.

See also [`InverseProblem.Regul`](@ref)

"""
function fit_spectrum!(z::AbstractVector{T},
                       F::SparseInterpolator{T},
                       H::AbstractArray{T,N},
                       D::CalibratedData{T},
                       Reg::Regularization;
                       kwds...) where {T,N}

    if InverseProblem.use_direct_inversion(Reg)
        copyto!(z, solve_analytic!(F, H, D, Reg))#FIXME: no `!` here?
    else
        copyto!(z, solve_vmlmb!(z, Diag(H) * F, D, Reg; kwds...))#FIXME: no `!` here?
    end
    
    return z
end




"""
    solve_analytic!(F, H, D, Reg)

yields the estimator by directly inverting the Normal equations:
`
        (F'*H'*Diag(D.w)*H*F + get_grad_op(R)) * z = F'.H*Diag(D.w)*D.d
`
where `D` is a `CalibratedData`, `get_grad_op` is a method returning the
operator of the gradient of the `Regularization` `Reg` and `z` is the estimator
the user is looking for.

See also [`InversePbm.get_grad_op`](@ref)

#FIXME: no `!` here?
"""
function solve_analytic!(F::SparseInterpolator{T},
                         H::AbstractArray{T,N},
                         D::CalibratedData,
                         Reg::Regularization) where {T,N}
    
    d, w = D.d, D.w
    sz=size(d)
    n=length(z)

    @assert sz == size(w)
    @assert sz == size(H)
    @assert μ >=0
    
    A=spzeros(n,n)
    
    A = LinearInterpolators.SparseInterpolators.AtWA(F,H.*w.*H)
    A += InversePbm.get_grad_op(Reg)

    b=zeros(n)
    b = F'*(H.*w.*d)

    return A\b
end


"""
    solve_vmlmb!(z0, A, D, Reg; [nonnegative=false,] kwds...)

uses the `vmlmb!` method of `OptimPackNextGen` to estimate the solution of:
`
            argmin (A*z - D.d)'*Diag(D.w)*(A*z - D.d) + Reg(z)
               z
`
with `D` a `CalibratedData`. The result is stored in `z0`.

It is possible to indicate to the method `vmlmb` a positivity constraint for
`z` by using the `nonnegative` keyword, as well as indicate more keywords to
constrain the optimization.

See also [`OptimPackNextGen.vmlbm!`](@ref)

#FIXME: no `!` here?
"""
function solve_vmlmb!(z0::AbstractVector{T},
                      A::Mapping,
                      D::CalibratedData{T},
                      Reg::Regularization;
                      nonnegative::Bool = false,
                      kwds...) where {T}
    
    d, w = D.d, D.w

    @assert is_linear(A)

    IP = InvProblem(A, d, w, Reg)
    function fg_solve!(z, g)
        return call!(IP, z, g)
    end
    
    if nonnegative
        vmlmb!(fg_solve!, z0; lower=T(0), kwds...)
    else
        vmlmb!(fg_solve!, z0; kwds...)
    end
    
    return z0
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
function fit_psf_params(psf::AbstractPSF,
                        psf_center::AbstractVector{T},
                        z::AbstractVector{T},
                        F::SparseInterpolator{T},
                        D::CalibratedData{T};
                        psf_params_bnds::AbstractVector{Tuple{T,T}} = [(0.,0.) for i in 1:length(psf[:])],
                        rho_tol::Union{UndefInitializer,Tuple{Real,Real}} = undef,
                        kwds...) where {T,N}
    
    @assert length(psf[:]) == length(psf_params_bnds)

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

    # create a vector containing the parameters of the PSF and calling Bobyqa
    psf_param = collect(psf[:])                                  
    Bobyqa.minimize!(p -> likelihood!(p), psf_param, bnd_min, bnd_max, 
                     rho_tol[1], rho_tol[2]; kwds...)[2]
    
    return typeof(psf)(psf_param)
end




"""
    fit_psf_center(psf_center, psf, z, F, D; kwds...)

yields the center of the psf along the spectral axis, stored in `psf_center`,
estimated by minimizing the likelihood formed with the `CalibratedData` `D`, the
spectrum of the object `z` and the `SparseInterpolator` `F`.

#Keywords
 - `psf_center_bnds` : (a vector of zero-values Tuple by default) defines the
   boundaries of the center of the PSF along the spectral axis. The user must
   specify them.
 - `rho_tol` : (`undef` by default) defines the size of the trust-region used to
   estimate the parameters of the PSF.
 - Other keywords can be given which are forwarded to the Bobyqa nethod.

See also [`OptimPackNextGen.Powell.Bobyqa`](@ref)

"""
function fit_psf_center!(psf_center::AbstractVector{T},
                         psf::AbstractPSF,
                         z::AbstractVector{T},
                         F::SparseInterpolator{T},
                         D::CalibratedData{T};
                         psf_center_bnds::AbstractVector{Tuple{T,T}} = [(0.,0.) for i in 1:length(psf_center)],
                         rho_tol::Union{UndefInitializer,Tuple{Real,Real}} = undef,
                         kwds...) where {T,N}
    
    @assert length(psf_center) == length(psf_center_bnds)

    # define the size of the trust region for bobyqa
    if rho_tol == undef
        rhobeg=(psf_center_bnds[1][2] - psf_center_bnds[1][1])/2
        for i = 2:length(psf_center_bnds)
            rhobeg = min(rhobeg, (psf_center_bnds[i][2] - psf_center_bnds[i][1])/2)
        end
        rho_tol = (0.99*rhobeg, 1e-3*rhobeg)
    end
    # define boundaries to use for bobyqa
    bnd_min = zeros(length(psf_center))
    bnd_max = zeros(length(psf_center))
    for p in eachindex(psf_center_bnds)
        bnd = psf_center_bnds[p]
        bnd_min[p] = bnd[1]
        bnd_max[p] = bnd[2]
    end

    # define the likelihood to be minimized
    d, w, ρ_map, λ_map = D.d, D.w, D.ρ_map, D.λ_map
    ρ_map_c = zeros(T, size(d))
    H_c = zeros(T, size(d))
    function likelihood!(c)
        ρ_map_c = ρ_map .- reshape(c, 1, 1, length(psf_center))#FIXME: here c was replaced by psf_center!!!!! Check if still working
        psf_map!(H_c, psf, ρ_map_c, λ_map)
        L = Lkl(Diag(H_c) * F, d, w)
        return L(z)
    end


    return Bobyqa.minimize!(c -> likelihood!(c), psf_center, bnd_min, bnd_max, 
                            rho_tol[1], rho_tol[2]; kwds...)[2]
end



"""
    loss(D, H, F, z, Reg [, Bkg=undef])

yields the value of the loss function used to estimate the different parameters
of the problem:
`
    (Diag(H)*F*z + Bkg.b - D.d)'*Diag(D.w)*(Diag(H)*F*z + Bkg.b - D.d) + Reg(z) + regul(Bkg)
`
where `Bkg.b` and `regul(Bkg)` are not taken into account if `Bkg=undef`.

See also [`AbstractBkg`](@ref), [`InverseProblem.Regul`](@ref)

"""
function loss(D::CalibratedData{T},
              H::AbstractArray{T,N},
              F::SparseInterpolator{T},
              z::AbstractVector{T},
              Reg::Regularization,
              Bkg::Union{AbstractBkg,UndefInitializer} = undef) where {T,N}

    # Computing Likelihood
    wks = vcopy(D.d)
    wks -= H .* (F*z)
    if Bkg != undef
        wks .-= Bkg.b
    end
    lkl = vdot(wks, D.w .* wks)

    #Computing Regularization
    r = Reg(z)
    if Bkg != undef
        r += regul(Bkg)
    end

    return lkl + r
end
