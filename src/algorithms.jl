#
# algorithms.jl
#
# ------------------------------------------------
#
# This file is part of ASSET


"""
    extract_spectrum!(z, F, psf, D, Reg [, Bkg]; kwds...)

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
there is no background in the data and will only call `fit_spectrum_and_psf!`. If there
is a `Bkg` given, the algorithm will alternate the estimation of the background
using `fit_bkg!` with the extraction of the spectrum (done in `fit_spectrum_and_psf!`).

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

# See also
- [`AbstractBkg`](@ref)
- [`PointSpreadFunctions.AbstractPSF`]
"""
function extract_spectrum!(z::AbstractVector{T},
    F::SparseInterpolator{T},
    psf::AbstractPSF,
    D::CalibratedData{T},
    Reg::Regularization,#TODO: multiple regul for multiple params?
    Bkg::Union{<:AbstractBkg,UndefInitializer} = undef;
    auto_calib::Val = Val(true),
    mask_width::Real = 3.0,
    max_iter::Int = 1000,
    loss_tol::Tuple{Real,Real} = (0.0, 1e-6),
    z_tol::Tuple{Real,Real} = (0.0, 1e-6),
    iterate_memory::Integer = 3, #the last iteration to compare to avoid eternal loops
    extract_kwds::K = (verb=true,)) where {T,K<:NamedTuple}
    
    @assert size(D.d) == LinearInterpolators.output_size(F)
    
    # Initialization
    iter = 0
    Res = CalibratedData(copy(D.d), copy(D.w), copy(D.ρ_map), copy(D.λ_map))
    mask_z = similar(Res.w)
    H = similar(Res.d)
    psf_map!(H, psf, Res.ρ_map, Res.λ_map)
    z_last = zeros(length(z),iterate_memory)
    z_last[:,iterate_memory] .= copy(z)
    loss_last = zeros(iterate_memory)
    loss_last[iterate_memory] = loss(CalibratedData(Res.d, Res.w, Res.ρ_map, Res.λ_map), psf, F, z, Reg, Bkg)
    while true
        if Bkg != undef
            copyto!(Res.d, D.d - H .* (F*z))
            # Estimate the background component in the data
            if iter == 0
                # Mask object for first iteration
                copyto!(mask_z, mask_object(Res.ρ_map, Res.λ_map, psf; 
                                            mask_width=mask_width))
                
                res_bkg_masked_object = CalibratedData(D.d, D.w.*mask_z, D.ρ_map, D.λ_map) 
                fit_bkg!(Bkg, res_bkg_masked_object)                                            
            else
                fit_bkg!(Bkg, Res)
            end
            # copy the background subtracted to data to Res for object spectrum extraction
            copyto!(Res.d, D.d - Bkg)
        end

        # Estimate object's spectrum and autocalibrate object psf parameters 
        psf = fit_spectrum_and_psf!(z, psf, F, Res, Reg; 
                           auto_calib=auto_calib, extract_kwds...)[2]
                           
        # re-center the spatial map with the new centers of the psf
        psf_map!(H, psf, Res.ρ_map, Res.λ_map)


        # Stop criterions
        loss_temp = loss(CalibratedData(Res.d, Res.w, Res.ρ_map, Res.λ_map), psf, F, z, Reg, Bkg)

       if (iter >= max_iter) || test_tol(loss_temp, loss_last, loss_tol) || 
                                 test_tol(z, z_last, z_tol) 
            if auto_calib == Val(:delay)
                auto_calib = Val(true)
            else
                break
            end
        end
        iter += 1
        loss_last[mod(iter-1,iterate_memory)+1] = loss_temp
        z_last[:,mod(iter-1,iterate_memory)+1] .= z
    end
    
    return z, psf
end




"""
    mask_object(ρ_map, λ_map, psf [; mask_width=3])

yields a mask array of same size than `ρ_map` and `λ_map`, where all the pixels
of distance less than the fwhm of `psf` times the `mask_width` are flagged as
zeros, while the rest are at unitary level.

To do so, the user needs to make sure that the `ρ_map` has its origin centered
on the object of interest.

# See also
- [`PointSpreadFunctions.AbstractPSF`]
- [`PointSpreadFunctions.get_fwhm`]
"""
function mask_object(ρ_map::AbstractArray{T,N},
                     λ_map::AbstractArray{T,N},
                     psf::AbstractPSF;
                     mask_width::Union{T,AbstractVector{T}} = 3.0) where {T,N}
    
    @assert axes(ρ_map) == axes(λ_map)
    
    mask_z = ones(T, size(ρ_map))
    map = psf_map(psf, ρ_map, λ_map)
    for i in eachindex(ρ_map, λ_map, map, mask_z)
        fwmh=getfwhm(psf, ρ_map[i], λ_map[i])
        if map[i] <= mask_width * fwmh 
            mask_z[i] = T(0)
        end
    end

    return mask_z
end




"""
    fit_spectrum!(z, F, H, D, Reg ; kwds...)

check if the `Regularization` `R` can be used in a direct inversion framework
before calling the solve function (either `solve_analytic` or `solve_vmlmb`).
The keywords `kwds` are forwarded to the method that solves the problem.

# See also
- [`InverseProblem.Regul`]
"""
function fit_spectrum!(z::AbstractVector{T},
                       F::SparseInterpolator{T},
                       H::AbstractArray{T,N},
                       D::CalibratedData{T},
                       Reg::Regularization;
                       kwds...) where {T,N}

    if InverseProblem.use_direct_inversion(Reg)
        copyto!(z, solve_analytic(F, H, D, Reg))
    else
        copyto!(z, solve_vmlmb(z, Diag(H) * F, D, Reg; kwds...))
    end
    
    return z
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

# See also
- [`OptimPackNextGen.Powell.Bobyqa`]
"""
function fit_psf_center!(psf_center::AbstractVector{T},
                         psf::AbstractPSF,
                         z::AbstractVector{T},
                         F::SparseInterpolator{T},
                         D::CalibratedData{T};
                         psf_center_bnds::AbstractVector{Tuple{T,T}} = [(0.,0.) for i in 1:length(psf_center)],
                         rho_tol::Union{UndefInitializer,Tuple{Real,Real}} = undef,
                         kwds...) where {T}
    
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
        if length(psf_center) > 1
        ρ_map_c = ρ_map .- reshape(c, 1, 1, length(psf_center))
        else
        ρ_map_c = ρ_map .- c
        end
        psf_map!(H_c, psf, ρ_map_c, λ_map)
        L = Lkl(Diag(H_c) * F, d, w)
        return L(z)
    end
    
    #info, pc, fp = bobyqa!(likelihood!, psf_center, bnd_min, bnd_max, 
    #                        rho_tol[1], rho_tol[2]; kwds...) FIXME: use PowellMethods instead of PRIMA, check bounds issues
    #psf_center .= pc    
    psf_center .= bobyqa(likelihood!, psf_center, xl=bnd_min, xu=bnd_max, 
                            rhobeg=rho_tol[1], rhoend=rho_tol[2]; kwds...)[1]
    return psf_center
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

# See also
- [`InversePbm.get_grad_op`]
"""
function solve_analytic(F::SparseInterpolator{T},
                         H::AbstractArray{T,N},
                         D::CalibratedData,
                         Reg::Regularization) where {T,N}
    
    d, w = D.d, D.w
    sz=size(d)
    n=F.ncols

    @assert sz == size(w)
    @assert sz == size(H)
    @assert Reg.Reg.mu >=0
    
    A=spzeros(n,n)
    
    A = LinearInterpolators.SparseInterpolators.AtWA(F,H.*w.*H)
    A += InverseProblem.get_grad_op(Reg)*spdiagm(0 => ones(n))

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

# See also
- [`OptimPackNextGen.vmlbm!`]
"""
function solve_vmlmb(z0::AbstractVector{T},
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
    loss(D, H, F, z, Reg [, Bkg=undef])

yields the value of the loss function used to estimate the different parameters
of the problem:
`
    (Diag(H)*F*z + Bkg.b - D.d)'*Diag(D.w)*(Diag(H)*F*z + Bkg.b - D.d) + Reg(z) + regul(Bkg)
`
where `Bkg.b` and `regul(Bkg)` are not taken into account if `Bkg=undef`.

# See also
- [`AbstractBkg`](@ref)
- [`InverseProblem.Regul`]
"""
function loss(D::CalibratedData{T},
              psf::AbstractPSF,
              F::SparseInterpolator{T},
              z::AbstractVector{T},
              Reg::Regularization,
              Bkg::Union{AbstractBkg,UndefInitializer} = undef) where {T}
    H = similar(D.d)  
    psf_map!(H, psf, D.ρ_map, D.λ_map)

    # Computing Likelihood
    wks = vcopy(D.d)
    wks -= H .* (F*z)
    if Bkg != undef
        wks .-= Bkg.b
    end
    lkl = vdot(wks, D.w .* wks)

    #Computing Regularization 
    r = Reg(z)
    if typeof(psf) == NonParametricPSF
        r += psf.R(psf.h[:])
    end
    if Bkg != undef
        r += regul(Bkg)
    end

    return lkl + r
end

