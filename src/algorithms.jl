

"""
    extract_spectrum!()

used to estimate the spectrum of an object observed with long-slit spectroscopy,
when the data is corrupted by a background component and noise. The direct model
of such data can be writen:
`
                    d = Bkg + H*(F.z) + n
`
where `*` denotes the element-wise multiplication and `.` the matrix product.
`z` is the spectrum of the object of interest and `F` its associated
interpolation operator, `H` the operator modeling the Point Spread Function of
the instrument, while `Bkg` is the background component to disentangle from the
object of interest and `n` accounts for noises.

FIXME: fit_bkg! must take in arg `AbstractBkg`, `AbstractArray{T,N}`,
`AbstractArray{T,N}`, `Regularization`, `AbstractArray{T,N}`,
`AbstractArray{T,N}`.

#FIXME: size(psf_center) = (1, 1, size(ρ_map,3))
#FIXME: keywords for psf_parameters bound

"""
function extract_spectrum!(z::AbstractVector{T},
    F::SparseInterpolator{T},
    psf_center::AbstractArray{T},
    psf::AbstractPSF,
    ρ_map::AbstractArray{T,N},
    λ_map::AbstractArray{T,N},
    d::AbstractArray{T,N},
    w::AbstractArray{T,N},
    Reg::Regularization,
    Bkg::Union{AbstractBkg,UndefInitializer} = undef,
    Reg_bkg::Union{Regul,UndefInitializer} = undef,
    fit_bkg!::Union{Function,UndefInitializer} = undef; #FIXME: not as an argument
    auto_calib::Val = Val(true),
    mask_width::Real = 3.0,
    max_iter::Int = 1000,
    loss_tol::Tuple{Real,Real} = (0.0, 1e-6),
    z_tol::Tuple{Real,Real} = (0.0, 1e-6),
    kwds...) where {T,N}
    
    @assert axes(d) == axes(w) == axes(ρ_map) == axes(λ_map)
    @assert axes(d) == LinearInterpolators.output_size(F)
    
    # Initialization
    iter = 0
    loss_last = loss(d, w, H, F, s, Reg, Bkg)
    z_last = copy(z)
    res = copy(d)
    ρ_map_centered = ρ_map .- reshape(psf_center, 1, 1, length(psf_center))
    H = similar(d)
    while true
        if fit_bkg! != undef #FIXME: if Bkg != undef
            # Mask object
            mask_z = mask_object(d, ρ_map_centered, λ_map; mask_width=mask_width)
            # Estimate background
            bkg_step!(Bkg, fit_bkg!, res, w, Reg_bkg, ρ_map, λ_map, mask_z; 
                      hide_object=(iter == 0)) #FIXME: fit_bkg! should not be an argument
            copyto!(res, d - Bkg)
        end
        # Estimate object's spectrum
        object_step!(z, psf, psf_center, F, res, w, Reg, ρ_map, λ_map; 
                     auto_calib=auto_calib, kwds...) #FIXME: keywords for psf_parameters bound
        copyto!(ρ_map_centered, ρ_map .- reshape(psf_center, 1, 1, length(psf_center)))
        psf_map!(H, psf, ρ_map_centered, λ_map)
        # Stop criterions
        loss = loss(d, w, H, F, z, Reg, Bkg)
        if (iter >= max_iter) || test_tol(loss, loss_last, loss_tol) || 
                                 test_tol(z, z_last, z_tol)
            if auto_calib == Val(:delay)
                auto_calib = Val(true)
            else
                break
            end
        end
        copyto!(res, d - H .* (F*z))
        iter += 1
        loss_last = lost
        copyto!(z_last, z)
    end

    return z
end




"""
"""
function mask_object(d::AbstractArray{T,N},
    ρ_map::AbstractArray{T,N},
    λ_map::AbstractArray{T,N};
    mask_width::Real = 3.0) where {T,N}
    
    @assert axes(d) == axes(ρ_map) == axes(λ_map)

    mask_z = ones(T, size(d))
    for f in 1:size(d, 3)
        loc_z = lambda_ref(λ_map[:,:,f]) ./ λ_map[:,:,f] .* abs.(ρ_map[:,:,f])
        frame_mask = mask_z[:,:,f]
        frame_mask[loc_z .< mask_width] .= T(0)
    end

    return mask_z
end




"""
"""
function bkg_step!(Bkg::AbstractBkg,
    fit_bkg!::Function,
    d::AbstractArray{T,N},
    w::AbstractArray{T,N},
    Reg_bkg::Regularization,
    ρ_map::AbstractArray{T,N},
    λ_map::AbstractArray{T,N},
    mask_z::AbstractArray{T,N} = ones(size(w));
    hide_object::Bool = false) where {T,N}

    # Hide object component for first estimation of background
    w_bkg = similar(w)
    if hide_object
        @inbounds for i in eachindex(w_bkg, mask_z, w)
            w_bkg[i] = mask_z[i]*w[i]
        end
    else
        copyto!(w_bkg, w)
    end
    # Fit background
    return fit_bkg!(Bkg, d, w_bkg, Reg_bkg, ρ_map, λ_map)
end




"""
"""
function object_step!(z::AbstractVector{T},
    psf::AbstractPSF,
    psf_center::AbstractArray{T},
    F::SparseInterpolator{T},
    d::AbstractArray{T,N},
    w::AbstractArray{T,N},
    Reg::Regularization,
    ρ_map::AbstractArray{T,N},
    λ_map::AbstractArray{T,N};
    auto_calib::Val = Val(true),
    psf_params_bnds::AbstractVector{Tuple{T}} = [(0,0) for i in 1:length(psf[:])],
    psf_center_bnds::AbstractVector{Tuple{T}} = [(0,0) for i in 1:length(psf_center)],
    max_iter::Int = 1000,
    loss_tol::Tuple{Real,Real} = (0.0, 1e-6),
    z_tol::Tuple{Real,Real} = (0.0, 1e-6),
    rho_params::Tuple{Real,Real} = (2.5e-1, 1e-5),
    rho_center::Tuple{Real,Real} = (5e-1, 1e-5),
    kwds...) where {T,N}
    
    # Initialization
    iter = 0
    loss_last = loss(d, w, H, F, z, Reg)
    z_last = copy(z)
    ρ_map_centered = ρ_map .- reshape(psf_center, 1, 1, length(psf_center))
    H = similar(d)
    psf_map!(H, psf, ρ_map_centered, λ_map)
    while true
        # Extract spectrum
        fit_spectrum!(z, F, H, d, w, Reg; kwds...)
        # Stop criterions
        (auto_calib == Val(true)) && (loss = loss(d, w, H, F, z, Reg))
        if (auto_calib != Val(true)) || (iter >= max_iter) || 
                                        test_tol(loss, loss_last, loss_tol) || 
                                        test_tol(z, z_last, z_tol)
            break
        end
        # Auto-calibration step
        if auto_calib == Val(true)
            fit_psf_params!(psf, z, F, res, w, ρ_map_centered, λ_map; 
                            psf_params_bnds=psf_params_bnds, rho_tol=rho_params)
            fit_psf_center!(psf_center, psf, z, F, res, w, ρ_map, λ_map;
                            psf_center_bnds=psf_center_bnds, rho_tol=rho_center)
            copyto!(ρ_map_centered, ρ_map .- reshape(psf_center, 1, 1, length(psf_center)))
            psf_map!(H, psf, ρ_map_centered, λ_map)
        end
    end

    return z, psf, psf_center
end




"""

#FIXME: psf_center is not an argument
"""
function psf_map!(map::AbstractArray{T,N},
                 h::AbstractPSF,
                 ρ::AbstractArray{T,N},
                 λ::AbstractArray{T,N}) where {N, T<:AbstractFloat}
    sz=size(map)
    @assert sz == size(ρ)
    @assert sz == size(λ)
    
    @inbounds for i in eachinde(map, ρ, λ)
        map[i] = h(ρ[i], λ[i])
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

"""

FIXME: must see if direct inversion or vmlmb

"""
function fit_spectrum!(z::AbstractVector{T},
    F::SparseInterpolator{T},
    H::AbstractArray{T,N},
    d::AbstractArray{T,N},
    w::AbstractArray{T,N},
    Reg::Regularization;
    kwds...) where {T,N}

    if Reg.direct_inversion
        copyto!(z, solve_analytic!(F, H, d, w, Reg))
    else
        copyto!(z, solve_vmlmb!(z, Diag(H) * F, d, w, Reg; kwds...))
    end
    
    return z
end




"""
"""
function solve_analytic!(F::SparseInterpolator{T},
    H::AbstractArray{T,N},
    d::AbstractArray{T,N},
    w::AbstractArray{T,N},
    Reg::Regularization) where {T,N}

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
"""
function solve_vmlmb!(z0::AbstractVector{T},
    A::Mapping,
    d::AbstractArray{T,N},
    w::AbstractArray{T,N},
    Reg::Regularization;
    nonnegative::Bool = false,
    kwds...) where {T,N}

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
"""
function fit_psf_params!(psf::AbstractPSF,
    z::AbstractVector{T},
    F::SparseInterpolator{T},
    d::AbstractArray{T,N},
    w::AbstractArray{T,N},
    ρ_map::AbstractArray{T,N},
    λ_map::AbstractArray{T,N};
    psf_params_bnds::AbstractVector{Tuple{T}} = [(0,0) for i in 1:length(psf[:])],
    rho_tol::Tuple{Real,Real} = (2.5e-1, 1e-5),
    kwds...) where {T,N}
    
    @assert length(psf[:]) == length(psf_params_bnds)

    bnd_min = zeros(length(psf_params_bnds))
    bnd_max = zeros(length(psf_params_bnds))
    for p in eachindex(psf_params_bnds)
        bnds_param = psf_params_bnds[p]
        bnd_min[p] = bnds_param[1]
        bnd_max[p] = bnds_param[2]
    end

    H_p = zeros(T, size(d))
    function likelihood!(p)
        psf_p = psf(p)
        psf_map!(H_p, psf, ρ_map, λ_map)
        L = Lkl(Diag(H_p) * F, d, w)

        return L(z)
    end

    return Bobyqa.minimize!(p -> likelihood!(p), psf[:], bnd_min, bnd_max, 
                            rho_tol[1], rho_tol[2]; kwds...)[2]
end


function fit_psf_center!(psf_center::AbstractVector{T},
    psf::AbstractPSF,
    z::AbstractVector{T},
    F::SparseInterpolator{T},
    d::AbstractArray{T,N},
    w::AbstractArray{T,N},
    ρ_map::AbstractArray{T,N},
    λ_map::AbstractArray{T,N};
    psf_center_bnds::AbstractVector{Tuple{T}} = [(0,0) for i in 1:length(psf_center)],
    rho_tol::Tuple{Real,Real} = (5e-1, 1e-5),
    kwds...) where {T,N}
    
    @assert length(psf_center) == length(psf_center_bnds)

    bnd_min = zeros(length(psf_center))
    bnd_max = zeros(length(psf_center))
    for p in eachindex(psf_center_bnds)
        bnd = psf_center_bnds[p]
        bnd_min[p] = bnd[1]
        bnd_max[p] = bnd[2]
    end

    ρ_map_c = zeros(T, size(d))
    H_c = zeros(T, size(d))
    function likelihood!(c)
        ρ_map_c = ρ_map .- reshape(psf_center, 1, 1, length(psf_center))
        psf_map!(H_c, psf, ρ_map_c, λ_map)
        L = Lkl(Diag(H_c) * F, d, w)

        return L(z)
    end


    return Bobyqa.minimize!(c -> likelihood!(c), psf_center, bnd_min, bnd_max, 
                            rho_tol[1], rho_tol[2]; kwds...)[2]
end



"""

FIXME: Bkg must contain the regul or result of regul for Bkg. And overload call
and call! with AbstractBkg?

"""
function loss(d::AbstractArray{T,N},
    w::AbstractArray{T,N},
    H::AbstractArray{T,N},
    F::SparseInterpolator{T},
    z::AbstractVector{T},
    Reg::Regularization,
    Bkg::Union{AbstractBkg,UndefInitializer} = undef) where {T,N}

    # Computing Likelihood
    wks = vcopy(d)
    if Bkg != undef
        @. wks -= H .* (F*z) + Bkg
    else
        @. wks -= H .* (F*z)
    end
    lkl = vdot(wks, w .* wks)

    #Computing Regularization
    r = Reg(z)
    if Bkg != undef
        r += regul(Bkg)
    end

    return lkl + r
end
