

"""
    extract_spectrum!()

used to estimate the spectrum of an object observed with long-slit spectroscopy,
when the data is corrupted by a background component and noise. The direct model
of such data can be writen:
`
                    d = bkg + H*(F.z) + n
`
where `*` denotes the element-wise multiplication and `.` the matrix product.
`z` is the spectrum of the object of interest and `F` its associated
interpolation operator, `H` the operator modeling the Point Spread Function of
the instrument, while `bkg` is the background component to disentangle from the
object of interest and `n` accounts for noises.

FIXME: fit_bkg! must take in arg `AbstractBkg`, `AbstractArray{T,M}`,
`AbstractArray{T,M}`, `Regularization`, `AbstractArray{T,M}`, `AbstractArray{T,M}`.


"""
function extract_spectrum!(z::AbstractVector{T},
    F::SparseInterpolator{T},
    psf_center::AbstractVector{T},
    psf::AbstractPSF,
    ρ_map::AbstractArray{T,M},
    λ_map::AbstractArray{T,M},
    d::AbstractArray{T,M},
    w::AbstractArray{T,M},
    Reg::Regularization,
    bkg::Union{AbstractBkg,UndefInitializer} = undef,
    Reg_bkg::Union{Regul,UndefInitializer} = undef,
    fit_bkg!::Union{Function,UndefInitializer} = undef;
    psf_center_bnds::AbstractVector{Tuple{T}} = [(0,0) for i in 1:length(psf_center)],
    auto_calib::Val = Val(true),
    mask_width::Real = 3.0,
    maxiter::Int = 1000,
    loss_tol::Tuple{Real,Real} = (0.0, 1e-6),
    z_tol::Tuple{Real,Real} = (0.0, 1e-6),
    rho_params::Tuple{Real,Real} = (2.5e-1, 1e-5),
    rho_center::Tuple{Real,Real} = (5e-1, 1e-5),
    kwds...) where {T,M}
    
    #FIXME: fill tests
    @assert axes(d) == axes(w)
    
    # Mask object
    mask_z = mask_object(psf_center, d, ρ_map, λ_map; mask_width=mask_width)

    # Initialization
    iter = 0
    loss_last = loss(d, w, H, F, s, Reg, bkg)
    z_last = copy(z)
    res = similar(d)
    w_bkg = similar(w)
    H = similar(d)
    while true
        if fit_bkg! != undef
            # Hide object component for first estimation of background
            if iter == 0
                @inbounds for i in eachindex(w_bkg, mask_z, w)
                    w_bkg[i] = mask_z[i]*w[i]
                end
            elseif iter == 1
                copyto!(w_bkg, w)
            end
            # Estimate background
            fit_bkg!(bkg, res, w_bkg, Reg_bkg, ρ_map, λ_map)
            copyto!(res, d - bkg)
        end
        # Estimate object's spectrum
        psf_map!(H, psf_center, psf, ρ_map, λ_map)
        fit_spectrum!(z, F, H, res, w, Reg; kwds...)
        # Stop criterions
        loss = loss(d, w, H, F, z, Reg, bkg)
        if (iter >= maxiter) || test_tol(loss, loss_last, loss_tol) || test_tol(z, z_last, z_tol)
            if auto_calib == Val(:delay)
                auto_calib = Val(true)
            else
                break
            end
        end
        # Auto-calibration step
        if auto_calib == Val(true)
            fit_psf_params!(psf_center, psf, z, F, res, w, ρ_map, λ_map; rho_params)
            fit_psf_center!(psf_center, psf, z, F, res, w, ρ_map, λ_map;
                            psf_center_bnds=psf_center_bnds, rho_center)
            psf_map!(H, psf_center, psf, ρ_map, λ_map)
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
function mask_object(psf_center::AbstractVector{T},
    d::AbstractArray{T,N},
    ρ_map::AbstractArray{T,N},
    λ_map::AbstractArray{T,N};
    mask_width::Real = 3.0) where {T,N}

    mask_z = ones(T, size(d))
    for f in 1:size(d, 3)
        loc_z = lambda_ref(λ_map[:,:,f]) ./ λ_map[:,:,f] .* abs.(ρ_map[:,:,f]
                                                                 .- psf_center[f])
        frame_mask = mask_z[:,:,f]
        frame_mask[loc_z .< mask_width] .= T(0)
    end

    return mask_z
end




function psf_map!()

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
    Reg::HomogeneousRegularization;
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
function fit_psf_params!(psf_center::AbstractVector{T},
    psf::AbstractPSF,
    z::AbstractVector{T},
    F::SparseInterpolator{T},
    d::AbstractArray{T,N},
    w::AbstractArray{T,N},
    ρ_map::AbstractArray{T,N},
    λ_map::AbstractArray{T,N};
    rho_tol::Tuple{Real,Real} = (2.5e-1, 1e-5),
    kwds...) where {T,N}
    
    params = parameters(psf)
    bnds = get_param_bnds(psf)
    bnd_min = zeros(length(bnds))
    bnd_max = zeros(length(bnds))
    for p in eachindex(bnds)
        bnds_param = bnds[p]
        bnd_min[p] = bnds_param[1]
        bnd_max[p] = bnds_param[2]
    end

    H_p = zeros(T, size(d))
    function likelihood!(p)
        psf_p = psf(p)
        psf_map!(H_p, psf_center, psf, ρ_map, λ_map)
        L = Lkl(Diag(H_p) * F, d, w)

        return L(z)
    end

    return Bobyqa.minimize!(p -> likelihood!(p), params, bnd_min, bnd_max, 
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
        for m in 1:size(ρ_map_c, 3)
            ρ_map_c[:,:,m] = ρ_map[:,:,m] .- psf_center[m]
        end
        psf_map!(H_c, psf_center, psf, ρ_map_c, λ_map)
        L = Lkl(Diag(H_c) * F, d, w)

        return L(z)
    end


    return Bobyqa.minimize!(c -> likelihood!(c), psf_center, bnd_min, bnd_max, 
                            rho_tol[1], rho_tol[2]; kwds...)[2]
end



"""

FIXME: bkg must contain the regul or result of regul for bkg. And overload call
and call! with AbstractBkg?

"""
function loss(d::AbstractArray{T,M},
    w::AbstractArray{T,M},
    H::AbstractArray{T,M},
    F::SparseInterpolator{T},
    z::AbstractVector{T},
    Reg::Regularization,
    bkg::Union{AbstractBkg,UndefInitializer} = undef) where {T,M}

    # Computing Likelihood
    wks = vcopy(d)
    if bkg != undef
        @. wks -= H .* (F*z) + bkg
    else
        @. wks -= H .* (F*z)
    end
    lkl = vdot(wks, w .* wks)

    #Computing Regularization
    r = Reg(z)
    if bkg != undef
        r += regul(bkg)
    end

    return lkl + r
end
