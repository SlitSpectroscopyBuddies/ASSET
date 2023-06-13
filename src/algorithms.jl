

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
    maxiter::Int = 1000) where {T,M}
    
    @assert axes(d) == axes(w)
    
    # Mask object
    mask_z = ones(T, size(d))
    #FIXME: pos_z = 

    iter = 0
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
        fit_spectrum!(z, F, H, res, w, Reg, ρ_map, λ_map)
        # Stop criterions
        #FIXME: loss = 
        if (iter >= maxiter)
            if auto_calib == Val(:delay)
                auto_calib = Val(true)
            else
                break
            end
        end
        # Auto-calibration step
        if auto_calib == Val(true)
            fit_psf_params!(psf_center, psf_params_bnds, psf, z, res, w, ρ_map, λ_map;
                            psf_center_bnds=psf_center_bnds)
            psf_map!(H, psf_center, psf, ρ_map, λ_map)
        end
        copyto!(res, d - H .* (F*z))
        iter += 1
    end
    return z
end


function psf_map!()

end

"""

FIXME: must see if direct inversion or vmlmb

"""
function fit_spectrum!()

end

function fit_psf_params!(;
    psf_center_bnds::AbstractVector{Tuple{T}} = [(0,0) for i in 1:length(psf_center)],)

end

"""

FIXME: bkg must contain the regul or result of regul for bkg. And overload call
and call! with AbstractBkg

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
        r += regul(bkg)(bkg)
    end

    return lkl + r
end