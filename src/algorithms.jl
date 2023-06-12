

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

"""
function extract_spectrum!(z::AbstractVector{T},
    psf::AbstractPSF,
    bkg::BkgModel,
    map_rho::AbstractArray{T,N},
    map_lambda::AbstractArray{T,N},
    d::AbstractArray{T,M},
    w::AbstractArray{T,M},
    Reg::Regul,
    solve_bkg!::Union{Function,UndefInitializer} = undef;
    auto_calib::Val = Val(true),
    maxiter::Int = 1000) where {T,M,N}
    
    @assert axes(d) == axes(w)
    
    # Mask object
    mask_z = ones(T, size(d))
    #FIXME: pos_z = 

    iter = 0
    res = similar(d)
    w_prime = similar(w)
    H = similar(d)
    while true
        if solve_bkg! != undef
            # Hide object component for first estimation of background
            if iter == 0
                @inbounds for i in eachindex(w_prime, mask_z, w)
                    w_prime[i] = mask_z[i]*w[i]
                end
            elseif iter == 1
                copyto!(w_prime, w)
            end
            # Estimate background
            solve_bkg!(bkg, res, w_prime, map_rho, map_lambda)
        end
        # Estimate object's spectrum
        copyto!(res, d - bkg)
        psf_map!(H, psf, map_rho, map_lambda)
        fit_spectrum!(z, H, res, w, Reg, map_rho, map_lambda)
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
            fit_psf_params!(psf, z, res, w, map_rho, map_lambda)
            psf_map!(H, psf, map_rho, map_lambda)
        end
        copyto!(res, d - mdl_obj(H, z, map_lambda))
        iter += 1
    end
    return z
end


function psf_map!()

end

function fit_spectrum!()

end

function fit_psf_params!()

end

function mdl_obj()

end
