#
# psf_models.jl
#
# Provide different types of parametric psf using PointSpreadFunction:
#
# https://github.com/emmt/PointSpreadFunctions.jl.git
# 
# The users of ASSET can use this file as a template to implement their 
# own psf model.
#
# ------------------------------------------------
#
# This file is part of ASSET



"""
        h=chromPSF(α)
    
gives an AbstractPSF, more precisely a gaussian chromatic 
PSF parametred by α, such that the gaussian variance is  σ² = αλ².

Then: 

        h(ρ,λ) 
    
gives the value of the gaussian chromatic PSF for:

 - ρ, a float giving the position in the slit (in pixels)
    
 - λ, a float giving the wavelength (in μm)
    
which is:

        exp(-1/2 ρ²/σ²)/sqrt(2πσ²)
    
""" chromPSF

struct chromPSF <: AbstractPSF{1}
    α::Float64
    chromPSF(α::Float64) = new(α)
end
 
function (P::chromPSF)(ρ::T,λ::T) where {T<:AbstractFloat}
    σ² = P.α*λ*λ
    return exp(-1/2 *ρ*ρ/σ²)/sqrt(2*pi*σ²)
end



"""
        h=chromwmwPSF(p)
    
gives an AbstractPSF, more precisely a gaussian chromatic 
PSF parametred by two variables given in a vector:

 - α the chromatic scalling
   
 - β the minimum width
   
such that the the variance of the PSF is  σ² = αλ² + β

Then: 

        h(ρ,λ) 
    
gives the value of the gaussian chromatic PSF for:

 - ρ, a float giving the position in the slit (in pixels)
   
 - λ, a float giving the wavelength (in μm)
    
which is:

        exp(-1/2 ρ²/σ²)/sqrt(2πσ²)

""" chromwmwPSF

struct chromwmwPSF <: AbstractPSF{1}
    α::Float64
    β::Float64
    chromwmwPSF(p::AbstractArray{Float64,1}) = new(p[1], p[2])
end
    
function (P::chromwmwPSF)(ρ::T,λ::T) where {T<:AbstractFloat}
    σ² = P.α*λ*λ + P.β
    return exp(-1/2 *ρ*ρ/σ²)/sqrt(2*pi*σ²)
end



