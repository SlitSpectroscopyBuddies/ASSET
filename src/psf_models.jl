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
        h=chromGaussianPSF(a)
    
gives an AbstractPSF, more precisely a gaussian chromatic 
PSF parametred by a, such that the gaussian variance is  σ² = aλ².

Then: 

        h(ρ,λ) 
    
gives the value of the gaussian chromatic PSF for:

 - ρ, a float giving the position in the slit (in pixels)
    
 - λ, a float giving the wavelength (in μm)
    
which is:

        exp(-1/2 ρ²/σ²)/sqrt(2πσ²)
    
""" chromGaussianPSF

struct chromGaussianPSF <: AbstractPSF{1}
    a::Float64
end
 
function (P::chromGaussianPSF)(ρ::T,λ::T) where {T<:AbstractFloat}
    σ² = P.a*λ*λ
    return exp(-1/2 *ρ*ρ/σ²)/sqrt(2*pi*σ²)
end

@inline parameters(P::chromGaussianPSF) = getfield(P, :a)

function getfwhm(P::chromGaussianPSF, ρ::T,λ::T) where {T<:AbstractFloat}
    σ² = P.a*λ*λ
    return 2*sqrt(2*log(2)*σ²)
end


"""
        h=chromwmwGaussianPSF(p)
    
gives an AbstractPSF, more precisely a gaussian chromatic 
PSF parametred by two variables given in a vector:

 - a the chromatic scalling
   
 - b the minimum width
   
such that the the variance of the PSF is  σ² = aλ² + b

Then: 

        h(ρ,λ) 
    
gives the value of the gaussian chromatic PSF for:

 - ρ, a float giving the position in the slit (in pixels)
   
 - λ, a float giving the wavelength (in μm)
    
which is:

        exp(-1/2 ρ²/σ²)/sqrt(2πσ²)

""" chromwmwGaussianPSF

struct chromwmwGaussianPSF <: AbstractPSF{1}
    a::Float64
    b::Float64
    chromwmwGaussianPSF(p::AbstractArray{Float64,1}) = new(p[1], p[2])
end
    
function (P::chromwmwGaussianPSF)(ρ::T,λ::T) where {T<:AbstractFloat}
    σ² = P.a*λ*λ + P.b
    return exp(-1/2 *ρ*ρ/σ²)/sqrt(2*pi*σ²)
end

@inline parameters(P::chromwmwGaussianPSF) = (getfield(P, :a),getfield(P, :b))

function getfwhm(P::chromwmwGaussianPSF, ρ::T,λ::T) where {T<:AbstractFloat}
    σ² = P.a*λ*λ +  P.b
    return 2*sqrt(2*log(2)*σ²)
end


"""
        h=chromMoffatPSF(p)
    
gives an AbstractPSF, more precisely a centrosymetric Moffat chromatic 
PSF parametred by two variables given in a vector:

 - a the chromatic scalling
   
 - β the Moffat parameter
   
such that the the variance of the PSF is  σ² = aλ²

Then: 

        h(ρ,λ) 
    
gives the value of the centrosymetric Moffat chromatic PSF for:

 - ρ, a float giving the position in the slit (in pixels)
   
 - λ, a float giving the wavelength (in μm)
    
which is:

        (1 + ρ²/σ²)^(-β)

""" chromMoffatPSF

struct chromMoffatPSF <: AbstractPSF{1}
    a::Float64
    β::Float64
    chromMoffatPSF(p::AbstractArray{Float64,1}) = new(p[1], p[2])
end

 
function (P::chromMoffatPSF)(ρ::T,λ::T) where {T<:AbstractFloat}
    σ² = P.a*λ*λ
    return (P.β-1)/(pi*σ²) * (1 + ρ*ρ/σ²)^(-P.β)
end

@inline parameters(P::chromMoffatPSF) = (getfield(P, :a),getfield(P, :β))

function getfwhm(P::chromMoffatPSF, ρ::T,λ::T) where {T<:AbstractFloat}
    σ² = P.a*λ*λ 
    return 2*sqrt(σ²*(2^(1/P.β) -1))
end


"""
        h=chromwmwMoffatPSF(p)
    
gives an AbstractPSF, more precisely a centrosymetric Moffat chromatic 
PSF parametred by two variables given in a vector:

 - a the chromatic scalling
   
 - b the minimum width
 
 - β the Moffat parameter
   
such that the the variance of the PSF is  σ² = aλ² + b

Then: 

        h(ρ,λ) 
    
gives the value of the centrosymetric Moffat chromatic PSF for:

 - ρ, a float giving the position in the slit (in pixels)
   
 - λ, a float giving the wavelength (in μm)
    
which is:

        (1 + ρ²/σ²)^(-β)

""" chromwmwMoffatPSF

struct chromwmwMoffatPSF <: AbstractPSF{1}
    a::Float64
    b::Float64
    β::Float64
    chromwmwMoffatPSF(p::AbstractArray{Float64,1}) = new(p[1], p[2], p[3])
end

 
function (P::chromwmwMoffatPSF)(ρ::T,λ::T) where {T<:AbstractFloat}
    σ² = P.a*λ*λ +P.b
    return (P.β-1)/(pi*σ²) * (1 + ρ*ρ/σ²)^(-P.β)
end

@inline parameters(P::chromwmwMoffatPSF) = (getfield(P, :a),
                                            getfield(P, :b),
                                            getfield(P, :β))


function getfwhm(P::chromwmwMoffatPSF, ρ::T,λ::T) where {T<:AbstractFloat}
    σ² = P.a*λ*λ +  P.b
    return 2*sqrt(σ²*(2^(1/P.β) -1))
end


