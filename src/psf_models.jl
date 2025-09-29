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
        chromGaussianPSF(a) -> h
    
Yields an AbstractPSF `h`, more precisely a gaussian chromatic 
PSF parametred by `a`, such that the gaussian variance is  σ² = aλ².

Then: 

        h(ρ,λ) 
    
gives the value of the gaussian chromatic PSF for:

 - `ρ`, a float giving the position in the slit (in pixels)
    
 - `λ`, a float giving the wavelength (in μm)
    
which is:

        exp(-1/2 ρ²/σ²)/sqrt(2πσ²)
    
""" chromGaussianPSF

struct chromGaussianPSF <: ParametricPSF{1}
    a::Float64
end
 
function (P::chromGaussianPSF)(ρ::T,λ::T) where {T<:AbstractFloat}
    if iszero(λ)
        return  0.
    else
        σ² = P.a*λ*λ
        return exp(-1/2 *ρ*ρ/σ²)/sqrt(2*pi*σ²)
    end
end

@inline parameters(P::chromGaussianPSF) = getfield(P, :a)

function getfwhm(P::chromGaussianPSF, ρ::T,λ::T) where {T<:AbstractFloat}
    if iszero(λ)
        return  0.
    else
        σ² = P.a*λ*λ
        return 2*sqrt(2*log(2)*σ²)
    end
end


"""
        chromwmwGaussianPSF(p) -> h
    
Yields an AbstractPSF `h`, more precisely a gaussian chromatic 
PSF parametred by two variables given in a vector `p = [a,b]`:

 - `a` the chromatic scalling
   
 - `b` the minimum width
   
such that the the variance of the PSF is  σ² = aλ² + b

Then: 

        h(ρ,λ) 
    
gives the value of the gaussian chromatic PSF for:

 - `ρ`, a float giving the position in the slit (in pixels)
   
 - `λ`, a float giving the wavelength (in μm)
    
which is:

        exp(-1/2 ρ²/σ²)/sqrt(2πσ²)

""" chromwmwGaussianPSF

struct chromwmwGaussianPSF <: ParametricPSF{2}
    a::Float64
    b::Float64
    chromwmwGaussianPSF(p::AbstractArray{Float64,1}) = new(p[1], p[2])
end
    
function (P::chromwmwGaussianPSF)(ρ::T,λ::T) where {T<:AbstractFloat}
    if iszero(λ)
        return  0.
    else
        σ² = P.a*λ*λ + P.b
        return exp(-1/2 *ρ*ρ/σ²)/sqrt(2*pi*σ²)
    end
end

@inline parameters(P::chromwmwGaussianPSF) = (getfield(P, :a),getfield(P, :b))

function getfwhm(P::chromwmwGaussianPSF, ρ::T,λ::T) where {T<:AbstractFloat}
    if iszero(λ)
        return  0.
    else
        σ² = P.a*λ*λ +  P.b
        return 2*sqrt(2*log(2)*σ²)
    end
end


"""
        chromMoffatPSF(p) -> h
    
Yields an AbstractPSF `h`, more precisely a centrosymetric Moffat chromatic 
PSF parametred by two variables given in a vector `p=[a, β]`:

 - `a the chromatic scalling
   
 - `β` the Moffat parameter
   
such that the the variance of the PSF is  σ² = aλ²

Then: 

        h(ρ,λ) 
    
gives the value of the centrosymetric Moffat chromatic PSF for:

 - `ρ`, a float giving the position in the slit (in pixels)
   
 - `λ`, a float giving the wavelength (in μm)
    
which is:

        (1 + ρ²/σ²)^(-β)

""" chromMoffatPSF

struct chromMoffatPSF <: ParametricPSF{2}
    a::Float64
    β::Float64
    chromMoffatPSF(p::AbstractArray{Float64,1}) = new(p[1], p[2])
end

 
function (P::chromMoffatPSF)(ρ::T,λ::T) where {T<:AbstractFloat}
    if iszero(λ)
        return  0.
    else
        σ² = P.a*λ*λ
        A = sqrt(P.β/(σ²*pi*exp(1)))*(1 - 1/(2* P.β)) ^(1-P.β) #using Stirling Approximation
        return A * (1 + ρ*ρ/σ²)^(-P.β)
    end
end

@inline parameters(P::chromMoffatPSF) = (getfield(P, :a),getfield(P, :β))

function getfwhm(P::chromMoffatPSF, ρ::T,λ::T) where {T<:AbstractFloat}
    if iszero(λ)
        return  0.
    else
        σ² = P.a*λ*λ 
        return 2*sqrt(σ²*(2^(1/P.β) -1))
    end
end


"""
        chromwmwMoffatPSF(p) -> h
    
Yields an AbstractPSF `h`, more precisely a centrosymetric Moffat chromatic 
PSF parametred by two variables given in a vector `p=[a, b, β]`:

 - `a` the chromatic scalling
   
 - `b` the minimum width
 
 - `β` the Moffat parameter
   
such that the the variance of the PSF is  σ² = aλ² + b

Then: 

        h(ρ,λ) 
    
gives the value of the centrosymetric Moffat chromatic PSF for:

 - `ρ`, a float giving the position in the slit (in pixels)
   
 - `λ`, a float giving the wavelength (in μm)
    
which is:

        (1 + ρ²/σ²)^(-β)

""" chromwmwMoffatPSF

struct chromwmwMoffatPSF <: ParametricPSF{3}
    a::Float64
    b::Float64
    β::Float64
    chromwmwMoffatPSF(p::AbstractArray{Float64,1}) = new(p[1], p[2], p[3])
end

 
function (P::chromwmwMoffatPSF)(ρ::T,λ::T) where {T<:AbstractFloat}
    if iszero(λ)
        return  0.
    else
        σ² = P.a*λ*λ +P.b 
        A = sqrt(P.β/(σ²*pi*exp(1)))*(1 - 1/(2* P.β)) ^(1-P.β) #using Stirling Approximation
        return A * (1 + ρ*ρ/σ²)^(-P.β)
    end
end

@inline parameters(P::chromwmwMoffatPSF) = (getfield(P, :a),
                                            getfield(P, :b),
                                            getfield(P, :β))


function getfwhm(P::chromwmwMoffatPSF, ρ::T,λ::T) where {T<:AbstractFloat}
    if iszero(λ)
        return  0.
    else
        σ² = P.a*λ*λ +  P.b
        return 2*sqrt(σ²*(2^(1/P.β) -1))
    end
end


"""
        SeriesExpansionPSF(x) -> h
       
""" SeriesExpansionPSF

struct SeriesExpansionPSF{T<:AbstractFloat,
                          V<:AbstractArray{Float64,2},
                          S<:AbstractRange{Float64},
                          K<:Kernel,
                          R<:Regularization} <: NonParametricPSF{3}
    h::V
    x::S
    a::T
    ker::K
    R::R
end

function SeriesExpansionPSF(h::AbstractVector{AbstractArray{T,N}};
                              x::AbstractRange{T},
                              a=0.,
                              ker::Kernel = CatmullRomSpline(T, Flat),
                              R::Regularization = tikhonov()) where {T<: AbstractFloat,N}
        @assert length(x) == size(h)[2]                      
        return SeriesExpansionPSF(h,x, a, ker, R)
end

function (P::SeriesExpansionPSF)(ρ::AbstractArray{T,N},
                                λ::AbstractArray{T,N};
                                λref = maximum(λ)) where {T<:AbstractFloat,N}
    order=size(P.h)[2]
    γ = λref ./ λ
    γ[λ .== 0] .= 0.
    return ChromaticSeriesExpansionsInterpolator(P.ker, ρ, λ, P.x, order=order, λref=λref, a=P.a)
end

@inline parameters(P::SeriesExpansionPSF) = (getfield(P, :h))
@inline shift(P::SeriesExpansionPSF)=getfield(P,:a)
@inline coord(P::SeriesExpansionPSF)=get_field(P,:x)

function getfwhm(P::SeriesExpansionPSF, ρ::T,λ::T) where {T<:AbstractFloat}
    # @error "Not implented yet"
    σ²=0.017*λ*λ #FIXME:use something less arbitrary
    return 2*sqrt(2*log(2)*σ²)#one(T) #FIXME: default value adjusted by mask_width
end


