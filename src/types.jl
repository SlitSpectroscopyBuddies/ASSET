

"""

"""
abstract type AbstractBkg end

regul(B::AbstractBkg) = B.R



"""

#FIXME: overload regul(bkg)
#FIXME: needs to work for multiple regul at a time

"""
struct ParametrizedBkgMdl <: AbstractBkg where {T}
    θ::AbstractVector{T}
    f::Function
    R::Regularization
end

Base.:+(m::AbstractArray{T,M}, B::ParametrizedBkgMdl) = m + f(B.θ)



"""

#FIXME: overload regul(bkg)
"""
struct BkgMdl <: AbstractBkg where {T,N}
    b::AbstractArray{T,N}
    R::Regularization
end



