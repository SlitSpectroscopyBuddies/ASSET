

"""

"""
abstract type AbstractBkg end




"""

#FIXME: overload regul(bkg)
"""
struct BkgMdl <: AbstractBkg where {T,N}
    b::AbstractArray{T,N}
    R::Regularization
end

Base.:+(m::AbstractArray{T,M}, B::ParametrizedBkgMdl) = m + B.b

regul(B::BkgMdl) = B.R(B.b)




"""

#FIXME: overload regul(bkg)

"""
struct ParametrizedBkgMdl <: AbstractBkg where {T}
    θ::AbstractVector{T}
    f::Function
    R::Regularization
end

Base.:+(m::AbstractArray{T,M}, B::ParametrizedBkgMdl) = m + B.f(B.θ)

regul(B::ParametrizedBkgMdl) = B.R(B.f(B.θ))

