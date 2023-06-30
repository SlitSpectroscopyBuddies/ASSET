

"""

"""
abstract type AbstractBkg end




"""

#FIXME: overload regul(bkg)
"""
struct BkgMdl{T,N} <: AbstractBkg where {T,N}
    b::AbstractArray{T,N}
    R::Regularization
end

Base.:-(m::AbstractArray{T,N}, B::BkgMdl{T,N}) where {T,N} = m - B.b

regul(B::BkgMdl) = B.R(B.b)




"""

#FIXME: overload regul(bkg)

"""
struct ParametrizedBkgMdl{T} <: AbstractBkg where {T}
    θ::AbstractVector{T}
    f::Function
    R::Regularization
end

Base.:-(m::AbstractArray{T,N}, B::ParametrizedBkgMdl{T}) where {T,N} = m - B.f(B.θ)

regul(B::ParametrizedBkgMdl) = B.R(B.f(B.θ))

