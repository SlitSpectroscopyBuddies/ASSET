

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

function fit_bkg!(B::BkgMdl, 
                  D::CalibratedData; 
                  nonnegative = true,
                  kwds...)
    function fg_solve!(x, g) where {T <: AbstractFloat}
        fill!(g,0.0);
        if μ>0
            f += R(x,g)
        end
        r = (x .+ M - d)
        wr = w .*r
        vupdate!(g, 1, sum(wr, dims=3)[:,:])
        f =  vdot(r,wr)/2
        return f 
    end
    
    if nonnegative
        vmlmb!(fg_solve!, B.b; lower=T(0), kwds...)
    else
        vmlmb!(fg_solve!, B.b; kwds...)
    end
    
    return B.b
end


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

