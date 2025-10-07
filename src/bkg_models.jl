#
# bkg_models.jl
# 
# The users of ASSET can use this file as a template to implement their 
# own background model.
#
# ------------------------------------------------
#
# This file is part of ASSET


"""
    BkgMdl(b, R)

yields a structure of type `AbstractBkg`, composed of an `AbstractArray` `b` and
an associated regularization `R` of type `Regularization`.

# Example
```julia
# create a BkgMdl structure
julia> B = BkgMdl(b, R)
# add the structure to an AbstractArray
julia> A = ones(size(b))
julia> A + B
```
To apply the regularization to the background array:
```julia
julia> regul(B)
```
To fit the background to some `CalibratedData` `D`:
```julia
julia> fit_bkg!(B, D, true, kwds...)
```
where `kwds` specifies all the keywords used by the optimization method to
estimate the background.

# See also
- [`AbstractBkg`](@ref)
"""
struct BkgMdl{T,N,K} <: AbstractBkg where {T,N,K}
    b::AbstractArray{T,N}
    R::Regularization
    kwds::K
end

BkgMdl(b::AbstractArray{T,N}, R::Regularization) where {T,N} = BkgMdl(b, R, (lower=T(0),))

get_bkg(B::BkgMdl) = B.b

regul(B::BkgMdl) = B.R(B.b)

function fit_bkg!(B::BkgMdl{T}, 
                  D::CalibratedData) where {T <: AbstractFloat}
    
    function fg_solve!(x, g)
        fill!(g,0.0);
        f = T(0)
        if  InverseProblem.multiplier(B.R)>0 
            f += B.R(x,g)
        end
        r = (x .- D.d)
        wr = D.w .*r
        vupdate!(g, 1, sum(wr, dims=3)[:,:])
        f +=  vdot(r,wr)/2
        return f 
    end
    
    vmlmb!(fg_solve!, B.b; B.kwds...)
    
    return B.b
end

