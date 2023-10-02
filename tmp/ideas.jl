

abstract type AbstractBkg end

#TODO: so it applies to every possible structure of type AbstractBkg (see def of
#get_bkg method)
Base.:+(m::AbstractArray{T,N}, B::AbstractBkg) where {T,N} = m .+ get_bkg(B)
Base.:-(m::AbstractArray{T,N}, B::AbstractBkg) where {T,N} = m .- get_bkg(B)




struct BkgMdl{T,N} <: AbstractBkg where {T,N,K}
    b::AbstractArray{T,N}
    R::Regularization
    kwds::K #TODO: defining the keywords used for the optimization of the background
end

#TODO: to uniformize the use of AbstractBkg structures. Each structure needs to
#overload these three methods.
get_bkg(B::Bkg_mdl) = B.b
regul(B::BkgMdl) = B.R(B.b)

function fit_bkg!(B::BkgMdl{T}, 
                  D::CalibratedData) where {T <: AbstractFloat}#TODO: no need for kwds here as they are in the structures
                #   D::CalibratedData,
                #   nonnegative = true,
                #   kwds...) where {T <: AbstractFloat}
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
    
    #TODO: no need to have the if loop here then
    # if nonnegative
    #     vmlmb!(fg_solve!, B.b; lower=T(0), kwds...)
    # else
        # vmlmb!(fg_solve!, B.b; kwds...)
    # end
    vmlmb!(fg_solve!, B.b; B.kwds...)
    
    return B.b
end


################################################### (Sam) I'll do this after ############


"""

#FIXME: overload regul(bkg)

"""
struct ParametrizedBkgMdl{T} <: AbstractBkg where {T}
    θ::AbstractVector{T}
    f::Function
    R::Regularization
    # kwds
end

Base.:-(m::AbstractArray{T,N}, B::ParametrizedBkgMdl{T}) where {T,N} = m - B.f(B.θ)

regul(B::ParametrizedBkgMdl) = B.R(B.f(B.θ))

function fit_bkg!()
# needs to update the calib maps after autocalib
end

