
"""
    CalibratedData(d, w, ρ_map, λ_map) -> D

Yields a structure `D` containing the data `d` and their respective weights 
`w`, the angular separation `ρ_map` and wavelength `λ_map` maps calibrating the
detector.

"""
struct CalibratedData{T<:AbstractFloat,N,D<:AbstractArray{T,N}}#,
                         #R<:AbstractMatrix{T},
                         #L<:AbstractMatrix{T}} TODO: keep the possibility to have only one map for ρ_map and/or λ_map for all the exposures. 
    d::D
    w::D
    ρ_map::D
    λ_map::D

    function CalibratedData(d::D,
                            w::D,
                            ρ_map::D,
                            λ_map::D) where {T<:Real,N,
                                              D<:AbstractArray{T,N}}#,
                                              #R<:AbstractMatrix{T},
                                              #L<:AbstractMatrix{T}}
        @assert axes(d) == axes(w)== axes(ρ_map) == axes(λ_map)
        #@assert (axes(d,1), axes(d,2)) 
        
        return new{T,N,D}(d, w, ρ_map, λ_map)
    end
end

Base.show(io::IO, D::CalibratedData{T}) where {T} = begin
    print(io,"CalibratedData{$T}:")
    print(io,"\n - scientific data `d` : ",typeof(D.d))
    print(io,"\n - weight of data `w` : ",typeof(D.w))
    print(io,"\n - spatial map `ρ_map` : ",typeof(D.ρ_map))
    print(io,"\n - spectral map `λ_map` : ",typeof(D.λ_map))
end




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

Base.:+(m::AbstractArray{T,N}, B::BkgMdl{T,M}) where {T,N,M} = m .+ B.b
Base.:-(m::AbstractArray{T,N}, B::BkgMdl{T,M}) where {T,N,M} = m .- B.b



regul(B::BkgMdl) = B.R(B.b)

function fit_bkg!(B::BkgMdl{T}, 
                  D::CalibratedData,
                  nonnegative = true,
                  kwds...) where {T <: AbstractFloat}
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

