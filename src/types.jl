
"""
    CalibratedData(d, w, ρ_map, λ_map) -> D

Yields a structure `D` containing the data `d` and their respective weights 
`w`, the angular separation `ρ_map` and wavelength `λ_map` maps calibrating the
detector.

The `Base` methods `axes`, `size`, `eltype` and `show` have been overload to be
used with such a structure.

"""
struct CalibratedData{T<:AbstractFloat,N,D<:AbstractArray{T,N}}
    d::D
    w::D
    ρ_map::D
    λ_map::D

    function CalibratedData(d::D,
                            w::D,
                            ρ_map::D,
                            λ_map::D) where {T<:Real,N,D<:AbstractArray{T,N}}
        @assert axes(d) == axes(w)== axes(ρ_map) == axes(λ_map)
        
        return new{T,N,D}(d, w, ρ_map, λ_map)
    end
end

Base.axes(D::CalibratedData) = axes(D.d)
Base.size(D::CalibratedData) = size(D.d)
Base.eltype(D::CalibratedData) = eltype(typeof(D))
Base.eltype(::Type{<:CalibratedData{T}}) where {T} = T
Base.show(io::IO, D::CalibratedData{T}) where {T} = begin
    print(io,"CalibratedData{$T}:")
    print(io,"\n - scientific data `d` : ",typeof(D.d))
    print(io,"\n - weight of data `w` : ",typeof(D.w))
    print(io,"\n - spatial map `ρ_map` : ",typeof(D.ρ_map))
    print(io,"\n - spectral map `λ_map` : ",typeof(D.λ_map))
end




"""
    AbstractBkg

Type that defines a background structure, when one want to both estimate the
background sources of the data and extracting the spectrum of the object of
interest.

An instance `B` of `AbstractBkg` must contain the parameters necessary to form
the background model as well as a regularization structure. The method `regul`
must be overload to yield the result of the regularization function applied to
the model, or to parameters of the model, of the background.

To estimate the instance `B`, the method `fit_bkg!` must be overload, to be used
as:
```julia
julia> fit_bkg!(B, D, kwds...)
```
with `D` an instance of `CalibratedData`.

Each new structure of type `AbstractBkg` must overload the addition and
subtraction `Base` operators. 

FIXME: See the ideas.jl file

"""
abstract type AbstractBkg end




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

See also [`AbstractBkg`](@ref)

FIXME: See the ideas.jl file
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

