#
# types.jl
#
# ------------------------------------------------
#
# This file is part of ASSET


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

Type that defines a background structure, that is used to specify a certain
model of the background present in the data and its associated regularization.
This structure is used when one want to both estimate the background sources of
the data and extracting the spectrum of the object of interest.

To define a structure `B` of type `AbstractBkg`, the user needs to overload three
methods:

```julia
julia> get_bkg(B)
```
which yields the background as an `AbstractArray` from parameters contained in
`B`.

```julia
julia> regul(B)
```
which yields the result of the regularization term(s), applied to the curent
background model, that are used in the a posteriori likelihood minimization.

```julia
julia> fit_bkg!(B, D)
```
where `D` is a `CalibratedData`. This last function is the one called by
`extract_spectrum` to estimate the background when assuming the object spectrum
known (alternate estimation). 

See also [`BkgMdl`](@ref)

"""
abstract type AbstractBkg end

Base.:+(m::AbstractArray{T,N}, B::AbstractBkg) where {T,N} = m .+ get_bkg(B)
Base.:-(m::AbstractArray{T,N}, B::AbstractBkg) where {T,N} = m .- get_bkg(B)



"""
"""
abstract type ParametricPSF{N} <: AbstractPSF{N} end




"""
"""
abstract type NonParametricPSF{N} <: AbstractPSF{N} end

