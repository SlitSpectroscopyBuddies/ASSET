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


"""
"""
struct ChromaticSeriesExpansionsInterpolator{T<:AbstractFloat,
                                             N,
                                  Trows<:NTuple{N,Int},#Union{NTuple{2,Int},NTuple{3,Int}},
                                  Tker<:Union{Kernel{T,2},Kernel{T,4}},
                                  Tmap<:AbstractArray{T},#Union{AbstractArray{T,2},AbstractArray{T,3}},
                                  Tx<:AbstractRange{T}} <: LinearMapping
    cols::NTuple{2,Int} 
    rows::Trows 
    ker::Tker
    X::Tmap
    Λ::Tmap
    x::Tx
    λref::Real
    a::Real
end

function ChromaticSeriesExpansionsInterpolator{T}(ker::Kernel,
                                X::AbstractArray{T,N},#Union{AbstractArray{<:Any,2},AbstractArray{<:Any,3}},
                                Λ::AbstractArray{T,N},#Union{AbstractArray{<:Any,2},AbstractArray{<:Any,3}},
                                x::AbstractRange;
                                order=1,
                                λref = maximum(Λ),
                                a=0.) where {T <: AbstractFloat,N}
    @assert size(X) == size(Λ)
    ChromaticSeriesExpansionsInterpolator((length(x), order), 
                               size(X), 
                               convert(Kernel{T},ker), 
                               T.(X),
                               T.(Λ),
                               x,
                               λref,
                               a)
end
function ChromaticSeriesExpansionsInterpolator(ker::Kernel,
                                X::AbstractArray{<:Any,N},
                                Λ::AbstractArray{<:Any,N},
                                x::AbstractRange;
                                order=1,
                                λref = maximum(Λ),
                                a=0.) where{N}
    T=Float64
    ChromaticSeriesExpansionsInterpolator{T}(ker, X, Λ, x; order=order, λref=λref, a=a)
end


function vcreate(::Type{LazyAlgebra.Direct}, A::ChromaticSeriesExpansionsInterpolator{T},
                 x::AbstractArray{T,2}, scratch::Bool = false) where {T <: AbstractFloat}
    @assert !Base.has_offset_axes(x)
    @assert size(x) == A.cols
    N=length(A.rows)
    Array{T,N}(undef, A.rows)
end

function vcreate(::Type{LazyAlgebra.Adjoint}, A::ChromaticSeriesExpansionsInterpolator{T},
                 x::AbstractArray{T,N}, scratch::Bool = false) where {T <: AbstractFloat,N}
    @assert !Base.has_offset_axes(x)
    @assert size(x) == A.rows
    Array{T,2}(undef, A.cols)
end


function apply!(α::Real,
                ::Type{LazyAlgebra.Direct},
                R::ChromaticSeriesExpansionsInterpolator{T},
                src::AbstractArray{T,2},
                scratch::Bool,
                β::Real,
                dst::AbstractArray{T,N}) where {T<:AbstractFloat,N}
    @assert size(src) == R.cols
    @assert size(dst) == R.rows
    #Direct!(dst,src,R.ker,R.X,R.Y, R.Λ, R.x, R.y, R.λ);
    
    fill!(dst,zero(T))
    xs=step(R.x)
    
    cx = (first(R.x)+ last(R.x))/2
    qx = 1/xs
    cqrx =  (1+ length(R.x))/2- qx*cx 
    xlim = LinearInterpolators.limits(R.ker,length(R.x))
      
    S=4 #FIXME get directly the size of the kernel support
    Λc = zeros(Int64,S)
    Λw = zeros(S)
    Xc = zeros(Int64,S)
    Xw = zeros(S)
    
    λmax=R.λref   
    
    @inbounds for n=1:length(dst)
        iszero(R.Λ[n]) && continue
            nsum = zero(T)
            Λn = R.Λ[n]
            Xn = R.X[n]
            γ = sqrt((R.a^2 + λmax/Λn)/(R.a^2+1))
            Xn *= γ
            Xpos = qx*Xn  + cqrx 
                                 
            Xc[1],Xc[2],Xc[3],Xc[4],
            Xw[1],Xw[2],Xw[3],Xw[4] = LinearInterpolators.getcoefs(R.ker,xlim,Xpos)  
            for o=R.cols[2]:-1:1
            
                @simd for i=1:S
                    I=Xc[i]
                    nsum += Xw[i] * src[I,o]
                end
                
                nsum *= γ  
            end

            dst[n] = nsum 
    end
    
    dst .*=α
    dst .+= β
    
    return dst

end

function apply!(α::Real,
                ::Type{LazyAlgebra.Adjoint},
                R::ChromaticSeriesExpansionsInterpolator{T},
                src::AbstractArray{T,N},
                scratch::Bool,
                β::Real,
                dst::AbstractArray{T,2}) where {T<:AbstractFloat,N}
    @assert size(src) == R.rows
    @assert size(dst) == R.cols
    
    fill!(dst,zero(T))
    xs=step(R.x)
    

    cx = (first(R.x)+ last(R.x))/2
    qx = 1/xs
    cqrx =  (1+ length(R.x))/2- qx*cx 
    xlim = LinearInterpolators.limits(R.ker,length(R.x))
    
    S=4 #FIXME get directly the size of the kernel support
    Xc = zeros(Int64,S)
    Xw = zeros(S)
    
    λmax= R.λref
    
    @inbounds for n=1:length(src)
        iszero(R.Λ[n]) && continue
            sn = src[n]
            
            Λn = R.Λ[n]
            Xn = R.X[n]
            γ = sqrt((R.a^2 + λmax/Λn)/(R.a^2+1))
            
            Xn *= γ
            
            Xpos = qx*Xn + cqrx 
            Xc[1],Xc[2],Xc[3],Xc[4],
            Xw[1],Xw[2],Xw[3],Xw[4] = LinearInterpolators.getcoefs(R.ker,xlim,Xpos)
       
            for o=1:R.cols[2]
                sn *= γ
   
                @simd for i=1:S
                    I=Xc[i]
                    dst[I,o] += Xw[i] * sn 
                end
            end
    end
    
    dst .*=α
    dst .+= β
    
    return dst
end

