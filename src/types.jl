

"""


"""
abstract type BkgModel end



"""
    Cost

Abstract type which behaves like a `Mapping` type of `LazyAlgebra` and which 
refers to `call` and `call!` functions when applied to an element as an 
operator.

"""
abstract type Cost <: Mapping end

(C::Cost)(a, x) = call(a, C, x)
(C::Cost)(x) = call(1.0, C, x)
(C::Cost)(a, x::D, g::D; kwds...) where {D} = call!(a, C, x, g; kwds...)
(C::Cost)(x::D, g::D; kwds...) where {D} = call!(1.0, C, x, g; kwds...)



"""
    call!([a::Real=1,] f, x, g; incr::Bool = false)

gives back the value of a*f(x) while updating (in place) its gradient in g. 
incr is a boolean indicating if the gradient needs to be incremanted or 
reseted.

"""
call!(f, x, g; kwds...) = call!(1.0, f, x, g; kwds...)

"""
    call([a::Real=1,] f, x)

yields the value of a*f(x).

"""
call(f, x) = call(1.0, f, x)




"""
    Regularization

Abstract type which obeys the `call!` and `call` laws, and behaves as a `Cost` type. 

"""
abstract type Regularization <: Cost end




"""
    Regul([mu::Real=1,] f)

yields a structure of supertype `Regularization` containing the element necessary 
to compute a regularization term `f` with it's multiplier `mu`. Creating an instance 
of Regul permits to use the `call` and `call!` functions specialized.

Getters of an instance `R` can be imported with `InversePbm.` before:
 - multiplier(R) # gets the multiplier `mu`.
 - func(R)       # gets the computing process `f`.

# Example
```julia
julia> R = Regul(a, MyRegul)
Regul:
 - level `mu` : a
 - function `func` : MyRegul
```
with `MyRegul` a strucure which defines the behaviors of `call([a,] MyRegul, x)` and 
`call!([a,] MyRegul, x, g [;incr])` on an `AbstractArray` `x` and its gradient `g`.

To apply it to `x`, use:
```julia
julia> R([a,] x)            # apply call([a,] R, x)
julia> R([a,] x, g [;incr]) # apply call!([a,] R, x, g [;incr])
```

Applying a scalar `b` to an insance `R` yields a new instance of `Regul`:
```julia
julia> R = b*Regul(a, MyRegul)
Regul:
 - level `mu` : b*a
 - function `func` : MyRegul
```

"""
struct Regul{T<:Real,F} <: Regularization
    mu::T # multiplier
    f::F # function
end
multiplier(R::Regul) = R.mu
func(R::Regul) = R.f

Regul(f) = Regul(1.0, f)

*(a::Real, R::Regul) = Regul(a*multiplier(R), func(R))
Base.show(io::IO, R::Regul) = begin
    print(io,"Regul:")
    print(io,"\n - level `mu` : ",multiplier(R))
    print(io,"\n - function `func` : ",func(R))
end

function call!(a::Real,
    R::Regul,
    x::AbstractArray{T,N},
    g::AbstractArray{T,N};
    incr::Bool = false) where {T,N}

    return call!(a*multiplier(R), func(R), x, g; incr=incr)
end

function call!(R::Regul,
    x::AbstractArray{T,N},
    g::AbstractArray{T,N};
    incr::Bool = false) where {T,N}

    return call!(multiplier(R), func(R), x, g; incr=incr)
end

function call(a::Real,
    R::Regul,
    x::AbstractArray{T,N}) where {T,N}
    
    return call(a*multiplier(R), func(R), x)
end

function call(R::Regul,
    x::AbstractArray{T,N}) where {T,N}
    
    return call(multiplier(R), func(R), x)
end




"""
    HomogeneousRegularization(mu, deg, f)

yields a structure representing an homogeneous regularization and containing
the multiplier which tunes it, its degree and the computing process/function
which gives its value.

"""
struct HomogeneousRegularization{F}
    mu::Float64
    deg::Float64
    f::F # function
end
multiplier(R::HomogeneousRegularization) = R.mu
degree(R::HomogeneousRegularization) = R.deg
func(R::HomogeneousRegularization) = R.f

HomogeneousRegularization(mu::Real, f) = HomogeneousRegularization(mu, degree(f), f)
HomogeneousRegularization(mu::Real, R::HomogeneousRegularization) = HomogeneousRegularization(mu, degree(R), func(R))
HomogeneousRegularization(f) = HomogeneousRegularization(1.0, f)


Base.:(*)(a::Real, R::HomogeneousRegularization) =
    HomogeneousRegularization(a*multiplier(R), degree(R), func(R))
Base.show(io::IO, R::HomogeneousRegularization) = begin
    print(io,"HomogeneousRegularization:")
    print(io,"\n - level `mu` : ",multiplier(R))
    print(io,"\n - degree `deg` : ",degree(R))
    print(io,"\n - function `func` : ",func(R))
end


function call!(α::Real,
               R::HomogeneousRegularization,
               x::AbstractArray{T,N},
               g::AbstractArray{T,N};
               incr::Bool = false) where {T,N}
    return call!(α*multiplier(R), func(R), x, g; incr=incr)
end

function call!(R::HomogeneousRegularization,
               x::AbstractArray{T,N},
               g::AbstractArray{T,N};
               incr::Bool = false) where {T,N}
    return call!(multiplier(R), func(R), x, g; incr=incr)
end

function call(α::Real,
              R::HomogeneousRegularization,
              x::AbstractArray{T,N}) where {T,N}
    return call(α*multiplier(R), func(R), x)
end

function call(R::HomogeneousRegularization,
              x::AbstractArray{T,N}) where {T,N}
    return call(multiplier(R), func(R), x)
end

