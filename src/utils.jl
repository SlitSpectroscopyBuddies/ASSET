#
# utils.jl
#
# ------------------------------------------------
#
# This file is part of ASSET

"""
    check_bnds(bnds::AbstractVector{Tuple{T,T}}) where {T}
    
Yields an error if lower bounds is upper or equal to the upper bound.
"""
function check_bnds(bnds::AbstractVector{Tuple{T,T}}) where {T}
    for b in bnds
        b[1] >= b[2] && error("bounds must be specified in keywords as (lower_bound, upper_bound)") 
    end
end




"""
    test_tol(temp::T, lasts::AbstractVector{T}, tol::Tuple{T,T}) where {T<:AbstractFloat}

Checks if the absolute difference between `temp` and each element in `lasts` is less than the specified tolerance.

# Arguments
- `temp::T`: The value to compare, where `T` is a subtype of `AbstractFloat`.
- `lasts::AbstractVector{T}`: A vector of previous values to compare against.
- `tol::Tuple{T,T}`: A tuple containing tolerance values. Only `tol[2]` is used as the threshold.

# Returns
- `Bool`: Returns `true` if the maximum of the comparison results is `true`, i.e., if at least one element in `lasts` is within the specified tolerance of `temp`.

---

    test_tol(temp::AbstractArray{T,1}, lasts::AbstractArray{T,2}, tol::Tuple{T,T}) where {T<:AbstractFloat}

Same as the vector version of the method but for arrays.
"""
function test_tol(temp::T, lasts::AbstractVector{T}, tol::Tuple{T,T}) where {T<:AbstractFloat}
    test=[(abs(temp - lasts[k]) < tol[2]) for k=1:length(lasts)]
    return maximum(test)
end


function test_tol(temp::AbstractArray{T,1}, lasts::AbstractArray{T,2}, tol::Tuple{T,T}) where {T<:AbstractFloat}
    test = [(sum(abs.(temp - lasts[:,k]))/sum(abs.(lasts[:,k])) < tol[2]) for k=1:size(lasts)[2]]
    return maximum(test)
end
