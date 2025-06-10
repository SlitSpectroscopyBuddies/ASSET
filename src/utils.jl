#
# utils.jl
#
# ------------------------------------------------
#
# This file is part of ASSET

"""
    check_bnds(bnds)
    
Yields an error if lower bounds is upper or equal to the upper bound.
"""
function check_bnds(bnds::AbstractVector{Tuple{T,T}}) where {T}
    for b in bnds
        b[1] >= b[2] && error("bounds must be specified in keywords as (lower_bound, upper_bound)") 
    end
end


function test_tol(temp::T, lasts::AbstractVector{T}, tol::Tuple{T,T}) where {T<:AbstractFloat}
    test=[(abs(temp - lasts[k]) < tol[2]) for k=1:length(lasts)]
    return maximum(test)
end


function test_tol(temp::AbstractArray{T,1}, lasts::AbstractArray{T,2}, tol::Tuple{T,T}) where {T<:AbstractFloat}
    test = [(sum(abs.(temp - lasts[:,k]))/sum(abs.(lasts[:,k])) < tol[2]) for k=1:size(lasts)[2]]
    return maximum(test)
end
