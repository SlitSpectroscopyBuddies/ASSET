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

