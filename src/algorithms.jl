module Algorithms

export escape_time, distance_estimation

"""
    Computes the escape time.

    Returns a number between 0 and 1 (1 being the maximum time).
"""
function escape_time(f::Function, z::Complex, max_iter::Int, radius::Float64)::Float64
    for n = 1:max_iter
        z = f(z)
        if abs(z) > radius
            return n / max_iter
        end
    end
    return 1.0
end

"""Estimates distance from point \$z\$ to the set."""
function distance_estimation(f::Function, df::Function, z::Complex, max_iter::Int, radius::Float64)::Float64
    dz = 1.0 + 0.0im
    for _ = 1:max_iter
        z, dz = f(z), df(z, dz)
        if abs(z) > radius
            break
        end
    end

    return abs(z) * log(abs(z)) / abs(dz)
end

end