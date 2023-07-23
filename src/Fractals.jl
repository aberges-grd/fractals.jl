module Fractals

include("algorithms.jl")
using .Algorithms

export
    JuliaSet,
    MandelbrotSet,
    Canvas,
    DistanceEstimation,
    EscapeTime,
    construct_mesh,
    compute

abstract type IterativeFractal end

struct JuliaSet <: IterativeFractal
    c::Complex
    f::Function
    df::Function # Only needed for the distance estimation algorithm
    max_iter::Int
    radius::Float64
end


struct MandelbrotSet <: IterativeFractal
    f::Function
    df::Function # Only needed for the distance estimation algorithm
    max_iter::Int
    radius::Float64
end

struct Canvas
    # Resolution
    width::Int
    height::Int
    # Complex limits
    topleft::Complex
    bottomright::Complex
end

struct DistanceEstimation end
struct EscapeTime end

compute(s::IterativeFractal, canvas::Canvas) = compute(s, canvas, EscapeTime())

function compute(s::JuliaSet, canvas::Canvas, ::EscapeTime)::Array{Float64}
    construct_mesh(canvas) .|> x -> escape_time(z -> s.f(z, s.c), x, s.max_iter, s.radius)
end

function compute(s::MandelbrotSet, canvas::Canvas, ::EscapeTime)::Array{Float64}
    construct_mesh(canvas) .|> c -> escape_time(z -> s.f(z, c), 0.0 + 0.0im, s.max_iter, s.radius)
end

function compute(s::JuliaSet, canvas::Canvas, ::DistanceEstimation)::Array{Float64}
    construct_mesh(canvas) .|> x -> distance_estimation(z -> s.f(z, s.c), s.df, x, s.max_iter, s.radius)
end

function compute(s::MandelbrotSet, canvas::Canvas, ::DistanceEstimation)::Array{Float64}
    construct_mesh(canvas) .|> c -> distance_estimation(z -> s.f(z, c), s.df, 0.0 + 0.0im, s.max_iter, s.radius)
end

function construct_mesh(canvas::Canvas)::Array{ComplexF64}
    re₀, re₁ = real(canvas.topleft), real(canvas.bottomright)
    im₀, im₁ = imag(canvas.topleft), imag(canvas.bottomright)
    return [complex(a, b)
            for a in range(re₀, length=canvas.width, stop=re₁),
            b in range(im₀, length=canvas.height, stop=im₁)]
end

end # module fractals
