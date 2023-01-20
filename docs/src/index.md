# TMatrix

[TMatrix.jl](https://github.com/JuliaRemoteSensing/TMatrix.jl) is a Julia package for the calculation of the T-matrix of a dielectric particle, using the Extended Boundary Condition Method (EBCM), also known as the Null-Field Method (NFM). The T-matrix is a compact representation of the scattering properties of a particle, and is widely used.

## Installation

```julia-repl
julia> using Pkg
julia> Pkg.add("TMatrix")
```

## Usage

First, one needs to define a scatterer. Currently, the following shapes are supported:

- Spheroids
- Cylinders
- Chebyshev particles
- Capsules, i.e., cylinders with half-spherical caps on both ends.
- Bi-cones

One can refer to the [API](@ref) for the detailed definitions of these shapes.

Take a spheroid as an example:

```julia
scatterer = TMatrix.Spheroid(a = 1.0, c = 2.0, m = 1.311, λ = 2π)
```

Then one can calculate the T-matrix of the scatterer:

```julia
T = TMatrix.tmatrix(scatterer)
```

With the T-matrix calculated, one can calculate the scattering properties of the scatterer.

For the fixed orientation case, the amplitude scattering matrix `S` and the phase matrix `Z` can be calculated as:

```julia
α = 30.0
β = 60.0
ϑᵢ = 30.0
ϑₛ = 60.0
φᵢ = 15.0
φₛ = 90.0
S, Z = TMatrix.calc_SZ(scatterer, α, β, ϑᵢ, ϑₛ, φᵢ, φₛ, T)
```

For the random orientation case, one can first calculate the expansion coefficients of the T-matrix:

```julia
coeffs = TMatrix.calc_expansion_coefficients(scatterer, T)
```

Here the expansion coefficients are given in the form `(α₁, α₂, α₃, α₄, β₁, β₂)`, where each element is a `Vector`.

Then the orientation-averaged scattering matrix can be evaluated at an arbitrary scattering angle:

```julia
Θ = 10.0
F = TMatrix.calc_F(coeffs..., Θ)
```

Here `F` is in the form of `(F₁₁, F₁₂, F₂₂, F₃₃, F₃₄, F₄₄)`. The scattering matrix can be constructed from `F` as:

```math
\begin{bmatrix}
F_{11} & F_{12} & 0 & 0 \\
F_{12} & F_{22} & 0 & 0 \\
0 & 0 & F_{33} & F_{34} \\
0 & 0 & -F_{34} & F_{44}
\end{bmatrix}
```
