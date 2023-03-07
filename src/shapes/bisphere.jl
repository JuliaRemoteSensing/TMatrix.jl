@doc raw"""
A bisphere scatterer, which consists of two identical spheres touching each other. The coordinate origin is set at their tangency point.

Attributes:

- `m`: The complex refractive index.
- `r`: The radius of the sphere.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Bisphere{T <: Real, CT <: Number, RV, RM, CV, CM} <: AbstractScatterer{T, CT}
    m::CT
    r::T
    λ::T
    info::ScattererInfo{RV, RM, CV, CM}
end

@doc raw"""
Construct a bisphere.
"""
function Bisphere(T::Type{<:Real} = Float64; r::Real, m::Number, λ::Real)
    m = Complex{T}(m)
    r = T(r)
    λ = T(λ)
    return Bisphere(m, r, λ, ScattererInfo(T))
end

has_symmetric_plane(bicone::Bisphere) = true

volume_equivalent_radius(bicone::Bisphere) = ∛2 * bicone.r

function calc_r!(scatterer::Bisphere{T},
                 ngauss::Int64,
                 x::AbstractArray,
                 w::AbstractArray,
                 r::AbstractArray,
                 dr::AbstractArray) where {T <: Real}
    theta_split!(scatterer, ngauss, x, w)

    @simd for i in 1:(ngauss ÷ 2)
        cosθ = abs(x[i])
        sinθ = √(1 - cosθ^2)
        r[i] = 2scatterer.r * cosθ
        r[ngauss + 1 - i] = r[i]
        dr[i] = -2scatterer.r * sinθ
        dr[ngauss + 1 - i] = -dr[i]
    end
end

function theta_split!(scatterer::Bisphere{T}, ngauss::Int64, x::AbstractArray,
                      w::AbstractArray) where {T <: Real}
    ng = ngauss ÷ 2
    x1, w1 = gausslegendre(T, ng)
    @. x[1:ng] = 0.5(x1 - 1)
    @. w[1:ng] = 0.5w1
    @. x[(ng + 1):ngauss] = (-1) * x[ng:-1:1]
    @. w[(ng + 1):ngauss] = w[ng:-1:1]
end
