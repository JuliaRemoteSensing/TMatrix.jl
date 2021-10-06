@doc raw"""
A spheroid scatterer.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `a_to_c`: The ratio $a/c$ of the horizontal to rotational axes.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Spheroid{T<:Real,CT<:Number,RV,RM,CV,CM} <: AbstractScatterer{T,CT}
    rev::T
    m::CT
    a_to_c::T
    λ::T
    info::ScattererInfo{RV,RM,CV,CM}
end

has_symmetric_plane(spheroid::Spheroid) = true

function calc_r!(
    scatterer::Spheroid{T},
    ngauss::Int64,
    x::AbstractArray,
    w::AbstractArray,
    r::AbstractArray,
    dr::AbstractArray,
) where {T<:Real}
    theta_split!(scatterer, ngauss, x, w)
    rev = scatterer.rev
    e = scatterer.a_to_c
    a = rev * ∛e

    @simd for i in 1:(ngauss ÷ 2)
        cosθ = x[i]
        sinθ = √(1.0 - cosθ^2)
        r[i] = a * √(1.0 / (e^2 * cosθ^2 + sinθ^2))
        r[ngauss + 1 - i] = r[i]
        dr[i] = r[i]^3 * cosθ * sinθ * (e^2 - 1.0) / a^2
        dr[ngauss + 1 - i] = -dr[i]
    end
end

function calc_r!(
    scatterer::Spheroid{Arb},
    ngauss::Int64,
    x::Arblib.ArbVectorLike,
    w::Arblib.ArbVectorLike,
    r::Arblib.ArbVectorLike,
    dr::Arblib.ArbVectorLike,
)
    theta_split!(scatterer, ngauss, x, w)
    rev = scatterer.rev

    e = scatterer.a_to_c
    e² = e^2
    a = rev * ∛e
    a² = a^2
    ratio = e² - 1

    for i in 1:(ngauss ÷ 2)
        cosθ = x[i]
        cos²θ = cosθ^2
        sin²θ = 1 - cos²θ
        sinθ = √sin²θ

        r[i] = e²
        Arblib.mul!(r[i], r[i], cos²θ)
        Arblib.add!(r[i], r[i], sin²θ)
        Arblib.div!(r[i], ARB_ONE, r[i])
        Arblib.sqrt!(r[i], r[i])
        Arblib.mul!(r[i], r[i], a)
        r[ngauss + 1 - i] = r[i]

        Arblib.pow!(dr[i], r[i], UInt64(3))
        Arblib.mul!(dr[i], dr[i], cosθ)
        Arblib.mul!(dr[i], dr[i], sinθ)
        Arblib.mul!(dr[i], dr[i], ratio)
        Arblib.div!(dr[i], dr[i], a²)
        Arblib.mul!(dr[ngauss + 1 - i], dr[i], ARB_ONE_NEG)
    end
end
