@doc raw"""
A spheroid scatterer.

Attributes:

- `m`: The complex refractive index.
- `a`: Length of the semi-major axis.
- `c`: Length of the semi-minor axis.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Spheroid{T<:Real,CT<:Number,RV,RM,CV,CM} <: AbstractScatterer{T,CT}
    m::CT
    a::T
    c::T
    λ::T
    info::ScattererInfo{RV,RM,CV,CM}
end

@doc raw"""
Construct a spheroid. `λ` and `m` must be provided, and the size parameters are considered according to the following priority:

1. `a` and `c`
2. `a` and `a_to_c`
3. `c` and `a_to_c`
4. `rev` and `a_to_c`
5. `rea` and `a_to_c`
6. `rmax` and `a_to_c`

If none of the above is hit, an `ArgumentError` will be thrown.

"""
function Spheroid(
    T::Type{<:Real} = Float64;
    a::Real = 0,
    c::Real = 0,
    rev::Real = 0,
    rea::Real = 0,
    rmax::Real = 0,
    a_to_c::Real = 0,
    m::Number,
    λ::Real,
)
    m = Complex{T}(m)
    a = T(a)
    c = T(c)
    rev = T(rev)
    rea = T(rea)
    a_to_c = T(a_to_c)
    λ = T(λ)
    info = ScattererInfo(T)

    if !iszero(a) && !iszero(c)
        nothing
    elseif !iszero(a) && !iszero(a_to_c)
        c = a / a_to_c
    elseif !iszero(c) && !iszero(a_to_c)
        a = c * a_to_c
    elseif !iszero(rev) && !iszero(a_to_c)
        a = rev * ∛a_to_c
        c = a / a_to_c
    elseif !iszero(rea) && !iszero(a_to_c)
        d = a_to_c

        if abs(d - 1.0) < eps(d)
            ratio = 1.0
        elseif d > 1.0
            e = √(1 - 1 / d^2)
            ratio = 1 / √(R_1_4 * (2d^R_2_3 + d^(-R_4_3) * log((1 + e) / (1 - e)) / e))
        elseif d < 1.0
            e = √(1 - d^2)
            ratio = 1 / √(R_1_2 * (d^R_2_3 + d^(-R_1_3) * asin(e) / e))
        end

        rev = rea * ratio
        a = rev * ∛a_to_c
        c = a / a_to_c
    elseif !iszero(rmax) && !iszero(a_to_c)
        rev = a_to_c > 1.0 ? rmax / ∛a_to_c : rmax * a_to_c^R_2_3
        a = rev * ∛a_to_c
        c = a / a_to_c
    else
        throw(ArgumentError("Cannot construct a valid spheroid with the given parameters."))
    end

    return Spheroid(m, a, c, λ, info)
end

has_symmetric_plane(spheroid::Spheroid) = true

volume_equivalent_radius(spheroid::Spheroid) = ∛(spheroid.a^2 * spheroid.c)

function calc_r!(
    spheroid::Spheroid{T},
    ngauss::Int64,
    x::AbstractArray,
    w::AbstractArray,
    r::AbstractArray,
    dr::AbstractArray,
) where {T<:Real}
    theta_split!(spheroid, ngauss, x, w)
    e = spheroid.a / spheroid.c
    a = spheroid.a

    @simd for i in 1:(ngauss ÷ 2)
        cosθ = x[i]
        sinθ = √(1.0 - cosθ^2)
        r[i] = a * √(1.0 / (e^2 * cosθ^2 + sinθ^2))
        r[ngauss + 1 - i] = r[i]
        dr[i] = r[i]^3 * cosθ * sinθ * (e^2 - 1.0) / a^2
        dr[ngauss + 1 - i] = -dr[i]
    end
end

function theta_split!(spheroid::Spheroid{T}, ngauss::Int64, x::AbstractArray, w::AbstractArray) where {T<:Real}
    x0, w0 = gausslegendre(T, ngauss)
    x .= x0
    w .= w0
    return
end

function calc_r!(
    spheroid::Spheroid{Arb},
    ngauss::Int64,
    x::Arblib.ArbVectorLike,
    w::Arblib.ArbVectorLike,
    r::Arblib.ArbVectorLike,
    dr::Arblib.ArbVectorLike,
)
    theta_split!(spheroid, ngauss, x, w)

    e = spheroid.a / spheroid.c
    e² = e^2
    a = spheroid.a
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

function theta_split!(spheroid::Spheroid{Arb}, ngauss::Int64, x::Arblib.ArbVectorLike, w::Arblib.ArbVectorLike)
    gausslegendre!(Arb, ngauss, x, w)
    return
end
