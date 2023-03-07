@doc raw"""
A Chebyshev scatterer defined by

    $r(\theta, \phi)=r_0(1+\varepsilon T_n(\cos\theta))$

in which $T_n(\cos\theta)=\cos n\theta$.

Attributes:

- `m`: The complex refractive index.
- `r₀`: The radius of the base sphere.
- `ε`: The deformation parameter, which satisfies $-1\le\varepsilon<1$.
- `n`: The degree of the Chebyshev polynomial.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Chebyshev{T <: Real, CT <: Number, RV, RM, CV, CM} <: AbstractScatterer{T, CT}
    m::CT
    r₀::T
    ε::T
    n::Int64
    λ::T
    info::ScattererInfo{RV, RM, CV, CM}
end

function Chebyshev(T::Type{<:Real} = Float64;
                   r₀::Real = 0,
                   ε::Real,
                   n::Int64,
                   rev::Real = 0,
                   rea::Real = 0,
                   rmax::Real = 0,
                   m::Number,
                   λ::Real)
    @assert -1 < ε < 1

    m = Complex{T}(m)
    r₀ = T(r₀)
    ε = T(ε)
    rev = T(rev)
    rea = T(rea)
    rmax = T(rmax)
    λ = T(λ)
    info = ScattererInfo(T)

    if !iszero(r₀)
        nothing
    else
        a = R_3_2 * ε^2 * (4n^2 - 2) / (4n^2 - 1) + 1
        if iseven(n)
            a -= 3ε * (1 + R_1_4 * ε^2) / (n^2 - 1) + R_1_4 * ε^3 / (9n^2 - 1)
        end

        if !iszero(rev)
            r₀ = rev / ∛a
        elseif !iszero(rea)
            en = ε * n
            x, w = gausslegendre(T, DEFAULT_NGCHEB[])
            s = zero(T)
            v = zero(T)
            @simd for i in 1:DEFAULT_NGCHEB[]
                θ = acos(x[i])
                nθ = n * θ
                sinθ = sin(θ)
                sinnθ = sin(nθ)
                cosnθ = cos(nθ)
                b = 1 + ε * cosnθ
                ens = en * sinnθ
                s += w[i] * b * √(b^2 + ens^2)
                v += w[i] * (sinθ * b + x[i] * ens) * sinθ * b^2
            end
            rs = √(R_1_2 * s)
            rv = ∛(R_3_4 * v)

            rev = rea / rs * rv
            r₀ = rev / ∛a
        elseif !iszero(rmax)
            r₀ = rmax / (1 + abs(ε))
        else
            throw(ArgumentError("Cannot construct a valid Chebyshev particle with the given parameters."))
        end
    end

    return Chebyshev(m, r₀, ε, n, λ, info)
end

has_symmetric_plane(chebyshev::Chebyshev) = iseven(chebyshev.n)

function volume_equivalent_radius(chebyshev::Chebyshev)
    ε = chebyshev.ε
    n = chebyshev.n

    a = R_3_2 * ε^2 * (4n^2 - 2) / (4n^2 - 1) + 1
    if iseven(n)
        a -= 3ε * (1 + R_1_4 * ε^2) / (n^2 - 1) + R_1_4 * ε^3 / (9n^2 - 1)
    end

    return ∛a * chebyshev.r₀
end

function calc_r!(scatterer::Chebyshev{T},
                 ngauss::Int64,
                 x::AbstractArray,
                 w::AbstractArray,
                 r::AbstractArray,
                 dr::AbstractArray) where {T <: Real}
    theta_split!(scatterer, ngauss, x, w)
    ε = scatterer.ε
    n = scatterer.n
    r₀ = scatterer.r₀

    @simd for i in 1:ngauss
        xi = acos(x[i]) * n
        r[i] = r₀ * (1.0 + ε * cos(xi))
        dr[i] = -r₀ * ε * n * sin(xi)
    end
end

function theta_split!(scatterer::Chebyshev{T}, ngauss::Int64, x::AbstractArray,
                      w::AbstractArray) where {T <: Real}
    x0, w0 = gausslegendre(T, ngauss)
    x .= x0
    w .= w0
    return
end

function calc_r!(scatterer::Chebyshev{Arb},
                 ngauss::Int64,
                 x::Arblib.ArbVectorLike,
                 w::Arblib.ArbVectorLike,
                 r::Arblib.ArbVectorLike,
                 dr::Arblib.ArbVectorLike)
    theta_split!(scatterer, ngauss, x, w)
    ε = scatterer.ε
    n = scatterer.n
    r₀ = scatterer.r₀
    r₀₋ = -r₀

    for i in 1:ngauss
        xi = Arb(x[i])
        Arblib.acos!(xi, xi)
        Arblib.mul!(xi, xi, n)

        r[i] = ε
        Arblib.mul!(r[i], r[i], cos(xi))
        Arblib.add!(r[i], r[i], ARB_ONE)
        Arblib.mul!(r[i], r[i], r₀)

        dr[i] = r₀₋
        Arblib.mul!(dr[i], dr[i], ε)
        Arblib.mul!(dr[i], dr[i], n)
        Arblib.mul!(dr[i], dr[i], sin(xi))
    end
end

function theta_split!(scatterer::Chebyshev{Arb}, ngauss::Int64, x::Arblib.ArbVectorLike,
                      w::Arblib.ArbVectorLike)
    gausslegendre!(Arb, ngauss, x, w)
    return
end
