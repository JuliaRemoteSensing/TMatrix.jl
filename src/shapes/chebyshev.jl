@doc raw"""
A Chebyshev scatterer defined by

    $r(\theta, \phi)=r_0(1+\varepsilon T_n(\cos\theta))$

in which $T_n(\cos\theta)=\cos n\theta$.

Attributes:

- `rev`: The equivalent volume radius.
- `m`: The complex refractive index.
- `ε`: The deformation parameter, which satisfies $0\le\varepsilon<1$.
- `n`: The degree of the Chebyshev polynomial.
- `λ`: The wavelength of the incident wave.
- `info`: The accompanied information.
"""
struct Chebyshev{T<:Real,CT<:Number,RV,RM,CV,CM} <: AbstractScatterer{T,CT}
    rev::T
    m::CT
    ε::T
    n::Int64
    λ::T
    info::ScattererInfo{RV,RM,CV,CM}
end

has_symmetric_plane(chebyshev::Chebyshev) = chebyshev.n % 2 == 0

function calc_r!(
    scatterer::Chebyshev{T},
    ngauss::Int64,
    x::AbstractArray,
    w::AbstractArray,
    r::AbstractArray,
    dr::AbstractArray,
) where {T<:Real}
    theta_split!(scatterer, ngauss, x, w)
    rev = scatterer.rev
    @assert typeof(scatterer) <: Chebyshev
    e = scatterer.ε
    n = scatterer.n
    dn = T(n)

    a = 1.5e^2 * (4.0dn^2 - 2.0) / (4.0dn^2 - 1.0) + 1.0

    if n % 2 == 0
        a -= 3.0e * (1.0 + 0.25e^2) / (dn^2 - 1.0) + 0.25e^3 / (9.0dn^2 - 1.0)
    end

    r0 = rev / ∛a

    @simd for i in 1:ngauss
        xi = acos(x[i]) * n
        r[i] = r0 * (1.0 + e * cos(xi))
        dr[i] = -r0 * e * n * sin(xi)
    end
end

function theta_split!(scatterer::Chebyshev{T}, ngauss::Int64, x::AbstractArray, w::AbstractArray) where {T<:Real}
    x0, w0 = gausslegendre(T, ngauss)
    x .= x0
    w .= w0
    return
end

function calc_r!(
    scatterer::Chebyshev{Arb},
    ngauss::Int64,
    x::Arblib.ArbVectorLike,
    w::Arblib.ArbVectorLike,
    r::Arblib.ArbVectorLike,
    dr::Arblib.ArbVectorLike,
)
    theta_split!(scatterer, ngauss, x, w)
    rev = scatterer.rev
    e = scatterer.ε
    n = scatterer.n
    a = 1.5e^2 * (4n^2 - 2) / (4n^2 - 1) + 1

    if n % 2 == 0
        a -= 3e * (1 + 0.25e^2) / (n^2 - 1) + 0.25e^3 / (9n^2 - 1)
    end

    r0 = rev / ∛a
    r0₋ = -r0

    for i in 1:ngauss
        xi = Arb(x[i])
        Arblib.acos!(xi, xi)
        Arblib.mul!(xi, xi, n)

        r[i] = e
        Arblib.mul!(r[i], r[i], cos(xi))
        Arblib.add!(r[i], r[i], ARB_ONE)
        Arblib.mul!(r[i], r[i], r0)

        dr[i] = r0₋
        Arblib.mul!(dr[i], dr[i], e)
        Arblib.mul!(dr[i], dr[i], n)
        Arblib.mul!(dr[i], dr[i], sin(xi))
    end
end

function theta_split!(scatterer::Chebyshev{Arb}, ngauss::Int64, x::Arblib.ArbVectorLike, w::Arblib.ArbVectorLike)
    gausslegendre!(Arb, ngauss, x, w)
    return
end
