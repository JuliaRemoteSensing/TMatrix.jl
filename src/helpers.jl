function calc_an(T::Type{<:Real}, nmax::Int64)
    an = zeros(T, nmax)
    for n in 1:nmax
        an[n] = n * (n + 1)
    end
    return an
end

function calc_an(T::Type{<:Arblib.ArbLike}, nmax::Int64)
    an = ArbRefVector(nmax)
    for n in 1:nmax
        an[n] = n * (n + 1)
    end
    return an
end

function calc_ann(T::Type{<:Real}, nmax::Int64)
    ann = zeros(T, nmax, nmax)
    for n1 in 1:nmax
        for n2 in 1:nmax
            ann[n1, n2] = √T((2n1 + 1) * (2n2 + 1) // (n1 * (n1 + 1) * n2 * (n2 + 1))) / 2
        end
    end
    return ann
end

function calc_ann(T::Type{<:Arblib.ArbLike}, nmax::Int64)
    ann = ArbRefMatrix(nmax, nmax)
    half = Arb(1 // 2)
    for n1 in 1:nmax
        for n2 in 1:nmax
            ann[n1, n2] = (2n1 + 1) * (2n2 + 1) // (n1 * (n1 + 1) * n2 * (n2 + 1))
            Arblib.sqrt!(ann[n1, n2], ann[n1, n2])
            Arblib.mul!(ann[n1, n2], ann[n1, n2], half)
        end
    end
    return ann
end

function calc_sig(T::Type{<:Real}, nmax::Int64)
    sig = zeros(T, nmax)
    for n in 1:nmax
        sig[n] = 1 - 2 * (n % 2)
    end
    return sig
end

function calc_sig(T::Type{<:Arblib.ArbLike}, nmax::Int64)
    sig = ArbRefVector(nmax)
    for n in 1:nmax
        sig[n] = 1 - 2 * (n % 2)
    end
    return sig
end

function vig!(nmax::Int64, m::Int64, x::T, dv1::AbstractArray{T},
              dv2::AbstractArray{T}) where {T <: Real}
    if x < -1 || x > 1 || abs(1.0 - abs(x)) < eps(x)
        error("Constraint violated: x ∈ (-1, 1)")
    end

    if m < 0
        error("Constraint violated: m >= 0")
    end

    if nmax <= 0
        error("Constraint violated: nmax >= 1")
    end

    NRT = nonref(x)

    a = one(x)
    qs = √(one(x) - x * x)
    qs1 = one(x) / qs
    if m == 0
        d1 = one(x)
        d2 = x
        for i in 1:nmax
            d3 = (NRT(2i + 1) * x * d2 - i * d1) / (i + 1)
            der = qs1 * i * (i + 1) / (2i + 1) * (d3 - d1)
            dv1[i] = d2
            dv2[i] = der
            d1, d2 = d2, d3
        end
    else
        for i in 1:m
            a *= √NRT((2i - 1) // 2i) * qs
        end
        d1 = zero(x)
        d2 = a
        for i in m:nmax
            qnm = √NRT(i^2 - m^2)
            qnm1 = √NRT((i + 1)^2 - m^2)
            d3 = ((2i + 1) * x * d2 - qnm * d1) / qnm1
            der = qs1 * (-(i + 1) * qnm * d1 + i * qnm1 * d3) / (2i + 1)
            dv1[i] = d2
            dv2[i] = der
            d1, d2 = d2, d3
        end
    end
end

function vig(nmax::Int64, m::Int64, x::T) where {T <: Real}
    dv1 = zeros(T, nmax)
    dv2 = zeros(T, nmax)
    vig!(nmax, m, x, dv1, dv2)
    return dv1, dv2
end

function vigampl(nmax::Int64, m::Int64, x::T) where {T <: Real}
    if abs(1.0 - abs(x)) < eps(x)
        if m != 1
            return zeros(T, nmax), zeros(T, nmax)
        else
            if x < 0.0
                dv1 = [(-1)^(i + 1) * T(0.5) * √T(i * (i + 1)) for i in 1:nmax]
                dv2 = -dv1
            else
                dv1 = [T(0.5) * √T(i * (i + 1)) for i in 1:nmax]
                dv2 = dv1
            end

            return dv1, dv2
        end
    end

    dv1, dv2 = vig(nmax, m, x)
    dv1 /= √(one(x) - x^2)

    return dv1, dv2
end

@doc raw"""
```
sphericalbesselj!(x::T, nmax::Int64, nnmax1::Int64, y::AbstractArray{T}, u::AbstractArray{T}, z::AbstractArray{T}) where {T <: Number}
```

Calculate spherical Bessel function $j_n(x)$ and $\frac{1}{x}\frac{\mathrm{d}}{\mathrm{d}x}[xj_n(x)]$ in place.
"""
function sphericalbesselj!(x::T,
                           nmax::Int64,
                           nnmax1::Int64,
                           y::AbstractArray{T},
                           u::AbstractArray{T},
                           z::AbstractArray{T}) where {T <: Number}
    l = nmax + nnmax1
    x1 = one(x) / x
    if length(z) < l
        resize!(z, l)
    end
    z[l] = x / (2l + 1)
    for i in (l - 1):-1:1
        z[i] = one(x) / ((2i + 1) * x1 - z[i + 1])
    end
    z0 = one(x) / (x1 - z[1])
    y0 = z0 * cos(x) * x1
    y[1] = y0 * z[1]
    u[1] = y0 - y[1] * x1
    for i in 2:nmax
        y[i] = y[i - 1] * z[i]
        u[i] = y[i - 1] - y[i] * x1 * i
    end

    return
end

function sphericalbesselj(x::T, nmax::Int64, nnmax1::Int64) where {T <: Number}
    y = zeros(T, nmax)
    u = zeros(T, nmax)
    z = zeros(T, nmax + nnmax1)
    sphericalbesselj!(x, nmax, nnmax1, y, u, z)
    return y, u
end

function sphericalbessely!(x::T, nmax::Int64, y::AbstractArray{T},
                           v::AbstractArray{T}) where {T <: Number}
    x1 = one(x) / x
    y[1] = -cos(x) * x1^2 - sin(x) * x1
    y[2] = (-3x1^3 + x1) * cos(x) - 3x1^2 * sin(x)
    for i in 2:(nmax - 1)
        y[i + 1] = (2i + 1) * x1 * y[i] - y[i - 1]
    end
    v[1] = -x1 * (cos(x) + y[1])
    for i in 2:nmax
        v[i] = y[i - 1] - i * x1 * y[i]
    end
end

function sphericalbessely(x::T, nmax::Int64) where {T <: Number}
    y = zeros(T, nmax)
    v = zeros(T, nmax)
    sphericalbessely!(x, nmax, y, v)
    return y, v
end

function cross_section(TT::Vector{<:AbstractMatrix}, λ::Real)
    nmax = length(TT) - 1

    Qsca = zero(λ)
    Qext = zero(λ)
    for n2 in 1:nmax
        nn2 = n2 + nmax
        for n1 in 1:nmax
            nn1 = n1 + nmax
            Qsca += TT[1][n1, n2] * TT[1][n1, n2]' +
                    TT[1][n1, nn2] * TT[1][n1, nn2]' +
                    TT[1][nn1, n2] * TT[1][nn1, n2]' +
                    TT[1][nn1, nn2] * TT[1][nn1, nn2]'
        end
    end
    for n in 1:(2nmax)
        Qext += real(TT[1][n, n])
    end

    for mm in 1:nmax
        nm = nmax - mm + 1
        for n2 in 1:nm
            nn2 = n2 + nm
            for n1 in 1:nm
                nn1 = n1 + nm
                Qsca += (TT[mm + 1][n1, n2] * TT[mm + 1][n1, n2]' +
                         TT[mm + 1][n1, nn2] * TT[mm + 1][n1, nn2]' +
                         TT[mm + 1][nn1, n2] * TT[mm + 1][nn1, n2]' +
                         TT[mm + 1][nn1, nn2] * TT[mm + 1][nn1, nn2]') * 2.0
            end
        end

        for n in 1:(2nm)
            Qext += real(TT[mm + 1][n, n]) * 2.0
        end
    end

    coeff = λ^2 / 2 / π
    Csca = abs(real(Qsca)) * coeff
    Cext = abs(Qext) * coeff
    ω = Csca / Cext

    if ω > 1.0
        @warn "ω is greater than 1.0"
    end

    return Csca, Cext, ω
end

gausslegendre(::Type{Float64}, n::Integer) = FastGaussQuadrature.gausslegendre(n)

function gausslegendre(T::Type{<:Real}, n::Integer)
    prec = precision(T)
    level = -((prec - 1) ÷ 64 + 1) * 16
    check = 10.0^level

    z = zeros(T, n)
    w = zeros(T, n)
    k = n ÷ 2 + n % 2
    for i in 1:k
        m = n + 1 - i
        if i == 1
            x = 1 - 2 / ((n + 1) * n)
        elseif i == 2
            x = (z[n] - 1) * 4 + z[n]
        elseif i == 3
            x = (z[n - 1] - z[n]) * 8 // 5 + z[n - 1]
        else
            x = (z[m + 1] - z[m + 2]) * 3 + z[m + 3]
        end
        if i == k && isodd(n)
            x = 0
        end

        check = T(check)
        pa = zero(T)
        pb = one(T)
        it = 0
        while abs(pb) > check * abs(x)
            it += 1
            if it > 100
                check *= 10
            end
            pa = zero(T)
            pb = one(T)
            pc = x
            for j in 2:n
                pa = pb
                pb = pc
                pc = x * pb + (x * pb - pa) * (j - 1) / j
            end
            pa = 1 / ((pb - x * pc) * n)
            pb = pa * pc * (1 - x * x)
            x -= pb
        end

        z[m] = x
        w[m] = 2 * pa * pa * (1 - x * x)
        if i != k || iseven(n)
            z[i] = -z[m]
            w[i] = w[m]
        end
    end

    return z, w
end

function gausslegendre(T::Type{<:Union{Arb, ArbRef}}, n::Integer)
    x = ArbRefVector(n)
    w = ArbRefVector(n)

    gausslegendre!(T, n, x, w)

    return x, w
end

function gausslegendre!(::Type{<:Union{Arb, ArbRef}}, n::Integer, x::Arblib.ArbVectorLike,
                        w::Arblib.ArbVectorLike)
    for i in 1:(n ÷ 2)
        Arblib.hypgeom_legendre_p_ui_root!(x[n + 1 - i], w[n + 1 - i], UInt64(n),
                                           UInt64(i - 1))
        x[i] = -x[n + 1 - i]
        w[i] = w[n + 1 - i]
    end

    if n % 2 == 1
        Arblib.hypgeom_legendre_p_ui_root!(x[n ÷ 2 + 1], w[n ÷ 2 + 1], UInt64(n),
                                           UInt64(n ÷ 2))
    end
end

function ccg(T::Type{<:Real}, n::Int64, n1::Int64, nmax::Int64, k1::Int64, k2::Int64)
    m_l = (k1 == 1 && k2 == 0) ? 0 : -n

    return [(m >= m_l && nn >= max(m * k1 + k2, n - n1)) ?
            clebschgordan(T, n, m, n1, m * (k1 - 1) + k2, nn) : zero(T)
            for
            m in (-n):n, nn in 0:min(n + n1, nmax)]
end
