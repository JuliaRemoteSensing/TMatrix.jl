function vig!(nmax::Int64, m::Int64, x::T, dv1::AbstractArray, dv2::AbstractArray) where {T<:Real}
    if x < -1 || x > 1 || abs(1.0 - abs(x)) < eps(x)
        error("Constraint violated: x ∈ (-1, 1)")
    end

    if m < 0
        error("Constraint violated: m >= 0")
    end

    if nmax <= 0
        error("Constraint violated: nmax >= 1")
    end

    a = one(x)
    qs = √(one(x) - x * x)
    qs1 = one(x) / qs
    if m == 0
        d1 = one(x)
        d2 = x
        for i in 1:nmax
            d3 = (T(2i + 1) * x * d2 - i * d1) / (i + 1)
            der = qs1 * T(i * (i + 1) / (2i + 1)) * (d3 - d1)
            dv1[i] = d2
            dv2[i] = der
            d1, d2 = d2, d3
        end
    else
        for i in 1:m
            a *= √T((2i - 1) / 2i) * qs
        end
        d1 = zero(x)
        d2 = a
        for i in m:nmax
            qnm = √T(i^2 - m^2)
            qnm1 = √T((i + 1)^2 - m^2)
            d3 = (T(2i + 1) * x * d2 - qnm * d1) / qnm1
            der = qs1 * (-(i + 1) * qnm * d1 + i * qnm1 * d3) / (2i + 1)
            dv1[i] = d2
            dv2[i] = der
            d1, d2 = d2, d3
        end
    end
end

function vig(nmax::Int64, m::Int64, x::T) where {T<:Real}
    dv1 = zeros(T, nmax)
    dv2 = zeros(T, nmax)
    vig!(nmax, m, x, dv1, dv2)
    return dv1, dv2
end

function vigampl(nmax::Int64, m::Int64, x::T) where {T<:Real}
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
function sphericalbesselj!(
    x::T,
    nmax::Int64,
    nnmax1::Int64,
    y::AbstractArray{T},
    u::AbstractArray{T},
    z::AbstractArray{T},
) where {T<:Number}
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

@doc raw"""
```
sphericalbesselj!(x::T, nmax::Int64, y::AbstractArray{T}, u::AbstractArray{T}) where {T <: Union{Arb, Acb}}
```

For `Arb` and `Acb`, use `Arblib.hypegeom_bessel_j!` instead.
"""
function sphericalbesselj!(
    x::T,
    nmax::Int64,
    jkr::AbstractArray{T},
    djkr::AbstractArray{T},
) where {T<:Union{Arb,Acb}}
    x1 = 1 / x
    y0 = zero(x)
    half = T(1//2)
    coeff = √(T(π) / 2x)
    Arblib.hypgeom_bessel_j!(y0, half, x)
    for i in 1:nmax
        jkr[i] = Arblib.hypgeom_bessel_j!(jkr[i], T(i) + half, x)
    end

    djkr[1] = y0 - x1 * jkr[1]
    for n in 2:nmax
        djkr[n] = jkr[n - 1] - n * x1 * jkr[n]
    end

    jkr .*= coeff
    djkr .*= coeff

    return
end

function sphericalbesselj(x::T, nmax::Int64, nnmax1::Int64) where {T<:Number}
    y = zeros(T, nmax)
    u = zeros(T, nmax)
    z = zeros(T, nmax + nnmax1)
    sphericalbesselj!(x, nmax, nnmax1, y, u, z)
    return y, u
end

function sphericalbessely!(x::T, nmax::Int64, y::AbstractArray{T}, v::AbstractArray{T}) where {T<:Number}
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

@doc raw"""
```
sphericalbessely!(x::T, nmax::Int64, y::AbstractArray{T}, u::AbstractArray{T}) where {T <: Union{Arb, Acb}}
```

For `Arb` and `Acb`, use `Arblib.hypegeom_bessel_y!` instead.
"""
function sphericalbessely!(x::T, nmax::Int64, ykr::AbstractArray{T}, dykr::AbstractArray{T}) where {T<:Union{Arb,Acb}}
    x1 = 1 / x
    y0 = zero(x)
    half = T(1//2)
    coeff = √(T(π) / 2x)
    Arblib.hypgeom_bessel_y!(y0, half, x)
    for i in 1:nmax
        ykr[i] = Arblib.hypgeom_bessel_y!(ykr[i], T(i) + half, x)
    end

    dykr[1] = y0 - x1 * ykr[1]
    for n in 2:nmax
        dykr[n] = ykr[n - 1] - n * x1 * ykr[n]
    end

    ykr .*= coeff
    dykr .*= coeff

    return
end

function sphericalbessely(x::T, nmax::Int64) where {T<:Number}
    y = zeros(T, nmax)
    v = zeros(T, nmax)
    sphericalbessely!(x, nmax, y, v)
    return y, v
end

function cross_section(TT::Vector{<:AbstractMatrix{CT}}, λ::T) where {T<:Real, CT<:Number}
    nmax = length(TT) - 1

    Qsca = zero(λ)
    Qext = zero(λ)
    for n2 in 1:nmax
        nn2 = n2 + nmax
        for n1 in 1:nmax
            nn1 = n1 + nmax
            Qsca +=
                TT[1][n1, n2] * TT[1][n1, n2]' +
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
                Qsca +=
                    (
                        TT[mm + 1][n1, n2] * TT[mm + 1][n1, n2]' +
                        TT[mm + 1][n1, nn2] * TT[mm + 1][n1, nn2]' +
                        TT[mm + 1][nn1, n2] * TT[mm + 1][nn1, n2]' +
                        TT[mm + 1][nn1, nn2] * TT[mm + 1][nn1, nn2]'
                    ) * 2.0
            end
        end

        for n in 1:(2nm)
            Qext += real(TT[mm + 1][n, n]) * 2.0
        end
    end

    coeff = 0.5 * λ^2 / π
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
    x = ArbVector(n, prec = prec)
    w = ArbVector(n, prec = prec)
    for i in 1:(n ÷ 2)
        Arblib.hypgeom_legendre_p_ui_root!(ref(x, i), ref(w, i), UInt64(n), UInt64(i - 1), prec = prec)
    end
    for i in (n - n ÷ 2 + 1):n
        x[i] = -ref(x, n + 1 - i)
        w[i] = ref(w, n + 1 - i)
    end

    return [T(ref(x, i)) for i in 1:n], [T(ref(w, i)) for i in 1:n]
end

function gausslegendre(::Type{<:Union{Arb, ArbRef}}, n::Integer)
    x = ArbVector(n)
    w = ArbVector(n)
    for i in 1:(n ÷ 2)
        Arblib.hypgeom_legendre_p_ui_root!(ref(x, i), ref(w, i), UInt64(n), UInt64(i - 1))
    end
    for i in (n - n ÷ 2 + 1):n
        x[i] = -ref(x, n + 1 - i)
        w[i] = ref(w, n + 1 - i)
    end
    return x, w
end
