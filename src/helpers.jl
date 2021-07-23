function vig(nmax::Int64, m::Int64, x::T) where {T<:Real}
    if x < -1 || x > 1 || abs(1.0 - abs(x)) < eps(x)
        error("Constraint violated: x ∈ (-1, 1)")
    end

    if m < 0
        error("Constraint violated: m >= 0")
    end

    if nmax <= 0
        error("Constraint violated: nmax >= 1")
    end

    a = 1.0
    qs = √(1.0 - x * x)
    qs1 = 1.0 / qs
    dv1 = zeros(T, nmax)
    dv2 = zeros(T, nmax)
    if m == 0
        d1 = 1.0
        d2 = x
        for i in 1:nmax
            d3 = ((2i + 1) * x * d2 - i * d1) / (i + 1)
            der = qs1 * (i * (i + 1) / (2i + 1)) * (d3 - d1)
            dv1[i] = d2
            dv2[i] = der
            d1, d2 = d2, d3
        end
    else
        for i in 1:m
            a *= √((2i - 1) / 2i) * qs
        end
        d1 = 0.0
        d2 = a
        for i in m:nmax
            qnm = √(i^2 - m^2)
            qnm1 = √((i + 1)^2 - m^2)
            d3 = ((2i + 1) * x * d2 - qnm * d1) / qnm1
            der = qs1 * (-(i + 1) * qnm * d1 + i * qnm1 * d3) / (2i + 1)
            dv1[i] = d2
            dv2[i] = der
            d1, d2 = d2, d3
        end
    end

    return dv1, dv2
end

function vigampl(nmax::Int64, m::Int64, x::T) where {T<:Real}
    if abs(1.0 - abs(x)) < eps(x)
        if m != 1
            return zeros(T, nmax), zeros(T, nmax)
        else
            if x < 0.0
                dv1 = [(-1)^(i + 1) * 0.5 * √(i * (i + 1)) for i in 1:nmax]
                dv2 = -dv1
            else
                dv1 = [0.5 * √(i * (i + 1)) for i in 1:nmax]
                dv2 = dv1
            end

            return dv1, dv2
        end
    end

    dv1, dv2 = vig(nmax, m, x)
    dv1 /= √(1.0 - x^2)

    return dv1, dv2
end

function sphericalbesselj(x::T, nmax::Int64, nnmax1::Int64) where {T<:Number}
    l = nmax + nnmax1
    z = zeros(T, l)
    x1 = 1.0 / x
    z[l] = x / (2l + 1)
    for i in (l - 1):-1:1
        z[i] = 1.0 / ((2i + 1) * x1 - z[i + 1])
    end
    z0 = 1.0 / (x1 - z[1])
    y0 = z0 * cos(x) * x1
    y = zeros(T, nmax)
    u = zeros(T, nmax)
    y[1] = y0 * z[1]
    u[1] = y0 - y[1] * x1
    for i in 2:nmax
        y[i] = y[i - 1] * z[i]
        u[i] = y[i - 1] - y[i] * x1 * i
    end

    return y, u
end

function sphericalbessely(x::T, nmax::Int64) where {T<:Number}
    y = zeros(T, nmax)
    v = zeros(T, nmax)
    x1 = 1.0 / x
    y[1] = -cos(x) * x1^2 - sin(x) * x1
    y[2] = (-3.0 * x1^3 + x1) * cos(x) - 3.0x1^2 * sin(x)
    for i in 2:(nmax - 1)
        y[i + 1] = (2i + 1) * x1 * y[i] - y[i - 1]
    end
    v[1] = -x1 * (cos(x) + y[1])
    for i in 2:nmax
        v[i] = y[i - 1] - i * x1 * y[i]
    end

    return y, v
end
