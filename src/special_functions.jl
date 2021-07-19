function vig(nmax::Int64, m::Int64, x::T) where {T <: Real}
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

function vigampl(nmax::Int64, m::Int64, x::T) where {T <: Real}
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