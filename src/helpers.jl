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

    a = 1.0
    qs = √(1.0 - x * x)
    qs1 = 1.0 / qs
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

function sphericalbesselj!(x::T, nmax::Int64, nnmax1::Int64, y::AbstractArray, u::AbstractArray) where {T<:Number}
    l = nmax + nnmax1
    z = zeros(T, l)
    x1 = 1.0 / x
    z[l] = x / (2l + 1)
    for i in (l - 1):-1:1
        z[i] = 1.0 / ((2i + 1) * x1 - z[i + 1])
    end
    z0 = 1.0 / (x1 - z[1])
    y0 = z0 * cos(x) * x1
    y[1] = y0 * z[1]
    u[1] = y0 - y[1] * x1
    for i in 2:nmax
        y[i] = y[i - 1] * z[i]
        u[i] = y[i - 1] - y[i] * x1 * i
    end
end

function sphericalbesselj(x::T, nmax::Int64, nnmax1::Int64) where {T<:Number}
    y = zeros(T, nmax)
    u = zeros(T, nmax)
    sphericalbesselj!(x, nmax, nnmax1, y, u)
    return y, u
end

function sphericalbessely!(x::T, nmax::Int64, y::AbstractArray, v::AbstractArray) where {T<:Real}
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
end

function sphericalbessely(x::T, nmax::Int64) where {T<:Real}
    y = zeros(T, nmax)
    v = zeros(T, nmax)
    sphericalbessely!(x, nmax, y, v)
    return y, v
end

function cross_section(T::Vector{Matrix{Complex{P}}}, λ::P) where {P<:Real}
    nmax = length(T) - 1

    Qsca = 0.0
    Qext = 0.0
    for n2 in 1:nmax
        nn2 = n2 + nmax
        for n1 in 1:nmax
            nn1 = n1 + nmax
            Qsca +=
                T[1][n1, n2] * T[1][n1, n2]' +
                T[1][n1, nn2] * T[1][n1, nn2]' +
                T[1][nn1, n2] * T[1][nn1, n2]' +
                T[1][nn1, nn2] * T[1][nn1, nn2]'
        end
    end
    for n in 1:(2nmax)
        Qext += real(T[1][n, n])
    end

    for mm in 1:nmax
        nm = nmax - mm + 1
        for n2 in 1:nm
            nn2 = n2 + nm
            for n1 in 1:nm
                nn1 = n1 + nm
                Qsca +=
                    (
                        T[mm + 1][n1, n2] * T[mm + 1][n1, n2]' +
                        T[mm + 1][n1, nn2] * T[mm + 1][n1, nn2]' +
                        T[mm + 1][nn1, n2] * T[mm + 1][nn1, n2]' +
                        T[mm + 1][nn1, nn2] * T[mm + 1][nn1, nn2]'
                    ) * 2.0
            end
        end

        for n in 1:(2nm)
            Qext += real(T[mm + 1][n, n]) * 2.0
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
