const ARB_ONE = Arb(1)
const ARB_ONE_NEG = Arb(-1)
const ARB_HALF = Arb(1 // 2)
const ACB_ONE = Acb(1)
const ACB_HALF = Acb(1 // 2)

@doc raw"""
Constructor of `ScattererInfo` for `Arb`. Pre-assignment is not used, since `SubArray` does not work harmoniously with `ArbMatrix`.
"""
function ScattererInfo(T::Type{Arb})
    return ScattererInfo(
        0,
        0,
        0,
        0,
        ArbRefVector(0),
        ArbRefMatrix(0, 0),
        ArbRefVector(0),
        ArbRefVector(0),
        ArbRefVector(0),
        ArbRefVector(0),
        ArbRefVector(0),
        ArbRefVector(0),
        ArbRefVector(0),
        ArbRefVector(0),
        ArbRefVector(0),
        AcbRefVector(0),
        ArbRefMatrix(0, 0),
        ArbRefMatrix(0, 0),
        ArbRefMatrix(0, 0),
        ArbRefMatrix(0, 0),
        ArbRefMatrix(0, 0),
        ArbRefVector(0),
        ArbRefMatrix(0, 0),
        ArbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefVector(0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
        AcbRefMatrix(0, 0),
    )
end

function theta_split!(
    scatterer::AbstractScatterer{Arb},
    ngauss::Int64,
    x::Arblib.ArbVectorLike,
    w::Arblib.ArbVectorLike,
)
    if typeof(scatterer) <: Cylinder
        ng = ngauss ÷ 2
        ng1 = ng ÷ 2
        ng2 = ng - ng1
        x1, w1 = gausslegendre(Arb, ng1)
        x2, w2 = gausslegendre(Arb, ng2)
        xx = -cos(atan(scatterer.d_to_h))
        x[1:ng1] .= 0.5(xx + 1.0) .* x1 .+ 0.5(xx - 1.0)
        w[1:ng1] .= 0.5(xx + 1.0) .* w1
        x[(ng1 + 1):ng] .= -0.5xx .* x2 .+ 0.5xx
        w[(ng1 + 1):ng] .= -0.5xx .* w2
        x[(ng + 1):ngauss] .= (-1.0) .* x[ng:-1:1]
        w[(ng + 1):ngauss] .= w[ng:-1:1]
    else
        gausslegendre!(Arb, ngauss, x, w)
    end
end

function calc_r!(
    scatterer::AbstractScatterer{Arb},
    ngauss::Int64,
    x::Arblib.ArbVectorLike,
    w::Arblib.ArbVectorLike,
    r::Arblib.ArbVectorLike,
    dr::Arblib.ArbVectorLike,
)
    theta_split!(scatterer, ngauss, x, w)
    rev = scatterer.rev

    if typeof(scatterer) <: Cylinder
        if ngauss % 2 != 0
            error("Constraint violated: ngauss should be even for cylinders")
        end

        e = scatterer.d_to_h
        h = rev * ∛(2 / (3e^2))
        d = h * e

        for i in 1:(ngauss ÷ 2)
            cosθ = abs(x[i])
            cos²θ = cosθ^2
            sin²θ = 1 - cos²θ
            sinθ = √sin²θ
            r₁ = h / cosθ
            r₂ = d / sinθ
            if r₁ < r₂
                r[i] = r₁
                dr[i] = h * sinθ / cos²θ
            else
                r[i] = r₂
                dr[i] = -d * cosθ / sin²θ
            end
            r[ngauss + 1 - i] = r[i]
            dr[ngauss + 1 - i] = dr[i]
            dr[i] = -dr[i]
        end

    elseif typeof(scatterer) <: Spheroid
        if ngauss % 2 != 0
            error("Constraint violated: ngauss should be even for spheroids")
        end

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
    else
        @assert typeof(scatterer) <: Chebyshev
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
end

function update!(scatterer::AbstractScatterer{Arb}, ngauss::Int64, nmax::Int64)
    info = scatterer.info

    # No need to recalculate if both `ngauss` and `nmax` remains the same.
    if ngauss == info.ngauss && nmax == info.nmax
        return
    end

    # Need to recalculate `an`, `ann` and `sig` if `nmax` changes.
    if nmax != info.nmax
        info.an = calc_an(Arb, nmax)
        info.ann = calc_ann(Arb, nmax)
        info.sig = calc_sig(Arb, nmax)
    end

    # Need to recalculate `x`, `w`, `s`, `r`, `dr`, `kr1` and `kr_s1` if `ngauss` changes.
    k = 2 * Arb(π) / scatterer.λ

    if ngauss != info.ngauss
        info.x = ArbRefVector(ngauss)
        info.w = ArbRefVector(ngauss)
        info.s = ArbRefVector(ngauss)
        info.r = ArbRefVector(ngauss)
        info.dr = ArbRefVector(ngauss)
        info.drr = ArbRefVector(ngauss)
        info.wr2 = ArbRefVector(ngauss)
        info.kr1 = ArbRefVector(ngauss)
        info.kr_s1 = AcbRefVector(ngauss)

        calc_r!(scatterer, ngauss, info.x, info.w, info.r, info.dr)
    end

    kr = ArbRefVector(ngauss)
    Arblib.mul!(kr, info.r, k)
    kr_s = AcbRefVector(ngauss)
    Arblib.mul!(kr_s, AcbRefVector(kr), scatterer.m)

    if ngauss != info.ngauss
        for i in 1:ngauss
            Arblib.div!(info.drr[i], info.dr[i], info.r[i])

            Arblib.sqr!(info.wr2[i], info.r[i])
            Arblib.mul!(info.wr2[i], info.wr2[i], info.w[i])

            info.s[i] = info.x[i]
            Arblib.acos!(info.s[i], info.s[i])
            Arblib.sin!(info.s[i], info.s[i])
            Arblib.div!(info.s[i], ARB_ONE, info.s[i])

            Arblib.div!(info.kr1[i], ARB_ONE, kr[i])
            Arblib.div!(info.kr_s1[i], ACB_ONE, kr_s[i])
        end
    end

    # The rest need to be recalculated when either `ngauss` or `nmax` changes.

    info.jkr = ArbRefMatrix(ngauss, nmax)
    info.djkr = ArbRefMatrix(ngauss, nmax)
    info.ykr = ArbRefMatrix(ngauss, nmax)
    info.dykr = ArbRefMatrix(ngauss, nmax)
    info.hkr = AcbRefMatrix(ngauss, nmax)
    info.dhkr = AcbRefMatrix(ngauss, nmax)
    info.jkr_s = AcbRefMatrix(ngauss, nmax)
    info.djkr_s = AcbRefMatrix(ngauss, nmax)

    Threads.@threads for i in 1:ngauss
        x1 = info.kr1[i]
        x_s1 = info.kr_s1[i]
        j0 = Arb(0)
        y0 = Arb(0)
        j_s0 = Acb(0)
        cr = √(x1 * π / 2)
        cc = √(x_s1 * π / 2)

        Arblib.hypgeom_bessel_jy!(j0, y0, ARB_HALF, kr[i])
        Arblib.hypgeom_bessel_j!(j_s0, ACB_HALF, kr_s[i])

        for n in 1:nmax
            Arblib.hypgeom_bessel_jy!(info.jkr[i, n], info.ykr[i, n], ARB_HALF + n, kr[i])
            Arblib.hypgeom_bessel_j!(info.jkr_s[i, n], ACB_HALF + n, kr_s[i])
        end

        info.djkr[i, 1] = j0 - x1 * info.jkr[i, 1]
        info.dykr[i, 1] = y0 - x1 * info.ykr[i, 1]
        info.djkr_s[i, 1] = j_s0 - x_s1 * info.jkr_s[i, 1]
        for n in 2:nmax
            info.djkr[i, n] = info.jkr[i, n - 1] - n * x1 * info.jkr[i, n]
            info.dykr[i, n] = info.ykr[i, n - 1] - n * x1 * info.ykr[i, n]
            info.djkr_s[i, n] = info.jkr_s[i, n - 1] - n * x_s1 * info.jkr_s[i, n]
        end

        for n in 1:nmax
            info.jkr[i, n] *= cr
            info.djkr[i, n] *= cr
            info.ykr[i, n] *= cr
            info.dykr[i, n] *= cr
            info.jkr_s[i, n] *= cc
            info.djkr_s[i, n] *= cc
        end

        for n in 1:nmax
            info.hkr[i, n] = Acb(info.jkr[i, n], info.ykr[i, n])
            info.dhkr[i, n] = Acb(info.djkr[i, n], info.dykr[i, n])
        end
    end

    info.ngauss = ngauss
    info.nmax = nmax
    info.ngcap = ngauss
    info.ncap = nmax

    return
end

function tmatr0!(scatterer::AbstractScatterer{Arb}, ngauss::Int64, nmax::Int64)
    sym = has_symmetric_plane(scatterer)
    update!(scatterer, ngauss, nmax)

    info = scatterer.info
    an = info.an
    ann = info.ann
    sig = info.sig
    x = info.x

    d = [ArbRefVector(nmax) for _ in 1:ngauss]
    τ = [ArbRefVector(nmax) for _ in 1:ngauss]
    Threads.@threads for i in (ngauss ÷ 2 + 1):ngauss
        ineg = ngauss + 1 - i
        vig!(nmax, 0, x[i], d[i], τ[i])
        @. d[ineg] = d[i] * sig
        @. τ[ineg] = τ[i] * sig * (-1)
    end

    ngss = sym ? (ngauss ÷ 2) : ngauss
    drr = info.drr
    kr1 = info.kr1
    kr_s1 = info.kr_s1
    wr2 = info.wr2

    jkr = info.jkr
    djkr = info.djkr
    hkr = info.hkr
    dhkr = info.dhkr
    jkr_s = info.jkr_s
    djkr_s = info.djkr_s

    J12 = AcbRefMatrix(nmax, nmax)
    J21 = AcbRefMatrix(nmax, nmax)
    RgJ12 = AcbRefMatrix(nmax, nmax)
    RgJ21 = AcbRefMatrix(nmax, nmax)

    Threads.@threads for nn in 0:(nmax * nmax - 1)
        n2 = nn ÷ nmax + 1
        n1 = nn % nmax + 1
        for i in 1:ngss
            τ₁τ₂ = τ[i][n1] * τ[i][n2]
            d₁τ₂ = d[i][n1] * τ[i][n2]
            d₂τ₁ = d[i][n2] * τ[i][n1]

            if !(sym && (n1 + n2) % 2 == 1)
                J12[n1, n2] +=
                    wr2[i] * jkr_s[i, n2] * (dhkr[i, n1] * τ₁τ₂ + drr[i] * an[n1] * hkr[i, n1] * kr1[i] * d₁τ₂)

                J21[n1, n2] +=
                    wr2[i] * hkr[i, n1] * (djkr_s[i, n2] * τ₁τ₂ + drr[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d₂τ₁)

                RgJ12[n1, n2] +=
                    wr2[i] * jkr_s[i, n2] * (djkr[i, n1] * τ₁τ₂ + drr[i] * an[n1] * jkr[i, n1] * kr1[i] * d₁τ₂)

                RgJ21[n1, n2] +=
                    wr2[i] * jkr[i, n1] * (djkr_s[i, n2] * τ₁τ₂ + drr[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d₂τ₁)
            end
        end
    end

    J12 .*= -1.0im * ann
    J21 .*= 1.0im * ann

    RgJ12 .*= -1.0im * ann
    RgJ21 .*= 1.0im * ann

    k = 2 * Arb(π) / scatterer.λ
    k_s = k * scatterer.m
    kk = Acb(k^2)
    kk_s = k * k_s

    Q11 = AcbRefMatrix(nmax, nmax)
    Q12 = AcbRefMatrix(nmax, nmax)
    Q21 = AcbRefMatrix(nmax, nmax)
    Q22 = AcbRefMatrix(nmax, nmax)
    RgQ11 = AcbRefMatrix(nmax, nmax)
    RgQ12 = AcbRefMatrix(nmax, nmax)
    RgQ21 = AcbRefMatrix(nmax, nmax)
    RgQ22 = AcbRefMatrix(nmax, nmax)

    # Since T = -RgQ⋅Q', the coefficient -i of Q and RgQ can be cancelled out.

    tmp = AcbRefMatrix(nmax, nmax)
    Arblib.mul!(Q11, J21, kk_s)
    Arblib.mul!(tmp, J12, kk)
    Arblib.add!(Q11, Q11, tmp)

    Arblib.mul!(Q22, J12, kk_s)
    Arblib.mul!(tmp, J21, kk)
    Arblib.add!(Q22, Q22, tmp)

    Arblib.mul!(RgQ11, RgJ21, kk_s)
    Arblib.mul!(tmp, RgJ12, kk)
    Arblib.add!(RgQ11, RgQ11, tmp)

    Arblib.mul!(RgQ22, RgJ12, kk_s)
    Arblib.mul!(tmp, RgJ21, kk)
    Arblib.add!(RgQ22, RgQ22, tmp)

    Q = vcat(hcat(Q11, Q12), hcat(Q21, Q22))
    RgQ = vcat(hcat(RgQ11, RgQ12), hcat(RgQ21, RgQ22))

    IQ = similar(Q)
    if ARB_APPROX_INV[]
        ret = Arblib.approx_inv!(IQ, Q)
        if ret == 0
            @warn Q_INVERSION_WARNING
        end
        T0 = -RgQ * IQ
    else
        ret = Arblib.inv!(IQ, Q)
        if ret == 0
            @warn Q_INVERSION_WARNING
        end
        T0 = -RgQ * IQ

        if COLLECT_ACCURACY_INFO[]
            @debug begin
                push!(ACCURACY_INFO[], ("d", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(d)))
                push!(ACCURACY_INFO[], ("τ", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(τ)))
                push!(ACCURACY_INFO[], ("J12", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(J12)))
                push!(ACCURACY_INFO[], ("J21", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(J21)))
                push!(ACCURACY_INFO[], ("RgJ12", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ12)))
                push!(ACCURACY_INFO[], ("RgJ21", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ21)))
                push!(ACCURACY_INFO[], ("RgQ", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgQ)))
                push!(ACCURACY_INFO[], ("Q", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(Q)))
                push!(ACCURACY_INFO[], ("inv(Q)", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(IQ)))
                push!(ACCURACY_INFO[], ("T", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(T0)))
                "Collecting accuracy info..."
            end
        end

        @debug "Accuracy of T is $(rel_accuracy_bits(T0))"
    end

    return T0, Q, RgQ
end

function tmatr!(scatterer::AbstractScatterer{Arb}, m::Int64, ngauss::Int64, nmax::Int64;)
    sym = has_symmetric_plane(scatterer)
    update!(scatterer, ngauss, nmax)

    info = scatterer.info
    mm = max(m, 1)
    an = info.an
    ann = info.ann
    sig = info.sig
    x = info.x
    s = info.s

    d = [ArbRefVector(nmax) for _ in 1:ngauss]
    p = [ArbRefVector(nmax) for _ in 1:ngauss]
    τ = [ArbRefVector(nmax) for _ in 1:ngauss]

    Threads.@threads for i in (ngauss ÷ 2 + 1):ngauss
        ineg = ngauss + 1 - i
        vig!(nmax, m, x[i], d[i], τ[i])
        @. p[i] = d[i] * s[i] * m
        @. d[ineg] = d[i] * sig
        @. τ[ineg] = τ[i] * sig * (-1)
        @. p[ineg] = p[i] * sig
    end

    ngss = sym ? (ngauss ÷ 2) : ngauss
    drr = info.drr
    kr1 = info.kr1
    kr_s1 = info.kr_s1
    wr2 = info.wr2

    jkr = info.jkr
    djkr = info.djkr
    hkr = info.hkr
    dhkr = info.dhkr
    jkr_s = info.jkr_s
    djkr_s = info.djkr_s

    nm = nmax - mm + 1
    J11 = AcbRefMatrix(nm, nm)
    J12 = AcbRefMatrix(nm, nm)
    J21 = AcbRefMatrix(nm, nm)
    J22 = AcbRefMatrix(nm, nm)
    RgJ11 = AcbRefMatrix(nm, nm)
    RgJ12 = AcbRefMatrix(nm, nm)
    RgJ21 = AcbRefMatrix(nm, nm)
    RgJ22 = AcbRefMatrix(nm, nm)

    Threads.@threads for nn in 0:(nm * nm - 1)
        n2 = nn ÷ nm + mm
        n1 = nn % nm + mm
        nn2 = n2 - mm + 1
        nn1 = n1 - mm + 1
        for i in 1:ngss
            if !(sym && (n1 + n2) % 2 == 0)
                pττp = p[i][n1] * τ[i][n2] + p[i][n2] * τ[i][n1]
                p₁d₂ = p[i][n1] * d[i][n2]

                J11[nn1, nn2] += wr2[i] * hkr[i, n1] * jkr_s[i, n2] * pττp

                J22[nn1, nn2] +=
                    wr2[i] * (
                        dhkr[i, n1] * djkr_s[i, n2] * pττp +
                        drr[i] *
                        (
                            an[n1] * hkr[i, n1] * kr1[i] * djkr_s[i, n2] +
                            an[n2] * jkr_s[i, n2] * kr_s1[i] * dhkr[i, n1]
                        ) *
                        p₁d₂
                    )

                RgJ11[nn1, nn2] += wr2[i] * jkr[i, n1] * jkr_s[i, n2] * pττp

                RgJ22[nn1, nn2] +=
                    wr2[i] * (
                        djkr[i, n1] * djkr_s[i, n2] * pττp +
                        drr[i] *
                        (
                            an[n1] * jkr[i, n1] * kr1[i] * djkr_s[i, n2] +
                            an[n2] * jkr_s[i, n2] * kr_s1[i] * djkr[i, n1]
                        ) *
                        p₁d₂
                    )
            end

            if !(sym && (n1 + n2) % 2 == 1)
                ppττ = p[i][n1] * p[i][n2] + τ[i][n1] * τ[i][n2]
                d₁τ₂ = d[i][n1] * τ[i][n2]
                d₂τ₁ = d[i][n2] * τ[i][n1]

                J12[nn1, nn2] +=
                    wr2[i] * jkr_s[i, n2] * (dhkr[i, n1] * ppττ + drr[i] * an[n1] * hkr[i, n1] * kr1[i] * d₁τ₂)

                J21[nn1, nn2] +=
                    wr2[i] * hkr[i, n1] * (djkr_s[i, n2] * ppττ + drr[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d₂τ₁)

                RgJ12[nn1, nn2] +=
                    wr2[i] * jkr_s[i, n2] * (djkr[i, n1] * ppττ + drr[i] * an[n1] * jkr[i, n1] * kr1[i] * d₁τ₂)

                RgJ21[nn1, nn2] +=
                    wr2[i] * jkr[i, n1] * (djkr_s[i, n2] * ppττ + drr[i] * an[n2] * jkr_s[i, n2] * kr_s1[i] * d₂τ₁)
            end
        end
    end

    ann = view(ann, mm:nmax, mm:nmax)
    J11 .*= -ann
    J12 .*= -1.0im * ann
    J21 .*= 1.0im * ann
    J22 .*= -ann

    RgJ11 .*= -ann
    RgJ12 .*= -1.0im * ann
    RgJ21 .*= 1.0im * ann
    RgJ22 .*= -ann

    k = 2 * Arb(π) / scatterer.λ
    k_s = k * scatterer.m
    kk = Acb(k^2)
    kk_s = k * k_s

    Q11 = AcbRefMatrix(nm, nm)
    Q12 = AcbRefMatrix(nm, nm)
    Q21 = AcbRefMatrix(nm, nm)
    Q22 = AcbRefMatrix(nm, nm)
    RgQ11 = AcbRefMatrix(nm, nm)
    RgQ12 = AcbRefMatrix(nm, nm)
    RgQ21 = AcbRefMatrix(nm, nm)
    RgQ22 = AcbRefMatrix(nm, nm)

    # Since T = -RgQ⋅Q', the coefficient -i of Q and RgQ can be cancelled out.

    tmp = AcbRefMatrix(nm, nm)
    Arblib.mul!(Q11, J21, kk_s)
    Arblib.mul!(tmp, J12, kk)
    Arblib.add!(Q11, Q11, tmp)

    Arblib.mul!(Q12, J11, kk_s)
    Arblib.mul!(tmp, J22, kk)
    Arblib.add!(Q12, Q12, tmp)

    Arblib.mul!(Q21, J22, kk_s)
    Arblib.mul!(tmp, J11, kk)
    Arblib.add!(Q21, Q21, tmp)

    Arblib.mul!(Q22, J12, kk_s)
    Arblib.mul!(tmp, J21, kk)
    Arblib.add!(Q22, Q22, tmp)

    Arblib.mul!(RgQ11, RgJ21, kk_s)
    Arblib.mul!(tmp, RgJ12, kk)
    Arblib.add!(RgQ11, RgQ11, tmp)

    Arblib.mul!(RgQ12, RgJ11, kk_s)
    Arblib.mul!(tmp, RgJ22, kk)
    Arblib.add!(RgQ12, RgQ12, tmp)

    Arblib.mul!(RgQ21, RgJ22, kk_s)
    Arblib.mul!(tmp, RgJ11, kk)
    Arblib.add!(RgQ21, RgQ21, tmp)

    Arblib.mul!(RgQ22, RgJ12, kk_s)
    Arblib.mul!(tmp, RgJ21, kk)
    Arblib.add!(RgQ22, RgQ22, tmp)

    Q = vcat(hcat(Q11, Q12), hcat(Q21, Q22))
    RgQ = vcat(hcat(RgQ11, RgQ12), hcat(RgQ21, RgQ22))

    IQ = similar(Q)
    if ARB_APPROX_INV[]
        ret = Arblib.approx_inv!(IQ, Q)
        if ret == 0
            @warn Q_INVERSION_WARNING
        end
        Tm = -RgQ * IQ
    else
        ret = Arblib.inv!(IQ, Q)
        if ret == 0
            @warn Q_INVERSION_WARNING
        end
        Tm = -RgQ * IQ

        if COLLECT_ACCURACY_INFO[]
            @debug begin
                push!(ACCURACY_INFO[], ("d", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(d)))
                push!(ACCURACY_INFO[], ("p", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(p)))
                push!(ACCURACY_INFO[], ("τ", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(τ)))
                push!(ACCURACY_INFO[], ("J11", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(J11)))
                push!(ACCURACY_INFO[], ("J12", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(J12)))
                push!(ACCURACY_INFO[], ("J21", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(J21)))
                push!(ACCURACY_INFO[], ("J22", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(J22)))
                push!(ACCURACY_INFO[], ("RgJ11", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ11)))
                push!(ACCURACY_INFO[], ("RgJ12", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ12)))
                push!(ACCURACY_INFO[], ("RgJ21", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ21)))
                push!(ACCURACY_INFO[], ("RgJ22", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ22)))
                push!(ACCURACY_INFO[], ("RgQ", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgQ)))
                push!(ACCURACY_INFO[], ("Q", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(Q)))
                push!(ACCURACY_INFO[], ("inv(Q)", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(IQ)))
                push!(ACCURACY_INFO[], ("T", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(Tm)))
                "Collecting accuracy info..."
            end

            @debug "Accuracy of T is $(rel_accuracy_bits(Tm))"
        end
    end

    return Tm, Q, RgQ
end
