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

    # Need to recalculate `x`, `w`, `s`, `r`, `dr`, `kr⁻¹` and `kₛr⁻¹` if `ngauss` changes.
    k = 2 * Arb(π) / scatterer.λ

    if ngauss != info.ngauss
        info.x = ArbRefVector(ngauss)
        info.w = ArbRefVector(ngauss)
        info.s = ArbRefVector(ngauss)
        info.r = ArbRefVector(ngauss)
        info.dr = ArbRefVector(ngauss)
        info.drr = ArbRefVector(ngauss)
        info.wr² = ArbRefVector(ngauss)
        info.kr⁻¹ = ArbRefVector(ngauss)
        info.kₛr⁻¹ = AcbRefVector(ngauss)

        calc_r!(scatterer, ngauss, info.x, info.w, info.r, info.dr)
    end

    kr = ArbRefVector(ngauss)
    Arblib.mul!(kr, info.r, k)
    kₛr = AcbRefVector(ngauss)
    Arblib.mul!(kₛr, AcbRefVector(kr), Acb(scatterer.m))

    if ngauss != info.ngauss
        for i in 1:ngauss
            Arblib.div!(info.drr[i], info.dr[i], info.r[i])

            Arblib.sqr!(info.wr²[i], info.r[i])
            Arblib.mul!(info.wr²[i], info.wr²[i], info.w[i])

            info.s[i] = info.x[i]
            Arblib.acos!(info.s[i], info.s[i])
            Arblib.sin!(info.s[i], info.s[i])
            Arblib.div!(info.s[i], ARB_ONE, info.s[i])

            Arblib.div!(info.kr⁻¹[i], ARB_ONE, kr[i])
            Arblib.div!(info.kₛr⁻¹[i], ACB_ONE, kₛr[i])
        end
    end

    # The rest need to be recalculated when either `ngauss` or `nmax` changes.

    info.jkr = ArbRefMatrix(ngauss, nmax)
    info.djkr = ArbRefMatrix(ngauss, nmax)
    info.ykr = ArbRefMatrix(ngauss, nmax)
    info.dykr = ArbRefMatrix(ngauss, nmax)
    info.hkr = AcbRefMatrix(ngauss, nmax)
    info.dhkr = AcbRefMatrix(ngauss, nmax)
    info.jkₛr = AcbRefMatrix(ngauss, nmax)
    info.djkₛr = AcbRefMatrix(ngauss, nmax)

    Threads.@threads for i in 1:ngauss
        x1 = info.kr⁻¹[i]
        x_s1 = info.kₛr⁻¹[i]
        j0 = Arb(0)
        y0 = Arb(0)
        j_s0 = Acb(0)
        cr = √(x1 * π / 2)
        cc = √(x_s1 * π / 2)

        Arblib.hypgeom_bessel_jy!(j0, y0, ARB_HALF, kr[i])
        Arblib.hypgeom_bessel_j!(j_s0, ACB_HALF, kₛr[i])

        for n in 1:nmax
            Arblib.hypgeom_bessel_jy!(info.jkr[i, n], info.ykr[i, n], ARB_HALF + n, kr[i])
            Arblib.hypgeom_bessel_j!(info.jkₛr[i, n], ACB_HALF + n, kₛr[i])
        end

        info.djkr[i, 1] = j0 - x1 * info.jkr[i, 1]
        info.dykr[i, 1] = y0 - x1 * info.ykr[i, 1]
        info.djkₛr[i, 1] = j_s0 - x_s1 * info.jkₛr[i, 1]
        for n in 2:nmax
            info.djkr[i, n] = info.jkr[i, n - 1] - n * x1 * info.jkr[i, n]
            info.dykr[i, n] = info.ykr[i, n - 1] - n * x1 * info.ykr[i, n]
            info.djkₛr[i, n] = info.jkₛr[i, n - 1] - n * x_s1 * info.jkₛr[i, n]
        end

        for n in 1:nmax
            info.jkr[i, n] *= cr
            info.djkr[i, n] *= cr
            info.ykr[i, n] *= cr
            info.dykr[i, n] *= cr
            info.jkₛr[i, n] *= cc
            info.djkₛr[i, n] *= cc
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
    @assert iseven(ngauss)

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
    kr⁻¹ = info.kr⁻¹
    kₛr⁻¹ = info.kₛr⁻¹
    wr² = info.wr²

    jkr = info.jkr
    djkr = info.djkr
    hkr = info.hkr
    dhkr = info.dhkr
    jkₛr = info.jkₛr
    djkₛr = info.djkₛr

    J₁₂ = AcbRefMatrix(nmax, nmax)
    J₂₁ = AcbRefMatrix(nmax, nmax)
    RgJ₁₂ = AcbRefMatrix(nmax, nmax)
    RgJ₂₁ = AcbRefMatrix(nmax, nmax)

    Threads.@threads for nn in 0:(nmax * nmax - 1)
        n₂ = nn ÷ nmax + 1
        n₁ = nn % nmax + 1
        if !(sym && (n₁ + n₂) % 2 == 1)
            for i in 1:ngss
                τ₁τ₂ = τ[i][n₁] * τ[i][n₂]
                d₁τ₂ = d[i][n₁] * τ[i][n₂]
                d₂τ₁ = d[i][n₂] * τ[i][n₁]

                J₁₂[n₁, n₂] +=
                    wr²[i] * jkₛr[i, n₂] * (dhkr[i, n₁] * τ₁τ₂ + drr[i] * an[n₁] * hkr[i, n₁] * kr⁻¹[i] * d₁τ₂)

                J₂₁[n₁, n₂] +=
                    wr²[i] * hkr[i, n₁] * (djkₛr[i, n₂] * τ₁τ₂ + drr[i] * an[n₂] * jkₛr[i, n₂] * kₛr⁻¹[i] * d₂τ₁)

                RgJ₁₂[n₁, n₂] +=
                    wr²[i] * jkₛr[i, n₂] * (djkr[i, n₁] * τ₁τ₂ + drr[i] * an[n₁] * jkr[i, n₁] * kr⁻¹[i] * d₁τ₂)

                RgJ₂₁[n₁, n₂] +=
                    wr²[i] * jkr[i, n₁] * (djkₛr[i, n₂] * τ₁τ₂ + drr[i] * an[n₂] * jkₛr[i, n₂] * kₛr⁻¹[i] * d₂τ₁)
            end
        end
    end

    J₁₂ .*= -1.0im * ann
    J₂₁ .*= 1.0im * ann

    RgJ₁₂ .*= -1.0im * ann
    RgJ₂₁ .*= 1.0im * ann

    k = 2 * Arb(π) / scatterer.λ
    kₛ = k * scatterer.m
    k² = Acb(k^2)
    kkₛ = Acb(k * kₛ)

    Q₁₁ = AcbRefMatrix(nmax, nmax)
    Q₁₂ = AcbRefMatrix(nmax, nmax)
    Q₂₁ = AcbRefMatrix(nmax, nmax)
    Q₂₂ = AcbRefMatrix(nmax, nmax)
    RgQ₁₁ = AcbRefMatrix(nmax, nmax)
    RgQ₁₂ = AcbRefMatrix(nmax, nmax)
    RgQ₂₁ = AcbRefMatrix(nmax, nmax)
    RgQ₂₂ = AcbRefMatrix(nmax, nmax)

    # Since T = -RgQ⋅Q', the coefficient -i of Q and RgQ can be cancelled out.

    tmp = AcbRefMatrix(nmax, nmax)
    Arblib.mul!(Q₁₁, J₂₁, kkₛ)
    Arblib.mul!(tmp, J₁₂, k²)
    Arblib.add!(Q₁₁, Q₁₁, tmp)

    Arblib.mul!(Q₂₂, J₁₂, kkₛ)
    Arblib.mul!(tmp, J₂₁, k²)
    Arblib.add!(Q₂₂, Q₂₂, tmp)

    Arblib.mul!(RgQ₁₁, RgJ₂₁, kkₛ)
    Arblib.mul!(tmp, RgJ₁₂, k²)
    Arblib.add!(RgQ₁₁, RgQ₁₁, tmp)

    Arblib.mul!(RgQ₂₂, RgJ₁₂, kkₛ)
    Arblib.mul!(tmp, RgJ₂₁, k²)
    Arblib.add!(RgQ₂₂, RgQ₂₂, tmp)

    Q = vcat(hcat(Q₁₁, Q₁₂), hcat(Q₂₁, Q₂₂))
    RgQ = vcat(hcat(RgQ₁₁, RgQ₁₂), hcat(RgQ₂₁, RgQ₂₂))

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
                push!(ACCURACY_INFO[], ("J₁₂", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(J₁₂)))
                push!(ACCURACY_INFO[], ("J₂₁", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(J₂₁)))
                push!(ACCURACY_INFO[], ("RgJ₁₂", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ₁₂)))
                push!(ACCURACY_INFO[], ("RgJ₂₁", 0, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ₂₁)))
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
    @assert iseven(ngauss)

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
    kr⁻¹ = info.kr⁻¹
    kₛr⁻¹ = info.kₛr⁻¹
    wr² = info.wr²

    jkr = info.jkr
    djkr = info.djkr
    hkr = info.hkr
    dhkr = info.dhkr
    jkₛr = info.jkₛr
    djkₛr = info.djkₛr

    nm = nmax - mm + 1
    J₁₁ = AcbRefMatrix(nm, nm)
    J₁₂ = AcbRefMatrix(nm, nm)
    J₂₁ = AcbRefMatrix(nm, nm)
    J₂₂ = AcbRefMatrix(nm, nm)
    RgJ₁₁ = AcbRefMatrix(nm, nm)
    RgJ₁₂ = AcbRefMatrix(nm, nm)
    RgJ₂₁ = AcbRefMatrix(nm, nm)
    RgJ₂₂ = AcbRefMatrix(nm, nm)

    Threads.@threads for nn in 0:(nm * nm - 1)
        n₂ = nn ÷ nm + mm
        n₁ = nn % nm + mm
        nn₂ = n₂ - mm + 1
        nn₁ = n₁ - mm + 1
        if !(sym && (n₁ + n₂) % 2 == 0)
            for i in 1:ngss
                pττp = p[i][n₁] * τ[i][n₂] + p[i][n₂] * τ[i][n₁]
                p₁d₂ = p[i][n₁] * d[i][n₂]

                J₁₁[nn₁, nn₂] += wr²[i] * hkr[i, n₁] * jkₛr[i, n₂] * pττp

                J₂₂[nn₁, nn₂] +=
                    wr²[i] * (
                        dhkr[i, n₁] * djkₛr[i, n₂] * pττp +
                        drr[i] *
                        (an[n₁] * hkr[i, n₁] * kr⁻¹[i] * djkₛr[i, n₂] + an[n₂] * jkₛr[i, n₂] * kₛr⁻¹[i] * dhkr[i, n₁]) *
                        p₁d₂
                    )

                RgJ₁₁[nn₁, nn₂] += wr²[i] * jkr[i, n₁] * jkₛr[i, n₂] * pττp

                RgJ₂₂[nn₁, nn₂] +=
                    wr²[i] * (
                        djkr[i, n₁] * djkₛr[i, n₂] * pττp +
                        drr[i] *
                        (an[n₁] * jkr[i, n₁] * kr⁻¹[i] * djkₛr[i, n₂] + an[n₂] * jkₛr[i, n₂] * kₛr⁻¹[i] * djkr[i, n₁]) *
                        p₁d₂
                    )
            end
        end

        if !(sym && (n₁ + n₂) % 2 == 1)
            for i in 1:ngss
                ppττ = p[i][n₁] * p[i][n₂] + τ[i][n₁] * τ[i][n₂]
                d₁τ₂ = d[i][n₁] * τ[i][n₂]
                d₂τ₁ = d[i][n₂] * τ[i][n₁]

                J₁₂[nn₁, nn₂] +=
                    wr²[i] * jkₛr[i, n₂] * (dhkr[i, n₁] * ppττ + drr[i] * an[n₁] * hkr[i, n₁] * kr⁻¹[i] * d₁τ₂)

                J₂₁[nn₁, nn₂] +=
                    wr²[i] * hkr[i, n₁] * (djkₛr[i, n₂] * ppττ + drr[i] * an[n₂] * jkₛr[i, n₂] * kₛr⁻¹[i] * d₂τ₁)

                RgJ₁₂[nn₁, nn₂] +=
                    wr²[i] * jkₛr[i, n₂] * (djkr[i, n₁] * ppττ + drr[i] * an[n₁] * jkr[i, n₁] * kr⁻¹[i] * d₁τ₂)

                RgJ₂₁[nn₁, nn₂] +=
                    wr²[i] * jkr[i, n₁] * (djkₛr[i, n₂] * ppττ + drr[i] * an[n₂] * jkₛr[i, n₂] * kₛr⁻¹[i] * d₂τ₁)
            end
        end
    end

    ann = view(ann, mm:nmax, mm:nmax)
    J₁₁ .*= -ann
    J₁₂ .*= -1.0im * ann
    J₂₁ .*= 1.0im * ann
    J₂₂ .*= -ann

    RgJ₁₁ .*= -ann
    RgJ₁₂ .*= -1.0im * ann
    RgJ₂₁ .*= 1.0im * ann
    RgJ₂₂ .*= -ann

    k = 2 * Arb(π) / scatterer.λ
    kₛ = k * scatterer.m
    k² = Acb(k^2)
    kkₛ = Acb(k * kₛ)

    Q₁₁ = AcbRefMatrix(nm, nm)
    Q₁₂ = AcbRefMatrix(nm, nm)
    Q₂₁ = AcbRefMatrix(nm, nm)
    Q₂₂ = AcbRefMatrix(nm, nm)
    RgQ₁₁ = AcbRefMatrix(nm, nm)
    RgQ₁₂ = AcbRefMatrix(nm, nm)
    RgQ₂₁ = AcbRefMatrix(nm, nm)
    RgQ₂₂ = AcbRefMatrix(nm, nm)

    # Since T = -RgQ⋅Q', the coefficient -i of Q and RgQ can be cancelled out.

    tmp = AcbRefMatrix(nm, nm)
    Arblib.mul!(Q₁₁, J₂₁, kkₛ)
    Arblib.mul!(tmp, J₁₂, k²)
    Arblib.add!(Q₁₁, Q₁₁, tmp)

    Arblib.mul!(Q₁₂, J₁₁, kkₛ)
    Arblib.mul!(tmp, J₂₂, k²)
    Arblib.add!(Q₁₂, Q₁₂, tmp)

    Arblib.mul!(Q₂₁, J₂₂, kkₛ)
    Arblib.mul!(tmp, J₁₁, k²)
    Arblib.add!(Q₂₁, Q₂₁, tmp)

    Arblib.mul!(Q₂₂, J₁₂, kkₛ)
    Arblib.mul!(tmp, J₂₁, k²)
    Arblib.add!(Q₂₂, Q₂₂, tmp)

    Arblib.mul!(RgQ₁₁, RgJ₂₁, kkₛ)
    Arblib.mul!(tmp, RgJ₁₂, k²)
    Arblib.add!(RgQ₁₁, RgQ₁₁, tmp)

    Arblib.mul!(RgQ₁₂, RgJ₁₁, kkₛ)
    Arblib.mul!(tmp, RgJ₂₂, k²)
    Arblib.add!(RgQ₁₂, RgQ₁₂, tmp)

    Arblib.mul!(RgQ₂₁, RgJ₂₂, kkₛ)
    Arblib.mul!(tmp, RgJ₁₁, k²)
    Arblib.add!(RgQ₂₁, RgQ₂₁, tmp)

    Arblib.mul!(RgQ₂₂, RgJ₁₂, kkₛ)
    Arblib.mul!(tmp, RgJ₂₁, k²)
    Arblib.add!(RgQ₂₂, RgQ₂₂, tmp)

    Q = vcat(hcat(Q₁₁, Q₁₂), hcat(Q₂₁, Q₂₂))
    RgQ = vcat(hcat(RgQ₁₁, RgQ₁₂), hcat(RgQ₂₁, RgQ₂₂))

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
                push!(ACCURACY_INFO[], ("J₁₁", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(J₁₁)))
                push!(ACCURACY_INFO[], ("J₁₂", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(J₁₂)))
                push!(ACCURACY_INFO[], ("J₂₁", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(J₂₁)))
                push!(ACCURACY_INFO[], ("J₂₂", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(J₂₂)))
                push!(ACCURACY_INFO[], ("RgJ₁₁", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ₁₁)))
                push!(ACCURACY_INFO[], ("RgJ₁₂", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ₁₂)))
                push!(ACCURACY_INFO[], ("RgJ₂₁", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ₂₁)))
                push!(ACCURACY_INFO[], ("RgJ₂₂", mm, nmax, ngauss, precision(Arb), rel_accuracy_bits(RgJ₂₂)))
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
