module Wrapper

using Libdl

const NPN1 = 100
const NPN2 = 200
const NPN3 = 101
const NPN4 = 100
const NPN5 = 200
const NPN6 = 101
const NPL = 201
const NPNG1 = 500
const NPNG2 = 1000

function sarea(tm, e::Float64)
    ratio = zeros(1)
    ccall(
        Libdl.dlsym(tm, :sarea_),
        Cvoid,
        (Ref{Float64}, Ptr{Float64}),
        e,
        ratio,
    )
    return ratio[1]
end

function sareac(tm, e::Float64)
    ratio = zeros(1)
    ccall(
        Libdl.dlsym(tm, :sareac_),
        Cvoid,
        (Ref{Float64}, Ptr{Float64}),
        e,
        ratio,
    )
    return ratio[1]
end

function surfch(tm, n::Int64, e::Float64)
    ratio = zeros(1)
    ccall(
        Libdl.dlsym(tm, :surfch_),
        Cvoid,
        (Ref{Int32}, Ref{Float64}, Ptr{Float64}),
        convert(Int32, n),
        e,
        ratio,
    )
    return ratio[1]
end

function gauss(tm, ngauss::Int64)
    z = zeros(ngauss)
    w = zeros(ngauss)
    ccall(
        Libdl.dlsym(tm, :gauss_),
        Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        convert(Int32, ngauss),
        0,
        0,
        z,
        w,
    )
    return z, w
end

function rsp1(tm, ngauss::Int, rev::Float64, e::Float64)
    r = zeros(ngauss)
    dr = zeros(ngauss)

    np = -1
    x, _ = const_(tm, ngauss, 1, np, e)

    ccall(
        Libdl.dlsym(tm, :rsp1_),
        Cvoid,
        (Ptr{Float64}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        x,
        convert(Int32, ngauss),
        convert(Int32, ngauss ÷ 2),
        rev,
        e,
        convert(Int32, np),
        r,
        dr,
    )
    return r, dr
end

function rsp2(tm, ngauss::Int64, rev::Float64, e::Float64, n::Int64)
    r = zeros(ngauss)
    dr = zeros(ngauss)

    np = n
    x, _ = const_(tm, ngauss, 1, np, e)

    ccall(
        Libdl.dlsym(tm, :rsp2_),
        Cvoid,
        (Ptr{Float64}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        x,
        convert(Int32, ngauss),
        rev,
        e,
        convert(Int32, n),
        r,
        dr,
    )
    return r, dr
end

function rsp3(tm, ngauss::Int64, rev::Float64, e::Float64)
    r = zeros(ngauss)
    dr = zeros(ngauss)

    np = -2
    x, _ = const_(tm, ngauss, 1, np, e)

    ccall(
        Libdl.dlsym(tm, :rsp3_),
        Cvoid,
        (Ptr{Float64}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}),
        x,
        convert(Int32, ngauss),
        convert(Int32, ngauss ÷ 2),
        rev,
        e,
        r,
        dr,
    )
    return r, dr
end

function const_(tm, ngauss::Int64, nmax::Int64, np::Int64, e::Float64)
    x = zeros(NPNG2)
    w = zeros(NPNG2)
    an = zeros(NPN1)
    ann = zeros(NPN1, NPN1)
    s = zeros(NPNG2)
    ss = zeros(NPNG2)
    ccall(
        Libdl.dlsym(tm, :const_),
        Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ref{Float64}),
        convert(Int32, ngauss ÷ 2),
        convert(Int32, nmax),
        convert(Int32, nmax),
        float(π),
        x,
        w,
        an,
        ann,
        s,
        ss,
        convert(Int32, np),
        e,
    )
    return x, w, an, ann, s, ss
end

function vary(tm, x::Array{Float64}, λ::Float64, m::ComplexF64, rev::Float64, e::Float64, np::Int64, ngauss::Int64, nmax::Int64)
    mrr = real(m)
    mri = imag(m)
    p = float(π)
    ppi = zeros(1)
    pir = zeros(1)
    pii = zeros(1)
    r = zeros(NPNG2)
    dr = zeros(NPNG2)
    ddr = zeros(NPNG2)
    drr = zeros(NPNG2)
    dri = zeros(NPNG2)
    ccall(
        Libdl.dlsym(tm, :vary_),
        Cvoid,
        (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Int32}),
        λ,
        mrr,
        mri,
        rev,
        e,
        convert(Int32, np),
        convert(Int32, ngauss ÷ 2),
        x,
        p,
        ppi,
        pir,
        pii,
        r,
        dr,
        ddr,
        drr,
        dri,
        convert(Int32, nmax),
    )

    return ppi, pir, pii, r, dr, ddr, drr, dri
end

function vig(tm, nmax::Int64, m::Int64, x::Float64)
    dv1 = zeros(nmax)
    dv2 = zeros(nmax)
    ccall(
        Libdl.dlsym(tm, :vig_),
        Cvoid,
        (Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        x,
        convert(Int32, nmax),
        convert(Int32, m),
        dv1,
        dv2,
    )
    dv1, dv2
end

function vigampl(tm, nmax::Int64, m::Int64, x::Float64)
    dv1 = zeros(nmax)
    dv2 = zeros(nmax)
    ccall(
        Libdl.dlsym(tm, :vigampl_),
        Cvoid,
        (Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        x,
        convert(Int32, nmax),
        convert(Int32, m),
        dv1,
        dv2,
    )
    dv1, dv2
end

function cjb(tm, x::Array{ComplexF64}, n::Int64, nn::Int64)
    xr = real.(x)
    xi = imag.(x)
    yr = zeros(n)
    yi = zeros(n)
    ur = zeros(n)
    ui = zeros(n)
    ccall(
        Libdl.dlsym(tm, :cjb_),
        Cvoid,
        (Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ref{Int32}),
        xr,
        xi,
        yr,
        yi,
        ur,
        ui,
        convert(Int32, n),
        convert(Int32, nn),
    )

    y = [yri + yii * 1.0im for (yri, yii) in zip(yr, yi)]
    u = [uri + uii * 1.0im for (uri, uii) in zip(ur, ui)]
    return y, u
end

function tmatr0(tm, ngauss::Int64, nmax::Int64, np::Int64, e::Float64, λ::Float64, m::ComplexF64, rev::Float64)
    x, w, an, ann, s, ss = const_(tm, ngauss, nmax, np, e)
    ppi, pir, pii, r, dr, ddr, drr, dri = vary(tm, x, λ, m, rev, e, np, ngauss, nmax)

    ncheck = zeros(Int32, 1)

    if np == -1 || np == -2 || (np > 0 && np % 2 == 0)
        ncheck[1] = 1
    end

    ccall(
        Libdl.dlsym(tm, :tmatr0_),
        Cvoid,
        (Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ptr{Int32}),
        convert(Int32, ngauss ÷ 2),
        x,
        w,
        an,
        ann,
        s,
        ss,
        ppi,
        pir,
        pii,
        r,
        dr,
        ddr,
        drr,
        dri,
        convert(Int32, nmax),
        ncheck,
    )

    ctt_ptr = cglobal(Libdl.dlsym(tm, :ctt_), UInt64)
    blk_len = NPN2 * NPN2 * 8
    qr = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr), (NPN2, NPN2))[1:2nmax,1:2nmax]
    qi = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr + blk_len), (NPN2, NPN2))[1:2nmax,1:2nmax]
    rgqr = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr + blk_len * 2), (NPN2, NPN2))[1:2nmax,1:2nmax]
    rgqi = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr + blk_len * 3), (NPN2, NPN2))[1:2nmax,1:2nmax]
    q = complex.(qr, qi)
    rgq = complex.(rgqr, rgqi)

    ct_ptr = cglobal(Libdl.dlsym(tm, :ct_), UInt64)
    tr1 = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ct_ptr), (NPN2, NPN2))[1:2nmax,1:2nmax]
    ti1 = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ct_ptr + blk_len), (NPN2, NPN2))[1:2nmax,1:2nmax]

    t = complex.(tr1, ti1)

    return t, q, rgq
end

function tmatr(tm, mm::Int64, ngauss::Int64, nmax::Int64, np::Int64, e::Float64, λ::Float64, m::ComplexF64, rev::Float64)
    x, w, an, ann, s, ss = const_(tm, ngauss, nmax, np, e)
    ppi, pir, pii, r, dr, ddr, drr, dri = vary(tm, x, λ, m, rev, e, np, ngauss, nmax)

    ncheck = zeros(Int32, 1)

    if np == -1 || np == -2 || (np > 0 && np % 2 == 0)
        ncheck[1] = 1
    end

    ccall(
        Libdl.dlsym(tm, :tmatr_),
        Cvoid,
        (Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ptr{Int32}),
        convert(Int32, mm),
        convert(Int32, ngauss ÷ 2),
        x,
        w,
        an,
        ann,
        s,
        ss,
        ppi,
        pir,
        pii,
        r,
        dr,
        ddr,
        drr,
        dri,
        convert(Int32, nmax),
        ncheck,
    )

    ctt_ptr = cglobal(Libdl.dlsym(tm, :ctt_), UInt64)
    blk_len = NPN2 * NPN2 * 8
    nm = nmax - mm + 1
    qr = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr), (NPN2, NPN2))[1:2nm,1:2nm]
    qi = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr + blk_len), (NPN2, NPN2))[1:2nm,1:2nm]
    rgqr = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr + blk_len * 2), (NPN2, NPN2))[1:2nm,1:2nm]
    rgqi = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr + blk_len * 3), (NPN2, NPN2))[1:2nm,1:2nm]

    q = complex.(qr, qi)
    rgq = complex.(rgqr, rgqi)

    ct_ptr = cglobal(Libdl.dlsym(tm, :ct_), UInt64)
    tr1 = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ct_ptr), (NPN2, NPN2))
    ti1 = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ct_ptr + blk_len), (NPN2, NPN2))

    t = complex.(tr1, ti1)

    return t, q, rgq
end

function cross_section(tm, ngauss::Int64, nmax::Int64, np::Int64, e::Float64, λ::Float64, m::ComplexF64, rev::Float64)
    Qext = 0.0
    Qsca = 0.0

    t0, _, _ = tmatr0(tm, ngauss, nmax, np, e, λ, m, rev)
    for n2 in 1:nmax
        nn2 = n2 + nmax
        for n1 in 1:nmax
            nn1 = n1 + nmax
            Qsca += t0[n1, n2] * t0[n1, n2]' + t0[n1, nn2] * t0[n1, nn2]' + t0[nn1, n2] * t0[nn1, n2]' + t0[nn1, nn2] * t0[nn1, nn2]'
        end
    end
    for n in 1:2nmax
        Qext += real(t0[n, n])
    end

    for mm in 1:nmax
        Qsc = 0.0
        Qex = 0.0

        tmm, _, _ = tmatr(tm, mm, ngauss, nmax, np, e, λ, m, rev)
        nm = nmax - mm + 1
        for n2 in 1:nm
            nn2 = n2 + nm
            for n1 in 1:nm
                nn1 = n1 + nm
                Qsca += (tmm[n1, n2] * tmm[n1, n2]' + tmm[n1, nn2] * tmm[n1, nn2]' + tmm[nn1, n2] * tmm[nn1, n2]' + tmm[nn1, nn2] * tmm[nn1, nn2]') * 2.0
            end
        end

        for n in 1:2nm
            Qext += real(tmm[n, n]) * 2.0
        end
    end

    coeff = 0.5λ^2 / π
    Csca = abs(real(Qsca)) * coeff
    Cext = abs(Qext) * coeff
    ω = Csca / Cext

    if ω > 1.0
        @warn "ω is greater than 1.0"
    end

    return Cext, Csca, ω
end

end
