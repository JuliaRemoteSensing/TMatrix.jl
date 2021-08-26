module Wrapper

using Libdl

const NPN1 = 100
const NPN2 = 2NPN1
const NPN3 = NPN1 + 1
const NPN4 = NPN1
const NPN5 = 2NPN4
const NPN6 = NPN4 + 1
const NPL = NPN2 + 1
const NPNG1 = 500
const NPNG2 = 2NPNG1
const TM_LIB = Ref{Ptr}()

function __init__()
    return TM_LIB[] = Libdl.dlopen(joinpath(@__DIR__, "..", "shared", "tmatrix.so"))
end

function sarea(e::Float64)
    ratio = zeros(1)
    ccall(Libdl.dlsym(TM_LIB[], :sarea_), Cvoid, (Ref{Float64}, Ptr{Float64}), e, ratio)
    return ratio[1]
end

function sareac(e::Float64)
    ratio = zeros(1)
    ccall(Libdl.dlsym(TM_LIB[], :sareac_), Cvoid, (Ref{Float64}, Ptr{Float64}), e, ratio)
    return ratio[1]
end

function surfch(n::Int64, e::Float64)
    ratio = zeros(1)
    ccall(Libdl.dlsym(TM_LIB[], :surfch_), Cvoid, (Ref{Int32}, Ref{Float64}, Ptr{Float64}), convert(Int32, n), e, ratio)
    return ratio[1]
end

function gauss(ngauss::Int64)
    z = zeros(ngauss)
    w = zeros(ngauss)
    ccall(
        Libdl.dlsym(TM_LIB[], :gauss_),
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

function rsp1(ngauss::Int, rev::Float64, e::Float64)
    r = zeros(ngauss)
    dr = zeros(ngauss)

    np = -1
    x, _ = const_(ngauss, 1, np, e)

    ccall(
        Libdl.dlsym(TM_LIB[], :rsp1_),
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

function rsp2(ngauss::Int64, rev::Float64, e::Float64, n::Int64)
    r = zeros(ngauss)
    dr = zeros(ngauss)

    np = n
    x, _ = const_(ngauss, 1, np, e)

    ccall(
        Libdl.dlsym(TM_LIB[], :rsp2_),
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

function rsp3(ngauss::Int64, rev::Float64, e::Float64)
    r = zeros(ngauss)
    dr = zeros(ngauss)

    np = -2
    x, _ = const_(ngauss, 1, np, e)

    ccall(
        Libdl.dlsym(TM_LIB[], :rsp3_),
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

function const_(ngauss::Int64, nmax::Int64, np::Int64, e::Float64)
    x = zeros(NPNG2)
    w = zeros(NPNG2)
    an = zeros(NPN1)
    ann = zeros(NPN1, NPN1)
    s = zeros(NPNG2)
    ss = zeros(NPNG2)
    ccall(
        Libdl.dlsym(TM_LIB[], :const_),
        Cvoid,
        (
            Ref{Int32},
            Ref{Int32},
            Ref{Int32},
            Ref{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ref{Int32},
            Ref{Float64},
        ),
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

function vary(
    x::Array{Float64},
    λ::Float64,
    m::ComplexF64,
    rev::Float64,
    e::Float64,
    np::Int64,
    ngauss::Int64,
    nmax::Int64,
)
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
        Libdl.dlsym(TM_LIB[], :vary_),
        Cvoid,
        (
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Int32},
            Ref{Int32},
            Ptr{Float64},
            Ref{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ref{Int32},
        ),
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

    cbess_ptr = cglobal(Libdl.dlsym(TM_LIB[], :cbess_), UInt64)
    blk_len = NPNG2 * NPN1 * 8
    jkr = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, cbess_ptr), (NPNG2, NPN1))[1:ngauss, 1:nmax]
    ykr = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, cbess_ptr + blk_len), (NPNG2, NPN1))[1:ngauss, 1:nmax]
    jkr_sr =
        unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, cbess_ptr + blk_len * 2), (NPNG2, NPN1))[1:ngauss, 1:nmax]
    jkr_si =
        unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, cbess_ptr + blk_len * 3), (NPNG2, NPN1))[1:ngauss, 1:nmax]
    djkr =
        unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, cbess_ptr + blk_len * 4), (NPNG2, NPN1))[1:ngauss, 1:nmax]
    dykr =
        unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, cbess_ptr + blk_len * 5), (NPNG2, NPN1))[1:ngauss, 1:nmax]
    djkr_sr =
        unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, cbess_ptr + blk_len * 6), (NPNG2, NPN1))[1:ngauss, 1:nmax]
    djkr_si =
        unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, cbess_ptr + blk_len * 7), (NPNG2, NPN1))[1:ngauss, 1:nmax]
    jkr_s = complex.(jkr_sr, jkr_si)
    djkr_s = complex.(djkr_sr, djkr_si)

    return ppi, pir, pii, r, dr, ddr, drr, dri, jkr, djkr, ykr, dykr, jkr_s, djkr_s
end

function vig(nmax::Int64, m::Int64, x::Float64)
    dv1 = zeros(nmax)
    dv2 = zeros(nmax)
    ccall(
        Libdl.dlsym(TM_LIB[], :vig_),
        Cvoid,
        (Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        x,
        convert(Int32, nmax),
        convert(Int32, m),
        dv1,
        dv2,
    )
    return dv1, dv2
end

function vigampl(nmax::Int64, m::Int64, x::Float64)
    dv1 = zeros(nmax)
    dv2 = zeros(nmax)
    ccall(
        Libdl.dlsym(TM_LIB[], :vigampl_),
        Cvoid,
        (Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        x,
        convert(Int32, nmax),
        convert(Int32, m),
        dv1,
        dv2,
    )
    return dv1, dv2
end

function rjb(x::Float64, n::Int64, nn::Int64)
    y = zeros(n)
    u = zeros(n)
    ccall(
        Libdl.dlsym(TM_LIB[], :rjb_),
        Cvoid,
        (Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ref{Int32}),
        x,
        y,
        u,
        convert(Int32, n),
        convert(Int32, nn),
    )

    return y, u
end

function ryb(x::Float64, n::Int64)
    y = zeros(n)
    u = zeros(n)
    ccall(
        Libdl.dlsym(TM_LIB[], :ryb_),
        Cvoid,
        (Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Int32}),
        x,
        y,
        u,
        convert(Int32, n),
    )

    return y, u
end

function cjb(x::ComplexF64, n::Int64, nn::Int64)
    xr = real(x)
    xi = imag(x)
    yr = zeros(n)
    yi = zeros(n)
    ur = zeros(n)
    ui = zeros(n)
    ccall(
        Libdl.dlsym(TM_LIB[], :cjb_),
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

function tmatr0(ngauss::Int64, nmax::Int64, np::Int64, e::Float64, λ::Float64, m::ComplexF64, rev::Float64)
    x, w, an, ann, s, ss = const_(ngauss, nmax, np, e)
    ppi, pir, pii, r, dr, ddr, drr, dri, _ = vary(x, λ, m, rev, e, np, ngauss, nmax)

    ncheck = zeros(Int32, 1)

    if np == -1 || np == -2 || (np > 0 && np % 2 == 0)
        ncheck[1] = 1
    end

    ccall(
        Libdl.dlsym(TM_LIB[], :tmatr0_),
        Cvoid,
        (
            Ref{Int32},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ref{Int32},
            Ptr{Int32},
        ),
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

    ctt_ptr = cglobal(Libdl.dlsym(TM_LIB[], :ctt_), UInt64)
    blk_len = NPN2 * NPN2 * 8
    QR = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr), (NPN2, NPN2))[1:(2nmax), 1:(2nmax)]
    QI = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr + blk_len), (NPN2, NPN2))[1:(2nmax), 1:(2nmax)]
    RgQR =
        unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr + blk_len * 2), (NPN2, NPN2))[1:(2nmax), 1:(2nmax)]
    RgQI =
        unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr + blk_len * 3), (NPN2, NPN2))[1:(2nmax), 1:(2nmax)]
    Q = complex.(QR, QI)
    RgQ = complex.(RgQR, RgQI)

    ct_ptr = cglobal(Libdl.dlsym(TM_LIB[], :ct_), UInt64)
    TR1 = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ct_ptr), (NPN2, NPN2))[1:(2nmax), 1:(2nmax)]
    TI1 = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ct_ptr + blk_len), (NPN2, NPN2))[1:(2nmax), 1:(2nmax)]

    T = complex.(TR1, TI1)

    return T, Q, RgQ
end

function tmatr(mm::Int64, ngauss::Int64, nmax::Int64, np::Int64, e::Float64, λ::Float64, m::ComplexF64, rev::Float64)
    x, w, an, ann, s, ss = const_(ngauss, nmax, np, e)
    ppi, pir, pii, r, dr, ddr, drr, dri, _ = vary(x, λ, m, rev, e, np, ngauss, nmax)

    ncheck = zeros(Int32, 1)

    if np == -1 || np == -2 || (np > 0 && np % 2 == 0)
        ncheck[1] = 1
    end

    ccall(
        Libdl.dlsym(TM_LIB[], :tmatr_),
        Cvoid,
        (
            Ref{Int32},
            Ref{Int32},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ref{Int32},
            Ptr{Int32},
        ),
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

    ctt_ptr = cglobal(Libdl.dlsym(TM_LIB[], :ctt_), UInt64)
    blk_len = NPN2 * NPN2 * 8
    nm = nmax - mm + 1
    QR = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr), (NPN2, NPN2))[1:(2nm), 1:(2nm)]
    QI = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr + blk_len), (NPN2, NPN2))[1:(2nm), 1:(2nm)]
    RgQR = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr + blk_len * 2), (NPN2, NPN2))[1:(2nm), 1:(2nm)]
    RgQI = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ctt_ptr + blk_len * 3), (NPN2, NPN2))[1:(2nm), 1:(2nm)]

    Q = complex.(QR, QI)
    RgQ = complex.(RgQR, RgQI)

    ct_ptr = cglobal(Libdl.dlsym(TM_LIB[], :ct_), UInt64)
    TR1 = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ct_ptr), (NPN2, NPN2))[1:(2nm), 1:(2nm)]
    TI1 = unsafe_wrap(Array{Float64,2}, convert(Ptr{Float64}, ct_ptr + blk_len), (NPN2, NPN2))[1:(2nm), 1:(2nm)]

    T = complex.(TR1, TI1)

    return T, Q, RgQ
end

function cross_section(ngauss::Int64, nmax::Int64, np::Int64, e::Float64, λ::Float64, m::ComplexF64, rev::Float64)
    T0, _ = tmatr0(ngauss, nmax, np, e, λ, m, rev)
    T = [T0]

    for mm in 1:nmax
        Tmm, _ = tmatr(mm, ngauss, nmax, np, e, λ, m, rev)
        push!(T, Tmm)
    end

    return TMatrix.cross_section(T, λ)
end

@doc raw"""
```
calc_tmatrix(axi::Float64, ratio::Float64, λ::Float64, m::ComplexF64, e::Float64, np::Int64, ddelt::Float64, ndgs::Int64)
```

Wrapper for the `CALCTMAT` function in Jussi Leinonen's modified version of `ampld.lp.f`.
"""
function calc_tmatrix(
    axi::Float64,
    ratio::Float64,
    λ::Float64,
    m::ComplexF64,
    e::Float64,
    np::Int64,
    ddelt::Float64,
    ndgs::Int64,
)
    nmax = zeros(Int32, 1)

    ccall(
        Libdl.dlsym(TM_LIB[], :calctmat_),
        Cvoid,
        (
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Int32},
            Ref{Float64},
            Ref{Int32},
            Ptr{Int32},
        ),
        axi,
        ratio,
        λ,
        real(m),
        imag(m),
        e,
        convert(Int32, np),
        ddelt,
        convert(Int32, ndgs),
        nmax,
    )

    tmat_ptr = cglobal(Libdl.dlsym(TM_LIB[], :tmat_), UInt64)
    blk_len = NPN6 * NPN4 * NPN4 * 4
    RT11 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr), (NPN6, NPN4, NPN4))
    RT12 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len), (NPN6, NPN4, NPN4))
    RT21 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len * 2), (NPN6, NPN4, NPN4))
    RT22 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len * 3), (NPN6, NPN4, NPN4))
    IT11 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len * 4), (NPN6, NPN4, NPN4))
    IT12 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len * 5), (NPN6, NPN4, NPN4))
    IT21 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len * 6), (NPN6, NPN4, NPN4))
    IT22 = unsafe_wrap(Array{Float32,3}, convert(Ptr{Float32}, tmat_ptr + blk_len * 7), (NPN6, NPN4, NPN4))

    nmax = Int64(nmax[1])

    T11 = complex.(RT11, IT11)
    T12 = complex.(RT12, IT12)
    T21 = complex.(RT21, IT21)
    T22 = complex.(RT22, IT22)

    T = []
    for mm in 0:nmax
        m1 = nmax - max(0, mm - 1)
        Tmm = zeros(ComplexF64, 2m1, 2m1)
        Tmm[1:m1, 1:m1] = T11[mm + 1, (nmax - m1 + 1):nmax, (nmax - m1 + 1):nmax]
        Tmm[1:m1, (m1 + 1):(2m1)] = T12[mm + 1, (nmax - m1 + 1):nmax, (nmax - m1 + 1):nmax]
        Tmm[(m1 + 1):(2m1), 1:m1] = T21[mm + 1, (nmax - m1 + 1):nmax, (nmax - m1 + 1):nmax]
        Tmm[(m1 + 1):(2m1), (m1 + 1):(2m1)] = T22[mm + 1, (nmax - m1 + 1):nmax, (nmax - m1 + 1):nmax]
        push!(T, Tmm)
    end

    return T, nmax
end

@doc raw"""
```
calc_amplitude(nmax::Int64, λ::Float64, α::Float64, β::Float64, ϑ_i::Float64, ϑ_s::Float64, φ_i::Float64, φ_s::Float64)
```

Wrapper for the `CALCAMPL` function in Jussi Leinonen's modified version of `ampld.lp.f`. Note that `calc_tmatrix` must be called first for this function to work.
"""
function calc_amplitude(
    nmax::Int64,
    λ::Float64,
    α::Float64,
    β::Float64,
    ϑ_i::Float64,
    ϑ_s::Float64,
    φ_i::Float64,
    φ_s::Float64,
)
    S = zeros(ComplexF64, 4)
    Z = zeros(16)

    ccall(
        Libdl.dlsym(TM_LIB[], :calcampl_),
        Cvoid,
        (
            Ref{Int32},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ref{Float64},
            Ptr{ComplexF64},
            Ptr{Float64},
        ),
        convert(Int32, nmax),
        λ,
        ϑ_i,
        ϑ_s,
        φ_i,
        φ_s,
        α,
        β,
        S,
        Z,
    )

    return reshape(S, (2, 2)), reshape(Z, (4, 4))
end

calc_SZ = calc_amplitude

end
