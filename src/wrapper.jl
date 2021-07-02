module Wrapper

using Libdl

function gauss(tm, n::Int64)
    z = zeros(n)
    w = zeros(n)
    ccall(
        Libdl.dlsym(tm, :gauss_),
        Cvoid,
        (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        convert(Int32, n),
        0,
        0,
        z,
        w,
    )
    return z, w
end


function rsp1(tm, n::Int, rev::Float64, eps::Float64)
    ngauss = n ÷ 2
    r = zeros(n)
    dr = zeros(n)
    np = -1
    x, _ = gauss(tm, n)
    ccall(
        Libdl.dlsym(tm, :rsp1_),
        Cvoid,
        (Ptr{Float64}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        x,
        convert(Int32, n),
        convert(Int32, ngauss),
        rev,
        eps,
        convert(Int32, np),
        r,
        dr,
    )
    return r, dr
end

function rsp2(tm, nmax::Int64, rev::Float64, eps::Float64, n::Int64)
    r = zeros(nmax)
    dr = zeros(nmax)
    x, _ = gauss(tm, nmax)
    ccall(
        Libdl.dlsym(tm, :rsp2_),
        Cvoid,
        (Ptr{Float64}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}),
        x,
        convert(Int32, nmax),
        rev,
        eps,
        convert(Int32, n),
        r,
        dr,
    )
    return r, dr
end

function rsp3(tm, n::Int64, rev::Float64, eps::Float64)
    ngauss = n ÷ 2
    r = zeros(n)
    dr = zeros(n)
    x, _ = gauss(tm, n)
    ccall(
        Libdl.dlsym(tm, :rsp3_),
        Cvoid,
        (Ptr{Float64}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}),
        x,
        convert(Int32, n),
        convert(Int32, ngauss),
        rev,
        eps,
        r,
        dr,
    )
    return r, dr
end

end