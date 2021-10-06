function theta_split!(
    scatterer::AbstractScatterer{T},
    ngauss::Int64,
    x::AbstractArray,
    w::AbstractArray,
) where {T<:Real}
    x0, w0 = gausslegendre(T, ngauss)
    x .= x0
    w .= w0
    return
end

function theta_split!(
    scatterer::AbstractScatterer{Arb},
    ngauss::Int64,
    x::Arblib.ArbVectorLike,
    w::Arblib.ArbVectorLike,
)
    gausslegendre!(Arb, ngauss, x, w)
    return
end
