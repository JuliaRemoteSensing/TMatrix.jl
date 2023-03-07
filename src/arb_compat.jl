import Base: abs2, floor, convert, complex, AbstractFloat, Complex

const ARF_PREC_EXACT = 9223372036854775807

abs2(x::Union{Acb, AcbRef}) = real(x) * real(x) + imag(x) * imag(x)
AbstractFloat(x::Union{Acb, AcbRef}) = AbstractFloat(Arblib.realref(x))
Complex{Arb}(x::Arblib.AcbLike) = x
Complex{Arb}(x::T) where {T <: Complex} = Acb(x.re, x.im)
complex(x::Arblib.ArbLike, y::Arblib.ArbLike) = Acb(x, y)
floor(x::Arblib.ArbLike) = floor(Float64(x))
nonref(::Union{Arb, ArbRef}) = Arb
nonref(::Union{Acb, AcbRef}) = Acb
nonref(x) = typeof(x)

function rel_accuracy_bits(A::AbstractArray{<:Union{Arb, ArbRef, Acb, AcbRef}})
    return minimum(Arblib.rel_accuracy_bits.(A); init = precision(Arb))
end

function rel_accuracy_bits(A::Vector{<:Union{ArbRefVector, AcbRefVector}})
    return minimum(minimum(Arblib.rel_accuracy_bits.(B); init = precision(Arb)) for B in A)
end
