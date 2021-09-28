import Base: abs2, floor, complex, AbstractFloat, Complex

import LinearAlgebra: inv

const ARF_PREC_EXACT = 9223372036854775807

abs2(x::Union{Acb,AcbRef}) = real(x) * real(x) + imag(x) * imag(x)
AbstractFloat(x::Union{Acb,AcbRef}) = AbstractFloat(Arblib.realref(x))
Complex{Arb}(x::Acb) = x
Complex{Arb}(x::Complex) = Acb(x.re, x.im)
complex(x::Arb, y::Arb) = Acb(x, y)
floor(x::Arb) = floor(Float64(x))
nonref(::Union{Arb,ArbRef}) = Arb
nonref(::Union{Acb,AcbRef}) = Acb
nonref(x) = typeof(x)

function rel_accuracy_bits(A::AbstractArray{<:Union{Arb,ArbRef,Acb,AcbRef}})
    return minimum(
        [Arblib.rel_accuracy_bits(a) for a in A if abs(Arblib.rel_accuracy_bits(a)) < ARF_PREC_EXACT],
        init = precision(Arb),
    )
end

function rel_accuracy_bits(A::Vector{<:Union{ArbRefVector,AcbRefVector}})
    return minimum(
        minimum(
            [Arblib.rel_accuracy_bits(a) for a in B if abs(Arblib.rel_accuracy_bits(a)) < ARF_PREC_EXACT],
            init = precision(Arb),
        ) for B in A
    )
end
