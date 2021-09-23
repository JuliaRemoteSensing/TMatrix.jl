import Base: floor, complex, AbstractFloat, Complex
import LinearAlgebra: inv

const ARF_PREC_EXACT = 9223372036854775807

AbstractFloat(x::Acb) = AbstractFloat(Arblib.imagref(x))
Complex{Arb}(x::Acb) = x
complex(x::Arb, y::Arb) = Acb(x, y)
floor(x::Arb) = floor(Float64(x))

function LinearAlgebra.inv(A::StridedMatrix{<:Union{Arb,ArbRef}})
    return (B = ArbMatrix(size(A)...); Arblib.inv!(B, ArbRefMatrix(A)); B)
end
function LinearAlgebra.inv(A::StridedMatrix{<:Union{Acb,AcbRef,Complex{Arb},Complex{ArbRef}}})
    return (B = AcbMatrix(size(A)...); Arblib.inv!(B, AcbRefMatrix(A)); B)
end

function approx_inv(A::StridedMatrix{<:Union{Arb,ArbRef}})
    return (B = ArbMatrix(size(A)...); Arblib.approx_inv!(B, ArbRefMatrix(A)); B)
end
function approx_inv(A::StridedMatrix{<:Union{Acb,AcbRef,Complex{Arb},Complex{ArbRef}}})
    return (B = AcbMatrix(size(A)...); Arblib.approx_inv!(B, AcbRefMatrix(A)); B)
end

function rel_accuracy_bits(A::AbstractArray{<:Union{Arb,ArbRef,Acb,AcbRef}})
    return minimum(
        [Arblib.rel_accuracy_bits(a) for a in A if abs(Arblib.rel_accuracy_bits(a)) < ARF_PREC_EXACT],
        init = precision(Arb),
    )
end
