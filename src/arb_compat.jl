import Base: floor, complex, AbstractFloat, Complex
import LinearAlgebra: inv

AbstractFloat(x::Acb) = AbstractFloat(Arblib.imagref(x))
Complex{Arb}(x::Acb) = x
complex(x::Arb, y::Arb) = Acb(x, y)
floor(x::Arb) = floor(Float64(x))

LinearAlgebra.inv(A::StridedMatrix{<:Union{Arb, ArbRef}}) = (B = ArbMatrix(size(A)...); Arblib.inv!(B, ArbRefMatrix(A)); B)
LinearAlgebra.inv(A::StridedMatrix{<:Union{Acb, AcbRef, Complex{Arb}, Complex{ArbRef}}}) = (B = AcbMatrix(size(A)...); Arblib.inv!(B, AcbRefMatrix(A)); B)

approx_inv(A::StridedMatrix{<:Union{Arb, ArbRef}}) = (B = ArbMatrix(size(A)...); Arblib.approx_inv!(B, ArbRefMatrix(A)); B)
approx_inv(A::StridedMatrix{<:Union{Acb, AcbRef, Complex{Arb}, Complex{ArbRef}}}) = (B = AcbMatrix(size(A)...); Arblib.approx_inv!(B, AcbRefMatrix(A)); B)
