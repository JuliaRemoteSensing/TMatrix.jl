var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API","title":"API","text":"CurrentModule = TMatrix","category":"page"},{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [TMatrix]","category":"page"},{"location":"api/#TMatrix.AbstractScatterer","page":"API","title":"TMatrix.AbstractScatterer","text":"Abstract type for all scatterers.\n\n\n\n\n\n","category":"type"},{"location":"api/#TMatrix.Bicone","page":"API","title":"TMatrix.Bicone","text":"Construct a bicone.\n\n\n\n\n\n","category":"type"},{"location":"api/#TMatrix.Bicone-2","page":"API","title":"TMatrix.Bicone","text":"A bicone scatterer.\n\nAttributes:\n\nm: The complex refractive index.\nr: The radius of the cone.\nh: The height of one cone.\nλ: The wavelength of the incident wave.\ninfo: The accompanied information.\n\n\n\n\n\n","category":"type"},{"location":"api/#TMatrix.Capsule","page":"API","title":"TMatrix.Capsule","text":"Construct a capsule.\n\n\n\n\n\n","category":"type"},{"location":"api/#TMatrix.Capsule-2","page":"API","title":"TMatrix.Capsule","text":"A capsule scatterer, which is a cylinder with half-spherical caps on both ends.\n\nAttributes:\n\nm: The complex refractive index.\nr: The radius of the cylinder base and the hemisphere.\nh: The height of the cylinder.\nλ: The wavelength of the incident wave.\ninfo: The accompanied information.\n\n\n\n\n\n","category":"type"},{"location":"api/#TMatrix.Chebyshev","page":"API","title":"TMatrix.Chebyshev","text":"A Chebyshev scatterer defined by\n\n$r(\\theta, \\phi)=r_0(1+\\varepsilon T_n(\\cos\\theta))$\n\nin which T_n(costheta)=cos ntheta.\n\nAttributes:\n\nm: The complex refractive index.\nr₀: The radius of the base sphere.\nε: The deformation parameter, which satisfies -1levarepsilon1.\nn: The degree of the Chebyshev polynomial.\nλ: The wavelength of the incident wave.\ninfo: The accompanied information.\n\n\n\n\n\n","category":"type"},{"location":"api/#TMatrix.Cylinder","page":"API","title":"TMatrix.Cylinder","text":"Construct a cylinder. λ and m must be provided, and the size parameters are considered according to the following priority:\n\nr and h\nr and r_to_h\nh and r_to_h\nrev and r_to_h\nrea and r_to_h\nrmax and r_to_h\n\nIf none of the above is hit, an ArgumentError will be thrown.\n\n\n\n\n\n","category":"type"},{"location":"api/#TMatrix.Cylinder-2","page":"API","title":"TMatrix.Cylinder","text":"A cylinder scatterer.\n\nAttributes:\n\nm: The complex refractive index.\nr: The radius of the cylinder base.\nh: The height of the cylinder.\nλ: The wavelength of the incident wave.\ninfo: The accompanied information.\n\n\n\n\n\n","category":"type"},{"location":"api/#TMatrix.ScattererInfo","page":"API","title":"TMatrix.ScattererInfo","text":"Accompanied information of a scatterer.\n\n\n\n\n\n","category":"type"},{"location":"api/#TMatrix.ScattererInfo-Tuple{Type{<:Real}}","page":"API","title":"TMatrix.ScattererInfo","text":"Constructor of ScattererInfo for general data types. Space is pre-allocated to reduce allocations.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.ScattererInfo-Tuple{Type{Arblib.Arb}}","page":"API","title":"TMatrix.ScattererInfo","text":"Constructor of ScattererInfo for Arb. Pre-assignment is not used, since SubArray does not work harmoniously with ArbMatrix.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.Spheroid","page":"API","title":"TMatrix.Spheroid","text":"Construct a spheroid. λ and m must be provided, and the size parameters are considered according to the following priority:\n\na and c\na and a_to_c\nc and a_to_c\nrev and a_to_c\nrea and a_to_c\nrmax and a_to_c\n\nIf none of the above is hit, an ArgumentError will be thrown.\n\n\n\n\n\n","category":"type"},{"location":"api/#TMatrix.Spheroid-2","page":"API","title":"TMatrix.Spheroid","text":"A spheroid scatterer.\n\nAttributes:\n\nm: The complex refractive index.\na: Length of the semi-major axis.\nc: Length of the semi-minor axis.\nλ: The wavelength of the incident wave.\ninfo: The accompanied information.\n\n\n\n\n\n","category":"type"},{"location":"api/#TMatrix.calc_SZ-Union{Tuple{T}, Tuple{TMatrix.AbstractScatterer, T, T, T, T, T, T, Vector{<:AbstractMatrix}}} where T<:Real","page":"API","title":"TMatrix.calc_SZ","text":"calc_SZ(scatterer::AbstractScatterer{T}, α::T, β::T, ϑ_i::T, ϑ_s::T, φ_i::T, φ_s::T, TT::Vector{<:AbstractMatrix}) where {T<:Real}\n\nCalculate the S matrix and the Z matrix sequentially.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.calc_Z","page":"API","title":"TMatrix.calc_Z","text":"calc_Z(S::AbstractMatrix)\n\nAlias of calc_phase\n\n\n\n\n\n","category":"function"},{"location":"api/#TMatrix.calc_amplitude-Union{Tuple{T}, Tuple{TMatrix.AbstractScatterer{T}, T, T, T, T, T, T, Vector{<:AbstractMatrix}}} where T<:Real","page":"API","title":"TMatrix.calc_amplitude","text":"calc_amplitude(scatterer::AbstractScatterer{T}, α::T, β::T, ϑ_i::T, ϑ_s::T, φ_i::T, φ_s::T, TT::Vector{<:AbstractMatrix}) where {T<:Real}\n\nCalculate the amplitude matrix and the phase matrix, given the scatterer and the geometry of the incident and the scattered beam. Use pre-computed T-Matrix when possible.\n\nParameters:\n\nscatterer: The scatterer.\nα, β: The Euler angle.\nϑ_i: The zenith angle of the incident beam.\nϑ_s: The zenith angle of the scattered beam.\nφ_i: The azimuth angle of the indicent beam.\nφ_s: The azimuth angle of the scatterer beam.\nTT: The pre-computed T-Matrix of the scatterer.\n\nAll the angles here are input in degrees.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.calc_expansion_coefficients-Tuple{Vector{<:AbstractMatrix}, Real, Real}","page":"API","title":"TMatrix.calc_expansion_coefficients","text":"calc_expansion_coefficients(TT::Vector{<:AbstractMatrix}, Csca::Real, λ::Real)\n\nCalculate the expansion coefficients from a given T-Matrix.\n\nParameters:\n\nTT: The precalculated T-Matrix of a scatterer.\nCsca: The scattering cross setction.\nλ: The wavelength.\n\nCsca and λ should have the same unit of length. E.g., if λ is in μm, Csca should be in μm².\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.calc_phase-Tuple{AbstractMatrix}","page":"API","title":"TMatrix.calc_phase","text":"calc_phase(S::AbstractMatrix)\n\nCalculate the phase matrix using the given amplitude matrix mathbfS.\n\nParameters:\n\nS, the precalculated amplitude matrix.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.calc_r!-Union{Tuple{T}, Tuple{TMatrix.Capsule{T}, Int64, AbstractArray, AbstractArray, AbstractArray, AbstractArray}} where T<:Real","page":"API","title":"TMatrix.calc_r!","text":"calc_r!(scatterer::Capsule{T}, ngauss::Int64, x::AbstractArray{T}, w::AbstractArray{T}, r::AbstractArray{T}, dr::AbstractArray{T}) where {T<:Real}\n\nCalculate r(theta) and fracmathrmdrmathrmdtheta at ngauss points for a given capsule, in place.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.calc_r!-Union{Tuple{T}, Tuple{TMatrix.Cylinder{T}, Int64, AbstractArray, AbstractArray, AbstractArray, AbstractArray}} where T<:Real","page":"API","title":"TMatrix.calc_r!","text":"calc_r!(scatterer::AbstractScatterer{T}, ngauss::Int64, x::AbstractArray{T}, w::AbstractArray{T}, r::AbstractArray{T}, dr::AbstractArray{T}) where {T<:Real}\n\nCalculate r(theta) and fracmathrmdrmathrmdtheta at ngauss points for a given scatterer, in place.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.calc_r-Union{Tuple{T}, Tuple{TMatrix.AbstractScatterer{T}, Int64}} where T<:Real","page":"API","title":"TMatrix.calc_r","text":"calc_r(scatterer::Scatterer, ngauss::Int64)\n\nCalculate r(theta) and fracmathrmdrmathrmdtheta at ngauss points for a given scatterer.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.calc_scattering_matrix-Tuple{TMatrix.AbstractScatterer, Vector{<:AbstractMatrix}, Integer}","page":"API","title":"TMatrix.calc_scattering_matrix","text":"calc_scattering_matrix(scatterer, TT, Nθ)\n\nCalculate the scatterering matrix elements from the given scatterer and precalculated T-Matrix.\n\nParameters:\n\nscatterer: The scatterer.\nTT: The T-Matrix.\nNθ: Number of θ intervals (so the result will have Nθ + 1 rows).\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.calc_scattering_matrix-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{T}, AbstractVector{T}, AbstractVector{T}, AbstractVector{T}, AbstractVector{T}, Real}} where T<:Real","page":"API","title":"TMatrix.calc_scattering_matrix","text":"calc_scattering_matrix(α₁, α₂, α₃, α₄, β₁, β₂, θ)\n\nCalculate the scatterering matrix elements from the given expansion coefficients.\n\nParameters:\n\nα₁, α₂, α₃, α₄, β₁, β₂: The precalculated expansion coefficients.\nθ: The scattering angle in degrees.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.calc_tmatrix!-Union{Tuple{TMatrix.AbstractScatterer{T}}, Tuple{T}} where T<:Real","page":"API","title":"TMatrix.calc_tmatrix!","text":"calc_tmatrix!(scatterer::AbstractScatterer{T}) where {T<:Real}\n\nCalculate the T-Matrix of the scatterer with default settings:\n\nUse tmatrix_routine_mishchenko to generate the routine function.\nUse ddelta = 0.001 and ndgs = 4.\n\nParameters:\n\nscatterer: The scatterer.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.calc_tmatrix!-Union{Tuple{T}, Tuple{TMatrix.AbstractScatterer{T}, Int64, Int64, Function}} where T<:Real","page":"API","title":"TMatrix.calc_tmatrix!","text":"calc_tmatrix(scatterer::AbstractScatterer{T}, ngstart::Int64, nstart::Int64, routine::Function) where {T<:Real}\n\nCalculate the T-Matrix of the scatterer.\n\nParameters:\n\nscatterer: The scatterer.\nngstart: The starting point of ngauss.\nnstart: The starting point of nmax.\nroutine: The iteration routine function generated by a routine generator, internally or customly implemented.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.hovenr-NTuple{6, Any}","page":"API","title":"TMatrix.hovenr","text":"hovenr(α₁, α₂, α₃, α₄, β₁, β₂)\n\nValidate the expansion coefficients with the test of Van der Mee & Hovenier.\n\nNote that the index of the coefficients should start from 0.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.set_arb_approx_inv","page":"API","title":"TMatrix.set_arb_approx_inv","text":"set_arb_approx_inv(t::Bool = true)\n\nTurn on/off Arblib's approx_inv for ArbMatrix and AcbMatrix.\n\n\n\n\n\n","category":"function"},{"location":"api/#TMatrix.set_default_ncap","page":"API","title":"TMatrix.set_default_ncap","text":"set_default_ncap(ncap::Int64 = 100)\n\nChange the default ncap value.\n\n\n\n\n\n","category":"function"},{"location":"api/#TMatrix.set_default_ngcap","page":"API","title":"TMatrix.set_default_ngcap","text":"set_default_ngcap(ngcap::Int64 = 500)\n\nChange the default ngcap value.\n\n\n\n\n\n","category":"function"},{"location":"api/#TMatrix.set_default_ngcheb","page":"API","title":"TMatrix.set_default_ngcheb","text":"set_default_ngcheb(ngcheb::Int64 = 60)\n\nChange the default ngcheb value.\n\n\n\n\n\n","category":"function"},{"location":"api/#TMatrix.sphericalbesselj!-Union{Tuple{T}, Tuple{T, Int64, Int64, AbstractArray{T}, AbstractArray{T}, AbstractArray{T}}} where T<:Number","page":"API","title":"TMatrix.sphericalbesselj!","text":"sphericalbesselj!(x::T, nmax::Int64, nnmax1::Int64, y::AbstractArray{T}, u::AbstractArray{T}, z::AbstractArray{T}) where {T <: Number}\n\nCalculate spherical Bessel function j_n(x) and frac1xfracmathrmdmathrmdxxj_n(x) in place.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.tmatrix_routine_mishchenko-Union{Tuple{T}, Tuple{TMatrix.AbstractScatterer{T}, T, Int64}} where T<:Real","page":"API","title":"TMatrix.tmatrix_routine_mishchenko","text":"tmatrix_routine_mishchenko(scatterer::AbstractScatterer{T}, ddelta::T, ndgs::Int64) where {T<:Real}\n\nGenerate the routine function which follows the iteration strategy used by M. Mishchenko.\n\nParameters:\n\nscatterer: The scatterer.\nddelta: The convergence threshold.\nndgs: Determines the initial value of ngauss as ndgs * nmax. The initial value of nmax is determined by the formula max(4 lfloor kr + 405 * sqrt3krrfloor).\nnstart: Optional for manual setting of the initial value of nmax.\nngstart: Optional for manual setting of the initial value of ngauss. This parameter takes effect only if nstart≠0.\n\n\n\n\n\n","category":"method"},{"location":"api/#TMatrix.tmatrix_routine_mishchenko_nmaxonly-Union{Tuple{T}, Tuple{TMatrix.AbstractScatterer{T}, T, Int64}} where T<:Real","page":"API","title":"TMatrix.tmatrix_routine_mishchenko_nmaxonly","text":"tmatrix_routine_mishchenko_nmaxonly(scatterer::AbstractScatterer{T}, ddelta::T, ndgs::Int64) where {T<:Real}\n\nGenerate the routine function which generally follows the iteration strategy used by M. Mishchenko, but does not modify the value of ngauss.\n\nParameters:\n\nscatterer: The scatterer.\nddelta: The convergence threshold.\nndgs: Determines the initial value of ngauss as ndgs * nmax. The initial value of nmax is determined by the formula max(4 lfloor kr + 405 * sqrt3krrfloor).\nnstart: Optional for manual setting of the initial value of nmax.\nngstart: Optional for manual setting of the initial value of ngauss. This parameter takes effect only if nstart≠0.\n\n\n\n\n\n","category":"method"},{"location":"#TMatrix","page":"Home","title":"TMatrix","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"TMatrix.jl is a Julia package for the calculation of the T-matrix of a dielectric particle, using the Extended Boundary Condition Method (EBCM), also known as the Null-Field Method (NFM). The T-matrix is a compact representation of the scattering properties of a particle, and is widely used.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg\njulia> Pkg.add(\"TMatrix\")","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"First, one needs to define a scatterer. Currently, the following shapes are supported:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Spheroids\nCylinders\nChebyshev particles\nCapsules, i.e., cylinders with half-spherical caps on both ends.\nBi-cones","category":"page"},{"location":"","page":"Home","title":"Home","text":"One can refer to the API for the detailed definitions of these shapes.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Take a spheroid as an example:","category":"page"},{"location":"","page":"Home","title":"Home","text":"scatterer = TMatrix.Spheroid(a = 1.0, c = 2.0, m = 1.311, λ = 2π)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Then one can calculate the T-matrix of the scatterer:","category":"page"},{"location":"","page":"Home","title":"Home","text":"T = TMatrix.tmatrix(scatterer)","category":"page"},{"location":"","page":"Home","title":"Home","text":"With the T-matrix calculated, one can calculate the scattering properties of the scatterer.","category":"page"},{"location":"","page":"Home","title":"Home","text":"For the fixed orientation case, the amplitude scattering matrix S and the phase matrix Z can be calculated as:","category":"page"},{"location":"","page":"Home","title":"Home","text":"α = 30.0\nβ = 60.0\nϑᵢ = 30.0\nϑₛ = 60.0\nφᵢ = 15.0\nφₛ = 90.0\nS, Z = TMatrix.calc_SZ(scatterer, α, β, ϑᵢ, ϑₛ, φᵢ, φₛ, T)","category":"page"},{"location":"","page":"Home","title":"Home","text":"For the random orientation case, one can first calculate the expansion coefficients of the T-matrix:","category":"page"},{"location":"","page":"Home","title":"Home","text":"coeffs = TMatrix.calc_expansion_coefficients(scatterer, T)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here the expansion coefficients are given in the form (α₁, α₂, α₃, α₄, β₁, β₂), where each element is a Vector.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Then the orientation-averaged scattering matrix can be evaluated at an arbitrary scattering angle:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Θ = 10.0\nF = TMatrix.calc_F(coeffs..., Θ)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here F is in the form of (F₁₁, F₁₂, F₂₂, F₃₃, F₃₄, F₄₄). The scattering matrix can be constructed from F as:","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginbmatrix\nF_11  F_12  0  0 \nF_12  F_22  0  0 \n0  0  F_33  F_34 \n0  0  -F_34  F_44\nendbmatrix","category":"page"}]
}
