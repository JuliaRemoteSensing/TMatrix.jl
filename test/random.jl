@testset "Random orientation" begin
    @testset "Calculate expansion coefficients" begin
        λ = 1.0
        m = 1.5 + 0.02im
        rev = 1.0
        ratio = 1.0
        ddelta = 0.001
        ndgs = 4

        @testset "for spheroids with a_to_c = $a_to_c" for a_to_c in [0.5, 1.01, 2.0]
            @test begin
                spheroid = TMatrix.Spheroid(m = m, a_to_c = a_to_c, rev = rev, λ = λ)

                T = TMatrix.calc_tmatrix!(spheroid, ddelta, ndgs)
                Csca, _ = TMatrix.cross_section(T, spheroid.λ)
                nmax = length(T) - 1

                np = -1
                TMatrix.Wrapper.Random.set_T!(T)

                α₁, α₂, α₃, α₄, β₁, β₂ = TMatrix.calc_expansion_coefficients(T, Csca, spheroid.λ)
                α₁′, α₂′, α₃′, α₄′, β₁′, β₂′, lmax = TMatrix.Wrapper.Random.gsp(nmax, Csca, spheroid.λ)

                all(isapprox.(α₁[0:lmax], α₁′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(α₂[0:lmax], α₂′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(α₃[0:lmax], α₃′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(α₄[0:lmax], α₄′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(β₁[0:lmax], β₁′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(β₂[0:lmax], β₂′[0:lmax], rtol = RTOL, atol = ATOL))
            end
        end

        @testset "for cylinders with r_to_h = $r_to_h" for r_to_h in [0.25, 0.5, 1.0]
            @test begin
                cylinder = TMatrix.Cylinder(m = m, r_to_h = r_to_h, rev = rev, λ = λ)
                T = TMatrix.calc_tmatrix!(cylinder, ddelta, ndgs)
                Csca, _ = TMatrix.cross_section(T, cylinder.λ)
                nmax = length(T) - 1

                np = -2
                TMatrix.Wrapper.Random.set_T!(T)

                α₁, α₂, α₃, α₄, β₁, β₂ = TMatrix.calc_expansion_coefficients(T, Csca, cylinder.λ)
                α₁′, α₂′, α₃′, α₄′, β₁′, β₂′, lmax = TMatrix.Wrapper.Random.gsp(nmax, Csca, cylinder.λ)

                all(isapprox.(α₁[0:lmax], α₁′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(α₂[0:lmax], α₂′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(α₃[0:lmax], α₃′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(α₄[0:lmax], α₄′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(β₁[0:lmax], β₁′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(β₂[0:lmax], β₂′[0:lmax], rtol = RTOL, atol = ATOL))
            end
        end

        @testset "for Chebyshev particles with ε = $ε and ncheb = $ncheb" for (ε, ncheb) in [
            (ε, ncheb) for ε in [-0.15, 0.01, 0.1], ncheb in [2, 3, 4]
        ]
            @test begin
                chebyshev = TMatrix.Chebyshev(m = m, ε = ε, n = ncheb, rev = rev, λ = λ)

                T = TMatrix.calc_tmatrix!(chebyshev, ddelta, ndgs)
                Csca, _ = TMatrix.cross_section(T, chebyshev.λ)
                nmax = length(T) - 1

                np = ncheb
                TMatrix.Wrapper.Random.set_T!(T)

                α₁, α₂, α₃, α₄, β₁, β₂ = TMatrix.calc_expansion_coefficients(T, Csca, chebyshev.λ)
                α₁′, α₂′, α₃′, α₄′, β₁′, β₂′, lmax = TMatrix.Wrapper.Random.gsp(nmax, Csca, chebyshev.λ)

                all(isapprox.(α₁[0:lmax], α₁′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(α₂[0:lmax], α₂′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(α₃[0:lmax], α₃′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(α₄[0:lmax], α₄′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(β₁[0:lmax], β₁′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(β₂[0:lmax], β₂′[0:lmax], rtol = RTOL, atol = ATOL))
            end
        end
    end
end
