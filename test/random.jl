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
                scatterer = TMatrix.Scatterer(
                    r = rev,
                    shape = TMatrix.SHAPE_SPHEROID,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = m,
                    axis_ratio = a_to_c,
                    λ = λ,
                )

                T = TMatrix.calc_tmatrix!(scatterer, ddelta, ndgs)
                Csca, _ = TMatrix.cross_section(T, scatterer.λ)
                nmax = length(T) - 1

                np = -1
                TMatrix.Wrapper.Random.set_T!(T)

                α₁, α₂, α₃, α₄, β₁, β₂ = TMatrix.calc_expansion_coefficients(T, Csca, scatterer.λ)
                α₁′, α₂′, α₃′, α₄′, β₁′, β₂′, lmax = TMatrix.Wrapper.Random.gsp(nmax, Csca, scatterer.λ)

                all(isapprox.(α₁[0:lmax], α₁′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(α₂[0:lmax], α₂′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(α₃[0:lmax], α₃′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(α₄[0:lmax], α₄′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(β₁[0:lmax], β₁′[0:lmax], rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(β₂[0:lmax], β₂′[0:lmax], rtol = RTOL, atol = ATOL))
            end
        end

        @testset "for cylinders with d_to_h = $d_to_h" for d_to_h in [0.5, 1.0, 2.0]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = rev,
                    shape = TMatrix.SHAPE_CYLINDER,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = m,
                    axis_ratio = d_to_h,
                    λ = λ,
                )

                T = TMatrix.calc_tmatrix!(scatterer, ddelta, ndgs)
                Csca, _ = TMatrix.cross_section(T, scatterer.λ)
                nmax = length(T) - 1

                np = -2
                TMatrix.Wrapper.Random.set_T!(T)

                α₁, α₂, α₃, α₄, β₁, β₂ = TMatrix.calc_expansion_coefficients(T, Csca, scatterer.λ)
                α₁′, α₂′, α₃′, α₄′, β₁′, β₂′, lmax = TMatrix.Wrapper.Random.gsp(nmax, Csca, scatterer.λ)

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
                scatterer = TMatrix.Scatterer(
                    r = rev,
                    shape = TMatrix.SHAPE_CHEBYSHEV,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = m,
                    axis_ratio = ε,
                    n = ncheb,
                    λ = λ,
                )

                T = TMatrix.calc_tmatrix!(scatterer, ddelta, ndgs)
                Csca, _ = TMatrix.cross_section(T, scatterer.λ)
                nmax = length(T) - 1

                np = ncheb
                TMatrix.Wrapper.Random.set_T!(T)

                α₁, α₂, α₃, α₄, β₁, β₂ = TMatrix.calc_expansion_coefficients(T, Csca, scatterer.λ)
                α₁′, α₂′, α₃′, α₄′, β₁′, β₂′, lmax = TMatrix.Wrapper.Random.gsp(nmax, Csca, scatterer.λ)

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
