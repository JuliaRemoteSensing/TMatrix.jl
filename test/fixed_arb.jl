TMatrix.set_arb_approx_inv(false)
TMatrix.collect_accuracy_info(false)
setprecision(Arb, 128)

@testset "Fixed orientation (for Arb)" begin
    @testset "Calculate amplitude" begin
        λ = 1
        mr = 3 // 2
        mi = 1 // 50
        m = complex(mr, mi)
        rev = 1
        ratio = 1
        ddelta = 1 // 1000
        ndgs = 4
        ϑ_i = 56
        ϑ_s = 65
        φ_i = 114
        φ_s = 128
        α = 145
        β = 52

        @testset "for spheroids with a_to_c = $a_to_c" for a_to_c in [1 // 2, 101 // 100, 2]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = rev,
                    shape = TMatrix.SHAPE_SPHEROID,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = m,
                    axis_ratio = a_to_c,
                    λ = λ,
                )

                T = TMatrix.calc_tmatrix!(scatterer, Float64(ddelta), ndgs)
                S, Z = TMatrix.calc_SZ(
                    scatterer,
                    Float64(α),
                    Float64(β),
                    Float64(ϑ_i),
                    Float64(ϑ_s),
                    Float64(φ_i),
                    Float64(φ_s),
                    T,
                )

                scatterer2 = TMatrix.Scatterer(
                    Arb,
                    r = rev,
                    shape = TMatrix.SHAPE_SPHEROID,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = m,
                    axis_ratio = a_to_c,
                    λ = λ,
                )
                Ta = TMatrix.calc_tmatrix!(scatterer2, Arb(ddelta), ndgs)
                Sa, Za = TMatrix.calc_SZ(scatterer2, Arb(α), Arb(β), Arb(ϑ_i), Arb(ϑ_s), Arb(φ_i), Arb(φ_s), Ta)

                all(isapprox.(T, Ta, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(S, Sa, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(Z, Za, rtol = RTOL, atol = ATOL))
            end
        end

        # FIXME: Currently, we have to use larger ndgs for cylinders so that Arb can converge.
        @testset "for cylinders with d_to_h = $d_to_h" for d_to_h in [1 // 2, 1, 2]
            local ndgs = 16

            @test begin
                scatterer = TMatrix.Scatterer(
                    r = rev,
                    shape = TMatrix.SHAPE_CYLINDER,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = m,
                    axis_ratio = d_to_h,
                    λ = λ,
                )

                T = TMatrix.calc_tmatrix!(scatterer, Float64(ddelta), ndgs)
                S, Z = TMatrix.calc_SZ(
                    scatterer,
                    Float64(α),
                    Float64(β),
                    Float64(ϑ_i),
                    Float64(ϑ_s),
                    Float64(φ_i),
                    Float64(φ_s),
                    T,
                )

                scatterer2 = TMatrix.Scatterer(
                    Arb,
                    r = rev,
                    shape = TMatrix.SHAPE_CYLINDER,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = m,
                    axis_ratio = d_to_h,
                    λ = λ,
                )
                Ta = TMatrix.calc_tmatrix!(scatterer2, Arb(ddelta), ndgs)
                Sa, Za = TMatrix.calc_SZ(scatterer2, Arb(α), Arb(β), Arb(ϑ_i), Arb(ϑ_s), Arb(φ_i), Arb(φ_s), Ta)

                all(isapprox.(T, Ta, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(S, Sa, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(Z, Za, rtol = RTOL, atol = ATOL))
            end
        end

        @testset "for Chebyshev particles with ε = $ε and ncheb = $ncheb" for (ε, ncheb) in [
            (ε, ncheb) for ε in [-3 // 20, 1 // 100, 1 // 10], ncheb in [2, 3, 4]
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

                T = TMatrix.calc_tmatrix!(scatterer, Float64(ddelta), ndgs)
                S, Z = TMatrix.calc_SZ(
                    scatterer,
                    Float64(α),
                    Float64(β),
                    Float64(ϑ_i),
                    Float64(ϑ_s),
                    Float64(φ_i),
                    Float64(φ_s),
                    T,
                )

                scatterer2 = TMatrix.Scatterer(
                    Arb,
                    r = rev,
                    shape = TMatrix.SHAPE_CHEBYSHEV,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = m,
                    axis_ratio = ε,
                    n = ncheb,
                    λ = λ,
                )
                Ta = TMatrix.calc_tmatrix!(scatterer2, Arb(ddelta), ndgs)
                Sa, Za = TMatrix.calc_SZ(scatterer2, Arb(α), Arb(β), Arb(ϑ_i), Arb(ϑ_s), Arb(φ_i), Arb(φ_s), Ta)

                all(isapprox.(T, Ta, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(S, Sa, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(Z, Za, rtol = RTOL, atol = ATOL))
            end
        end
    end
end
