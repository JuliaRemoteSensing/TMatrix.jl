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
                spheroid = TMatrix.Spheroid(; m = m, a_to_c = a_to_c, rev = rev, λ = λ)

                T = TMatrix.calc_tmatrix!(spheroid, Float64(ddelta), ndgs)
                S, Z = TMatrix.calc_SZ(spheroid,
                                       Float64(α),
                                       Float64(β),
                                       Float64(ϑ_i),
                                       Float64(ϑ_s),
                                       Float64(φ_i),
                                       Float64(φ_s),
                                       T)

                spheroid2 = TMatrix.Spheroid(Arb; m = m, a_to_c = a_to_c, rev = rev, λ = λ)
                Ta = TMatrix.calc_tmatrix!(spheroid2, Arb(ddelta), ndgs)
                Sa, Za = TMatrix.calc_SZ(spheroid2, Arb(α), Arb(β), Arb(ϑ_i), Arb(ϑ_s),
                                         Arb(φ_i), Arb(φ_s), Ta)

                all(isapprox.(T, Ta, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(S, Sa, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(Z, Za, rtol = RTOL, atol = ATOL))
            end
        end

        @testset "for cylinders with r_to_h = $r_to_h" for r_to_h in [1 // 4, 1 // 2, 1]
            @test begin
                cylinder = TMatrix.Cylinder(; m = m, r_to_h = r_to_h, rev = rev, λ = λ)

                T = TMatrix.calc_tmatrix!(cylinder, Float64(ddelta), ndgs)
                S, Z = TMatrix.calc_SZ(cylinder,
                                       Float64(α),
                                       Float64(β),
                                       Float64(ϑ_i),
                                       Float64(ϑ_s),
                                       Float64(φ_i),
                                       Float64(φ_s),
                                       T)

                cylinder2 = TMatrix.Cylinder(Arb; m = m, r_to_h = r_to_h, rev = rev, λ = λ)
                Ta = TMatrix.calc_tmatrix!(cylinder2, Arb(ddelta), ndgs)
                Sa, Za = TMatrix.calc_SZ(cylinder2, Arb(α), Arb(β), Arb(ϑ_i), Arb(ϑ_s),
                                         Arb(φ_i), Arb(φ_s), Ta)

                all(isapprox.(T, Ta, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(S, Sa, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(Z, Za, rtol = RTOL, atol = ATOL))
            end
        end

        @testset "for Chebyshev particles with ε = $ε and ncheb = $ncheb" for (ε, ncheb) in [(ε,
                                                                                              ncheb)
                                                                                             for ε in [
                                                                                                     -3 //
                                                                                                     20,
                                                                                                     1 //
                                                                                                     100,
                                                                                                     1 //
                                                                                                     10,
                                                                                                 ],
                                                                                                 ncheb in [
                                                                                                     2,
                                                                                                     3,
                                                                                                     4,
                                                                                                 ]]
            @test begin
                chebyshev = TMatrix.Chebyshev(; m = m, ε = ε, n = ncheb, rev = rev, λ = λ)

                T = TMatrix.calc_tmatrix!(chebyshev, Float64(ddelta), ndgs)
                S, Z = TMatrix.calc_SZ(chebyshev,
                                       Float64(α),
                                       Float64(β),
                                       Float64(ϑ_i),
                                       Float64(ϑ_s),
                                       Float64(φ_i),
                                       Float64(φ_s),
                                       T)

                chebyshev2 = TMatrix.Chebyshev(Arb; m = m, ε = ε, n = ncheb, rev = rev,
                                               λ = λ)
                Ta = TMatrix.calc_tmatrix!(chebyshev2, Arb(ddelta), ndgs)
                Sa, Za = TMatrix.calc_SZ(chebyshev2, Arb(α), Arb(β), Arb(ϑ_i), Arb(ϑ_s),
                                         Arb(φ_i), Arb(φ_s), Ta)

                all(isapprox.(T, Ta, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(S, Sa, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(Z, Za, rtol = RTOL, atol = ATOL))
            end
        end
    end
end
