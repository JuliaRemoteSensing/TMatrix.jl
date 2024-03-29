@testset "Fixed orientation" begin
    @testset "Calculate vig function" begin
        @testset "vig($nmax, $m, $x)" for (nmax, m, x) in [(nmax, m, x)
                                                           for nmax in [5, 10, 100],
                                                               m in [0, 1, 2, 3],
                                                               x in [
                                                                   -0.9,
                                                                   -0.5,
                                                                   -0.1,
                                                                   0.1,
                                                                   0.5,
                                                                   0.9,
                                                               ]]
            @test begin
                dv10, dv20 = TMatrix.Wrapper.Fixed.vig(nmax, m, x)
                dv1, dv2 = TMatrix.vig(nmax, m, x)
                all(dv1 .≈ dv10) && all(dv2 .≈ dv20)
            end
        end
    end

    @testset "Calculate vigampl function" begin
        @testset "vigampl($nmax, $m, $x)" for (nmax, m, x) in [(nmax, m, x)
                                                               for nmax in [5, 10, 100],
                                                                   m in [0, 1, 2, 3],
                                                                   x in [
                                                                       -1.0,
                                                                       -0.5,
                                                                       -0.1,
                                                                       0.1,
                                                                       0.5,
                                                                       1.0,
                                                                   ]]
            @test begin
                dv10, dv20 = TMatrix.Wrapper.Fixed.vigampl(nmax, m, x)
                dv1, dv2 = TMatrix.vigampl(nmax, m, x)
                all(dv1 .≈ dv10) && all(dv2 .≈ dv20)
            end
        end
    end

    @testset "Calculate r(θ) and dr/dθ" begin
        m = 1
        λ = 1

        @testset "for spheroids with rev = $rev, a_to_c = $a_to_c and ngauss = $ngauss" for (rev, a_to_c, ngauss) in [(rev,
                                                                                                                       a_to_c,
                                                                                                                       ngauss)
                                                                                                                      for rev in [
                                                                                                                              0.5,
                                                                                                                              1.0,
                                                                                                                              2.0,
                                                                                                                          ],
                                                                                                                          a_to_c in [
                                                                                                                              0.5,
                                                                                                                              1.0,
                                                                                                                              2.0,
                                                                                                                          ],
                                                                                                                          ngauss in [
                                                                                                                              4,
                                                                                                                              20,
                                                                                                                              1000,
                                                                                                                          ]]
            @test begin
                spheroid = TMatrix.Spheroid(; m = m, a_to_c = a_to_c, rev = rev, λ = λ)
                r, dr = TMatrix.calc_r(spheroid, ngauss)
                r², drr = TMatrix.Wrapper.Fixed.rsp1(ngauss, rev, a_to_c)

                all(r .^ 2 .≈ r²) && all(dr ./ r .≈ drr)
            end
        end

        @testset "for cylinders with rev = $rev, r_to_h = $r_to_h and ngauss = $ngauss" for (rev, r_to_h, ngauss) in [(rev,
                                                                                                                       r_to_h,
                                                                                                                       ngauss)
                                                                                                                      for rev in [
                                                                                                                              0.5,
                                                                                                                              1.0,
                                                                                                                              2.0,
                                                                                                                          ],
                                                                                                                          r_to_h in [
                                                                                                                              0.25,
                                                                                                                              0.5,
                                                                                                                              1.0,
                                                                                                                          ],
                                                                                                                          ngauss in [
                                                                                                                              4,
                                                                                                                              20,
                                                                                                                              1000,
                                                                                                                          ]]
            @test begin
                cylinder = TMatrix.Cylinder(; m = m, r_to_h = r_to_h, rev = rev, λ = λ)
                r, dr = TMatrix.calc_r(cylinder, ngauss)
                r², drr = TMatrix.Wrapper.Fixed.rsp3(ngauss, rev, 2r_to_h)

                all(r .^ 2 .≈ r²) && all(dr ./ r .≈ drr)
            end
        end

        @testset "for Chebyshev particles with rev = $rev, ε = $ε, ngauss = $ngauss and ncheb = $ncheb" for (rev,
        ε,
        ngauss,
        ncheb) in [(rev, ε, ngauss, ncheb)
                   for rev in [0.5, 1.0, 2.0], ε in [0.0, 0.1, 0.5],
                       ngauss in [4, 20, 1000],
                       ncheb in [2, 3, 4, 10]]
            @test begin
                chebyshev = TMatrix.Chebyshev(; m = m, ε = ε, rev = rev, n = ncheb, λ = λ)
                r, dr = TMatrix.calc_r(chebyshev, ngauss)
                r², drr = TMatrix.Wrapper.Fixed.rsp2(ngauss, rev, ε, ncheb)

                all(r .^ 2 .≈ r²) && all(dr ./ r .≈ drr)
            end
        end
    end

    @testset "Calculate rev from rea" begin
        m = 1
        λ = 1

        @testset "for spheroids with rea = $rea and a_to_c = $a_to_c" for (rea, a_to_c) in [(rea,
                                                                                             a_to_c)
                                                                                            for rea in [
                                                                                                    0.5,
                                                                                                    1.0,
                                                                                                    2.0,
                                                                                                ],
                                                                                                a_to_c in [
                                                                                                    0.5,
                                                                                                    1.000001,
                                                                                                    2.0,
                                                                                                ]]
            @test begin
                spheroid = TMatrix.Spheroid(; m = m, a_to_c = a_to_c, rea = rea, λ = λ)
                rev = TMatrix.volume_equivalent_radius(spheroid)
                ratio = TMatrix.Wrapper.Fixed.sarea(a_to_c)
                rev ≈ rea * ratio
            end
        end

        @testset "for cylinders with rea = $rea and r_to_h = $r_to_h" for (rea, r_to_h) in [(rea,
                                                                                             r_to_h)
                                                                                            for rea in [
                                                                                                    0.5,
                                                                                                    1.0,
                                                                                                    2.0,
                                                                                                ],
                                                                                                r_to_h in [
                                                                                                    0.25,
                                                                                                    0.5,
                                                                                                    1.0,
                                                                                                ]]
            @test begin
                cylinder = TMatrix.Cylinder(; m = m, r_to_h = r_to_h, rea = rea, λ = λ)
                rev = TMatrix.volume_equivalent_radius(cylinder)
                ratio = TMatrix.Wrapper.Fixed.sareac(2r_to_h)
                rev ≈ rea * ratio
            end
        end

        @testset "for Chebyshev particles with rea = $rea, ε = $ε and ncheb = $ncheb" for (rea, ε, ncheb) in [(rea,
                                                                                                               ε,
                                                                                                               ncheb)
                                                                                                              for rea in [
                                                                                                                      0.5,
                                                                                                                      1.0,
                                                                                                                      2.0,
                                                                                                                  ],
                                                                                                                  ε in [
                                                                                                                      0.0,
                                                                                                                      0.1,
                                                                                                                      0.5,
                                                                                                                  ],
                                                                                                                  ncheb in [
                                                                                                                      2,
                                                                                                                      3,
                                                                                                                      4,
                                                                                                                      10,
                                                                                                                  ]]
            @test begin
                chebyshev = TMatrix.Chebyshev(; m = m, ε = ε, rea = rea, n = ncheb, λ = λ)
                rev = TMatrix.volume_equivalent_radius(chebyshev)
                ratio = TMatrix.Wrapper.Fixed.surfch(ncheb, ε)
                rev ≈ rea * ratio
            end
        end
    end

    @testset "Calculate constants" begin
        m = 1
        λ = 1
        rev = 1

        @testset "for spheroids with a_to_c = $a_to_c and nmax = $nmax" for (a_to_c, nmax) in [(a_to_c,
                                                                                                nmax)
                                                                                               for a_to_c in [
                                                                                                       0.5,
                                                                                                       1.0,
                                                                                                       2.0,
                                                                                                   ],
                                                                                                   nmax in [
                                                                                                       4,
                                                                                                       10,
                                                                                                       20,
                                                                                                   ]]
            @test begin
                spheroid = TMatrix.Spheroid(; m = m, a_to_c = a_to_c, rev = rev, λ = λ)
                ngauss = nmax * 4
                x0, w0, an0, ann0, s0, ss0 = TMatrix.Wrapper.Fixed.constant(ngauss, nmax,
                                                                            -1, a_to_c)
                x0 = x0[1:ngauss]
                w0 = w0[1:ngauss]
                an0 = an0[1:nmax]
                ann0 = ann0[1:nmax, 1:nmax]
                s0 = s0[1:ngauss]
                ss0 = ss0[1:ngauss]
                x, w, an, ann, s, ss = TMatrix.constant(spheroid, ngauss, nmax)
                all(x .≈ x0) && all(w .≈ w0) && all(an .≈ an0) && all(ann .≈ ann0) &&
                    all(s .≈ s0) && all(ss .≈ ss0)
            end
        end

        @testset "for cylinders with r_to_h = $r_to_h and nmax = $nmax" for (r_to_h, nmax) in [(r_to_h,
                                                                                                nmax)
                                                                                               for r_to_h in [
                                                                                                       0.25,
                                                                                                       0.5,
                                                                                                       1.0,
                                                                                                   ],
                                                                                                   nmax in [
                                                                                                       4,
                                                                                                       10,
                                                                                                       20,
                                                                                                   ]]
            @test begin
                cylinder = TMatrix.Cylinder(; m = m, r_to_h = r_to_h, rev = rev, λ = λ)
                ngauss = nmax * 4
                x0, w0, an0, ann0, s0, ss0 = TMatrix.Wrapper.Fixed.constant(ngauss, nmax,
                                                                            -2, 2r_to_h)
                x0 = x0[1:ngauss]
                w0 = w0[1:ngauss]
                an0 = an0[1:nmax]
                ann0 = ann0[1:nmax, 1:nmax]
                s0 = s0[1:ngauss]
                ss0 = ss0[1:ngauss]
                x, w, an, ann, s, ss = TMatrix.constant(cylinder, ngauss, nmax)
                all(x .≈ x0) && all(w .≈ w0) && all(an .≈ an0) && all(ann .≈ ann0) &&
                    all(s .≈ s0) && all(ss .≈ ss0)
            end
        end

        @testset "for Chebyshev particles with ε = $ε, ncheb = $ncheb and nmax = $nmax" for (ε, ncheb, nmax) in [(ε,
                                                                                                                  ncheb,
                                                                                                                  nmax)
                                                                                                                 for ε in [
                                                                                                                         0.0,
                                                                                                                         0.1,
                                                                                                                         0.5,
                                                                                                                     ],
                                                                                                                     ncheb in [
                                                                                                                         2,
                                                                                                                         3,
                                                                                                                         4,
                                                                                                                         10,
                                                                                                                     ],
                                                                                                                     nmax in [
                                                                                                                         4,
                                                                                                                         10,
                                                                                                                         20,
                                                                                                                     ]]
            @test begin
                chebyshev = TMatrix.Chebyshev(; m = m, n = ncheb, ε = ε, rev = rev, λ = λ)
                ngauss = nmax * 4
                x0, w0, an0, ann0, s0, ss0 = TMatrix.Wrapper.Fixed.constant(ngauss, nmax,
                                                                            ncheb, ε)
                x0 = x0[1:ngauss]
                w0 = w0[1:ngauss]
                an0 = an0[1:nmax]
                ann0 = ann0[1:nmax, 1:nmax]
                s0 = s0[1:ngauss]
                ss0 = ss0[1:ngauss]
                x, w, an, ann, s, ss = TMatrix.constant(chebyshev, ngauss, nmax)
                all(x .≈ x0) && all(w .≈ w0) && all(an .≈ an0) && all(ann .≈ ann0) &&
                    all(s .≈ s0) && all(ss .≈ ss0)
            end
        end
    end

    @testset "Special functions" begin
        nmax = 50
        nnmax = nmax * 4
        @testset "Spherical Bessel j for real numbers (rjb)" for x in [
            -10.0,
            -1.5,
            -1.0,
            -0.5,
            0.01,
            0.5,
            1.0,
            1.5,
            10.0,
        ]
            @test begin
                y0, u0 = TMatrix.Wrapper.Fixed.rjb(x, nmax, nnmax)
                y, u = TMatrix.sphericalbesselj(x, nmax, nnmax)
                all(y .≈ y0) && all(u .≈ u0)
            end
        end

        @testset "Spherical Bessel j for complex numbers (cjb)" for x in [
            -10.0,
            -1.5 - 0.6im,
            -1.0 + 0.2im,
            -0.5 + 0.0im,
            0.0 - 0.1im,
            0.5 - 0.3im,
            1.0 + 2.1im,
            1.5 + 0.02im,
            10.0,
        ]
            @test begin
                y0, u0 = TMatrix.Wrapper.Fixed.cjb(x, nmax, nnmax)
                y, u = TMatrix.sphericalbesselj(x, nmax, nnmax)
                all(y .≈ y0) && all(u .≈ u0)
            end
        end

        @testset "Spherical Bessel y for real numbers (ryb)" for x in [
            -10.0,
            -1.5,
            -1.0,
            -0.5,
            0.01,
            0.5,
            1.0,
            1.5,
            10.0,
        ]
            @test begin
                y0, u0 = TMatrix.Wrapper.Fixed.ryb(x, nmax)
                y, u = TMatrix.sphericalbessely(x, nmax)
                all(y .≈ y0) && all(u .≈ u0)
            end
        end
    end

    @testset "Calculate Bessel variables" begin
        λ = 1.0
        m = 1.5 + 0.02im
        rev = 1.0
        nmax = 50
        ngauss = nmax * 4

        @testset "for spheroids with a_to_c = $a_to_c" for a_to_c in [0.5, 1.0, 2.0]
            @test begin
                spheroid = TMatrix.Spheroid(; m = m, a_to_c = a_to_c, rev = rev, λ = λ)

                np = -1
                x, _ = TMatrix.Wrapper.Fixed.constant(ngauss, nmax, np, a_to_c)
                _, _, _, _, _, ddr, drr, dri, jkr0, djkr0, ykr0, dykr0, jkₛr0, djkₛr0 = TMatrix.Wrapper.Fixed.vary(x,
                                                                                                                   λ,
                                                                                                                   m,
                                                                                                                   rev,
                                                                                                                   a_to_c,
                                                                                                                   np,
                                                                                                                   ngauss,
                                                                                                                   nmax)
                _, _, kr1, kₛr1, jkr, djkr, ykr, dykr, jkₛr, djkₛr = TMatrix.vary(spheroid,
                                                                                  ngauss,
                                                                                  nmax)

                all(ddr[1:ngauss] .≈ kr1) &&
                    all((drr + 1.0im * dri)[1:ngauss] .≈ kₛr1) &&
                    all(jkr0 .≈ jkr) &&
                    all(djkr0 .≈ djkr) &&
                    all(ykr0 .≈ ykr) &&
                    all(dykr0 .≈ dykr) &&
                    all(jkₛr0 .≈ jkₛr) &&
                    all(djkₛr0 .≈ djkₛr)
            end
        end

        @testset "for cylinders with r_to_h = $r_to_h" for r_to_h in [0.25, 0.5, 1.0]
            @test begin
                cylinder = TMatrix.Cylinder(; m = m, r_to_h = r_to_h, rev = rev, λ = λ)

                np = -2
                x, _ = TMatrix.Wrapper.Fixed.constant(ngauss, nmax, np, 2r_to_h)
                _, _, _, _, _, ddr, drr, dri, jkr0, djkr0, ykr0, dykr0, jkₛr0, djkₛr0 = TMatrix.Wrapper.Fixed.vary(x,
                                                                                                                   λ,
                                                                                                                   m,
                                                                                                                   rev,
                                                                                                                   2r_to_h,
                                                                                                                   np,
                                                                                                                   ngauss,
                                                                                                                   nmax)
                _, _, kr1, kₛr1, jkr, djkr, ykr, dykr, jkₛr, djkₛr = TMatrix.vary(cylinder,
                                                                                  ngauss,
                                                                                  nmax)

                all(ddr[1:ngauss] .≈ kr1) &&
                    all((drr + 1.0im * dri)[1:ngauss] .≈ kₛr1) &&
                    all(jkr0 .≈ jkr) &&
                    all(djkr0 .≈ djkr) &&
                    all(ykr0 .≈ ykr) &&
                    all(dykr0 .≈ dykr) &&
                    all(jkₛr0 .≈ jkₛr) &&
                    all(djkₛr0 .≈ djkₛr)
            end
        end

        @testset "for Chebyshev particles with ε = $ε and ncheb = $ncheb" for (ε, ncheb) in [(ε,
                                                                                              ncheb)
                                                                                             for ε in [
                                                                                                     0.0,
                                                                                                     0.1,
                                                                                                     0.5,
                                                                                                 ],
                                                                                                 ncheb in [
                                                                                                     2,
                                                                                                     3,
                                                                                                     4,
                                                                                                     10,
                                                                                                 ]]
            @test begin
                chebyshev = TMatrix.Chebyshev(; m = m, ε = ε, n = ncheb, rev = rev, λ = λ)

                np = ncheb
                x, _ = TMatrix.Wrapper.Fixed.constant(ngauss, nmax, np, ε)
                _, _, _, _, _, ddr, drr, dri, jkr0, djkr0, ykr0, dykr0, jkₛr0, djkₛr0 = TMatrix.Wrapper.Fixed.vary(x,
                                                                                                                   λ,
                                                                                                                   m,
                                                                                                                   rev,
                                                                                                                   ε,
                                                                                                                   np,
                                                                                                                   ngauss,
                                                                                                                   nmax)
                _, _, kr1, kₛr1, jkr, djkr, ykr, dykr, jkₛr, djkₛr = TMatrix.vary(chebyshev,
                                                                                  ngauss,
                                                                                  nmax)

                all(ddr[1:ngauss] .≈ kr1) &&
                    all((drr + 1.0im * dri)[1:ngauss] .≈ kₛr1) &&
                    all(jkr0 .≈ jkr) &&
                    all(djkr0 .≈ djkr) &&
                    all(ykr0 .≈ ykr) &&
                    all(dykr0 .≈ dykr) &&
                    all(jkₛr0 .≈ jkₛr) &&
                    all(djkₛr0 .≈ djkₛr)
            end
        end
    end

    @testset "Calculate sub T-Matrix" begin
        λ = 1.0
        m = 1.5 + 0.02im
        rev = 1.0
        nmax = 20
        ngauss = nmax * 4

        @testset "for spheroids with a_to_c = $a_to_c" for a_to_c in [0.5, 1.01, 2.0]
            @test begin
                spheroid = TMatrix.Spheroid(; m = m, a_to_c = a_to_c, rev = rev, λ = λ)

                np = -1
                valid = true
                for mm in 0:nmax
                    if mm == 0
                        T0, _ = TMatrix.Wrapper.Fixed.tmatr0(ngauss, nmax, np, a_to_c, λ, m,
                                                             rev)
                        T, _ = TMatrix.tmatr0!(spheroid, ngauss, nmax)
                    else
                        T0, _ = TMatrix.Wrapper.Fixed.tmatr(mm, ngauss, nmax, np, a_to_c, λ,
                                                            m, rev)
                        T, _ = TMatrix.tmatr!(spheroid, mm, ngauss, nmax)
                    end
                    if !all(isapprox.(T0, T, rtol = RTOL, atol = ATOL))
                        valid = false
                        break
                    end
                end

                valid
            end
        end

        @testset "for cylinders with r_to_h = $r_to_h" for r_to_h in [0.25, 0.5, 1.0]
            @test begin
                cylinder = TMatrix.Cylinder(; m = m, r_to_h = r_to_h, rev = rev, λ = λ)

                np = -2
                valid = true
                for mm in 0:nmax
                    if mm == 0
                        T0, _ = TMatrix.Wrapper.Fixed.tmatr0(ngauss, nmax, np, 2r_to_h, λ,
                                                             m, rev)
                        T, _ = TMatrix.tmatr0!(cylinder, ngauss, nmax)
                    else
                        T0, _ = TMatrix.Wrapper.Fixed.tmatr(mm, ngauss, nmax, np, 2r_to_h,
                                                            λ, m, rev)
                        T, _ = TMatrix.tmatr!(cylinder, mm, ngauss, nmax)
                    end
                    if !all(isapprox.(T0, T, rtol = RTOL, atol = ATOL))
                        valid = false
                        break
                    end
                end

                valid
            end
        end

        @testset "for Chebyshev particles with ε = $ε and ncheb = $ncheb" for (ε, ncheb) in [(ε,
                                                                                              ncheb)
                                                                                             for ε in [
                                                                                                     -0.15,
                                                                                                     0.01,
                                                                                                     0.1,
                                                                                                 ],
                                                                                                 ncheb in [
                                                                                                     2,
                                                                                                     3,
                                                                                                     4,
                                                                                                     10,
                                                                                                 ]]
            @test begin
                chebyshev = TMatrix.Chebyshev(; m = m, ε = ε, n = ncheb, rev = rev, λ = λ)

                np = ncheb
                valid = true
                for mm in 0:nmax
                    if mm == 0
                        T0, _ = TMatrix.Wrapper.Fixed.tmatr0(ngauss, nmax, np, ε, λ, m, rev)
                        T, _ = TMatrix.tmatr0!(chebyshev, ngauss, nmax)
                    else
                        T0, _ = TMatrix.Wrapper.Fixed.tmatr(mm, ngauss, nmax, np, ε, λ, m,
                                                            rev)
                        T, _ = TMatrix.tmatr!(chebyshev, mm, ngauss, nmax)
                    end
                    if !all(isapprox.(T0, T, rtol = RTOL, atol = ATOL))
                        valid = false
                        break
                    end
                end

                valid
            end
        end
    end

    @testset "Calculate amplitude" begin
        λ = 1.0
        m = 1.5 + 0.02im
        rev = 1.0
        ratio = 1.0
        ddelta = 0.001
        ndgs = 4
        ϑ_i = 56.0
        ϑ_s = 65.0
        φ_i = 114.0
        φ_s = 128.0
        α = 145.0
        β = 52.0

        @testset "for spheroids with a_to_c = $a_to_c" for a_to_c in [0.5, 1.01, 2.0]
            @test begin
                spheroid = TMatrix.Spheroid(; m = m, a_to_c = a_to_c, rev = rev, λ = λ)

                T = TMatrix.calc_tmatrix!(spheroid, ddelta, ndgs)
                S, Z = TMatrix.calc_SZ(spheroid, α, β, ϑ_i, ϑ_s, φ_i, φ_s, T)

                np = -1
                T0, nmax = TMatrix.Wrapper.Fixed.calc_tmatrix(rev, ratio, λ, m, a_to_c, np,
                                                              ddelta, ndgs ÷ 2)
                S0, Z0 = TMatrix.Wrapper.Fixed.calc_SZ(nmax, λ, α, β, ϑ_i, ϑ_s, φ_i, φ_s)

                all(isapprox.(T, T0, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(S, S0, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(Z, Z0, rtol = RTOL, atol = ATOL))
            end
        end

        @testset "for cylinders with r_to_h = $r_to_h" for r_to_h in [0.5, 1.0, 2.0]
            @test begin
                cylinder = TMatrix.Cylinder(; m = m, r_to_h = r_to_h, rev = rev, λ = λ)

                T = TMatrix.calc_tmatrix!(cylinder, ddelta, ndgs)
                S, Z = TMatrix.calc_SZ(cylinder, α, β, ϑ_i, ϑ_s, φ_i, φ_s, T)

                np = -2
                T0, nmax = TMatrix.Wrapper.Fixed.calc_tmatrix(rev, ratio, λ, m, 2r_to_h, np,
                                                              ddelta, ndgs ÷ 2)
                S0, Z0 = TMatrix.Wrapper.Fixed.calc_SZ(nmax, λ, α, β, ϑ_i, ϑ_s, φ_i, φ_s)

                all(isapprox.(T, T0, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(S, S0, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(Z, Z0, rtol = RTOL, atol = ATOL))
            end
        end

        @testset "for Chebyshev particles with ε = $ε and ncheb = $ncheb" for (ε, ncheb) in [(ε,
                                                                                              ncheb)
                                                                                             for ε in [
                                                                                                     -0.15,
                                                                                                     0.01,
                                                                                                     0.1,
                                                                                                 ],
                                                                                                 ncheb in [
                                                                                                     2,
                                                                                                     3,
                                                                                                     4,
                                                                                                 ]]
            @test begin
                chebyshev = TMatrix.Chebyshev(; m = m, ε = ε, n = ncheb, rev = rev, λ = λ)
                T = TMatrix.calc_tmatrix!(chebyshev, ddelta, ndgs)
                S, Z = TMatrix.calc_SZ(chebyshev, α, β, ϑ_i, ϑ_s, φ_i, φ_s, T)

                np = ncheb
                T0, nmax = TMatrix.Wrapper.Fixed.calc_tmatrix(rev, ratio, λ, m, ε, np,
                                                              ddelta, ndgs ÷ 2)
                S0, Z0 = TMatrix.Wrapper.Fixed.calc_SZ(nmax, λ, α, β, ϑ_i, ϑ_s, φ_i, φ_s)

                all(isapprox.(T, T0, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(S, S0, rtol = RTOL, atol = ATOL)) &&
                    all(isapprox.(Z, Z0, rtol = RTOL, atol = ATOL))
            end
        end
    end
end
