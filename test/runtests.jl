using TMatrix
using TMatrix.Wrapper
using Libdl
using Test

@testset "TMatrix.jl" begin
    tm = Libdl.dlopen("../shared/tmatrix.so")

    @testset "Calculate vig function" begin
        @testset "vig($nmax, $m, $x)" for (nmax, m, x) in [
            (nmax, m, x) for nmax in [5, 10, 100], m in [0, 1, 2, 3], x in [-0.9, -0.5, -0.1, 0.1, 0.5, 0.9]
        ]
            @test begin
                dv10, dv20 = TMatrix.Wrapper.vig(tm, nmax, m, x)
                dv1, dv2 = TMatrix.vig(nmax, m, x)
                dv1 ≈ dv10
                dv2 ≈ dv20
            end
        end
    end

    @testset "Calculate vigampl function" begin
        @testset "vigampl($nmax, $m, $x)" for (nmax, m, x) in [
            (nmax, m, x) for nmax in [5, 10, 100], m in [0, 1, 2, 3], x in [-1.0, -0.5, -0.1, 0.1, 0.5, 1.0]
        ]
            @test begin
                dv10, dv20 = TMatrix.Wrapper.vigampl(tm, nmax, m, x)
                dv1, dv2 = TMatrix.vigampl(nmax, m, x)
                dv1 ≈ dv10
                dv2 ≈ dv20
            end
        end
    end

    @testset "Calculate r(θ) and dr/dθ" begin
        @testset "for spheroids with rev = $rev, m = $m, a_to_c = $a_to_c and ngauss = $ngauss" for (
            rev,
            m,
            a_to_c,
            ngauss,
        ) in [
            (rev, m, a_to_c, ngauss) for rev in [0.5, 1.0, 2.0], m in [1.0, 0.9 + 0.001im, 1.1 - 0.001im],
            a_to_c in [0.5, 1.0, 2.0], ngauss in [4, 20, 1000]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = rev,
                    shape = TMatrix.SHAPE_SPHEROID,
                    refractive_index = m,
                    axis_ratio = a_to_c,
                )
                r2, drr = TMatrix.Wrapper.rsp1(tm, ngauss, rev, a_to_c)
                r, dr = TMatrix.calc_r(scatterer, ngauss)

                r .^ 2 ≈ r2 && dr ./ r ≈ drr
            end
        end

        @testset "for cylinders with rev = $rev, m = $m, d_to_h = $d_to_h and ngauss = $ngauss" for (
            rev,
            m,
            d_to_h,
            ngauss,
        ) in [
            (rev, m, d_to_h, ngauss) for rev in [0.5, 1.0, 2.0], m in [1.0, 0.9 + 0.001im, 1.1 - 0.001im],
            d_to_h in [0.5, 1.0, 2.0], ngauss in [4, 20, 1000]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = rev,
                    shape = TMatrix.SHAPE_CYLINDER,
                    refractive_index = m,
                    axis_ratio = d_to_h,
                )
                r2, drr = TMatrix.Wrapper.rsp3(tm, ngauss, rev, d_to_h)
                r, dr = TMatrix.calc_r(scatterer, ngauss)

                r .^ 2 ≈ r2 && dr ./ r ≈ drr
            end
        end

        @testset "for Chebyshev particles with rev = $rev, m = $m, ε = $ε, ngauss = $ngauss and ncheb = $ncheb" for (
            rev,
            m,
            ε,
            ngauss,
            ncheb,
        ) in [
            (rev, m, ε, ngauss, ncheb) for rev in [0.5, 1.0, 2.0], m in [1.0, 0.9 + 0.001im, 1.1 - 0.001im],
            ε in [0.0, 0.1, 0.5, 0.9], ngauss in [4, 20, 1000], ncheb in [2, 3, 4, 10]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = rev,
                    shape = TMatrix.SHAPE_CHEBYSHEV,
                    refractive_index = m,
                    axis_ratio = ε,
                    n = ncheb,
                )
                r2, drr = TMatrix.Wrapper.rsp2(tm, ngauss, rev, ε, ncheb)
                r, dr = TMatrix.calc_r(scatterer, ngauss)

                r .^ 2 ≈ r2 && dr ./ r ≈ drr
            end
        end
    end

    @testset "Calculate rev from rea" begin
        @testset "for spheroids with rea = $rea and a_to_c = $a_to_c" for (rea, a_to_c) in [
            (rea, a_to_c) for rea in [0.5, 1.0, 2.0], a_to_c in [0.5, 1.000001, 2.0]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = rea,
                    shape = TMatrix.SHAPE_SPHEROID,
                    radius_type = TMatrix.RADIUS_EQUAL_AREA,
                    refractive_index = 1.0 + 0.0im,
                    axis_ratio = a_to_c,
                )
                rev = scatterer.rev
                ratio = TMatrix.Wrapper.sarea(tm, a_to_c)
                rev ≈ rea * ratio
            end
        end

        @testset "for cylinders with rea = $rea and d_to_h = $d_to_h" for (rea, d_to_h) in [
            (rea, d_to_h) for rea in [0.5, 1.0, 2.0], d_to_h in [0.5, 1.0, 2.0]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = rea,
                    shape = TMatrix.SHAPE_CYLINDER,
                    radius_type = TMatrix.RADIUS_EQUAL_AREA,
                    refractive_index = 1.0 + 0.0im,
                    axis_ratio = d_to_h,
                )
                rev = scatterer.rev
                ratio = TMatrix.Wrapper.sareac(tm, d_to_h)
                rev ≈ rea * ratio
            end
        end

        @testset "for Chebyshev particles with rea = $rea, ε = $ε and ncheb = $ncheb" for (rea, ε, ncheb) in [
            (rea, ε, ncheb) for rea in [0.5, 1.0, 2.0], ε in [0.0, 0.1, 0.5, 0.9], ncheb in [2, 3, 4, 10]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = rea,
                    shape = TMatrix.SHAPE_CHEBYSHEV,
                    radius_type = TMatrix.RADIUS_EQUAL_AREA,
                    refractive_index = 1.0 + 0.0im,
                    axis_ratio = ε,
                    n = ncheb,
                )
                rev = scatterer.rev
                ratio = TMatrix.Wrapper.surfch(tm, ncheb, ε)
                rev ≈ rea * ratio
            end
        end
    end

    @testset "Calculate constants" begin
        @testset "for spheroids with a_to_c = $a_to_c and nmax = $nmax" for (a_to_c, nmax) in [
            (a_to_c, nmax) for a_to_c in [0.5, 1.0, 2.0], nmax in [4, 10, 20]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = 1.0,
                    shape = TMatrix.SHAPE_SPHEROID,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = 1.0 + 0.0im,
                    axis_ratio = a_to_c,
                )

                ngauss = nmax * 2
                x0, w0, an0, ann0, s0, ss0 = TMatrix.Wrapper.const_(tm, ngauss, nmax, -1, a_to_c)
                x0 = x0[1:ngauss]
                w0 = w0[1:ngauss]
                an0 = an0[1:nmax]
                ann0 = ann0[1:nmax, 1:nmax]
                s0 = s0[1:ngauss]
                ss0 = ss0[1:ngauss]
                x, w, an, ann, s, ss = TMatrix.const_(scatterer, ngauss, nmax)
                x ≈ x0 && w ≈ w0 && an ≈ an0 && ann ≈ ann0 && s ≈ s0 && ss ≈ ss0
            end
        end

        @testset "for cylinders with d_to_h = $d_to_h and nmax = $nmax" for (d_to_h, nmax) in [
            (d_to_h, nmax) for d_to_h in [0.5, 1.0, 2.0], nmax in [4, 10, 20]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = 1.0,
                    shape = TMatrix.SHAPE_CYLINDER,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = 1.0 + 0.0im,
                    axis_ratio = d_to_h,
                )

                ngauss = nmax * 2
                x0, w0, an0, ann0, s0, ss0 = TMatrix.Wrapper.const_(tm, ngauss, nmax, -2, d_to_h)
                x0 = x0[1:ngauss]
                w0 = w0[1:ngauss]
                an0 = an0[1:nmax]
                ann0 = ann0[1:nmax, 1:nmax]
                s0 = s0[1:ngauss]
                ss0 = ss0[1:ngauss]
                x, w, an, ann, s, ss = TMatrix.const_(scatterer, ngauss, nmax)
                x ≈ x0 && w ≈ w0 && an ≈ an0 && ann ≈ ann0 && s ≈ s0 && ss ≈ ss0
            end
        end

        @testset "for Chebyshev particles with ε = $ε, ncheb = $ncheb and nmax = $nmax" for (ε, ncheb, nmax) in [
            (ε, ncheb, nmax) for ε in [0.0, 0.1, 0.5, 0.9], ncheb in [2, 3, 4, 10], nmax in [4, 10, 20]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = 1.0,
                    shape = TMatrix.SHAPE_CHEBYSHEV,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = 1.0 + 0.0im,
                    axis_ratio = ε,
                    n = ncheb,
                )

                ngauss = nmax * 2
                x0, w0, an0, ann0, s0, ss0 = TMatrix.Wrapper.const_(tm, ngauss, nmax, ncheb, ε)
                x0 = x0[1:ngauss]
                w0 = w0[1:ngauss]
                an0 = an0[1:nmax]
                ann0 = ann0[1:nmax, 1:nmax]
                s0 = s0[1:ngauss]
                ss0 = ss0[1:ngauss]
                x, w, an, ann, s, ss = TMatrix.const_(scatterer, ngauss, nmax)
                x ≈ x0 && w ≈ w0 && an ≈ an0 && ann ≈ ann0 && s ≈ s0 && ss ≈ ss0
            end
        end
    end

    @testset "Special functions" begin
        nmax = 50
        nnmax = nmax * 4
        @testset "Spherical Bessel j for real numbers (rjb)" for x in
                                                                 [-10.0, -1.5, -1.0, -0.5, 0.01, 0.5, 1.0, 1.5, 10.0]
            @test begin
                y0, u0 = TMatrix.Wrapper.rjb(tm, x, nmax, nnmax)
                y, u = TMatrix.sphericalbesselj(x, nmax, nnmax)
                y ≈ y0 && u ≈ u0
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
                y0, u0 = TMatrix.Wrapper.cjb(tm, x, nmax, nnmax)
                y, u = TMatrix.sphericalbesselj(x, nmax, nnmax)
                y ≈ y0 && u ≈ u0
            end
        end

        @testset "Spherical Bessel y for real numbers (ryb)" for x in
                                                                 [-10.0, -1.5, -1.0, -0.5, 0.01, 0.5, 1.0, 1.5, 10.0]
            @test begin
                y0, u0 = TMatrix.Wrapper.ryb(tm, x, nmax)
                y, u = TMatrix.sphericalbessely(x, nmax)
                y ≈ y0 && u ≈ u0
            end
        end
    end

    @testset "Calculate Bessel variables" begin
        λ = 1.0
        m = 1.5 + 0.02im
        rev = 1.0
        nmax = 50
        ngauss = nmax * 2

        @testset "for spheroids with a_to_c = $a_to_c" for a_to_c in [0.5, 1.0, 2.0]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = rev,
                    shape = TMatrix.SHAPE_SPHEROID,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = m,
                    axis_ratio = a_to_c,
                )

                np = -1
                x, _ = TMatrix.Wrapper.const_(tm, ngauss, nmax, np, a_to_c)
                _, _, _, _, _, _, _, _, jkr0, djkr0, ykr0, dykr0, jkr_s0, djkr_s0 =
                    TMatrix.Wrapper.vary(tm, x, λ, m, rev, a_to_c, np, ngauss, nmax)
                _, _, jkr, djkr, ykr, dykr, jkr_s, djkr_s = TMatrix.vary(scatterer, ngauss, nmax)

                jkr0 ≈ jkr && djkr0 ≈ djkr && ykr0 ≈ ykr && dykr0 ≈ dykr && jkr_s0 ≈ jkr_s && djkr_s0 ≈ djkr_s
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
                )

                np = -2
                x, _ = TMatrix.Wrapper.const_(tm, ngauss, nmax, np, d_to_h)
                _, _, _, _, _, _, _, _, jkr0, djkr0, ykr0, dykr0, jkr_s0, djkr_s0 =
                    TMatrix.Wrapper.vary(tm, x, λ, m, rev, d_to_h, np, ngauss, nmax)
                _, _, jkr, djkr, ykr, dykr, jkr_s, djkr_s = TMatrix.vary(scatterer, ngauss, nmax)

                jkr0 ≈ jkr && djkr0 ≈ djkr && ykr0 ≈ ykr && dykr0 ≈ dykr && jkr_s0 ≈ jkr_s && djkr_s0 ≈ djkr_s
            end
        end

        @testset "for Chebyshev particles with ε = $ε and ncheb = $ncheb" for (ε, ncheb) in [
            (ε, ncheb) for ε in [0.0, 0.1, 0.5, 0.9], ncheb in [2, 3, 4, 10]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r = rev,
                    shape = TMatrix.SHAPE_CHEBYSHEV,
                    radius_type = TMatrix.RADIUS_EQUAL_VOLUME,
                    refractive_index = m,
                    axis_ratio = ε,
                    n = ncheb,
                )

                np = ncheb
                x, _ = TMatrix.Wrapper.const_(tm, ngauss, nmax, np, ε)
                _, _, _, _, _, _, _, _, jkr0, djkr0, ykr0, dykr0, jkr_s0, djkr_s0 =
                    TMatrix.Wrapper.vary(tm, x, λ, m, rev, ε, np, ngauss, nmax)
                _, _, jkr, djkr, ykr, dykr, jkr_s, djkr_s = TMatrix.vary(scatterer, ngauss, nmax)

                jkr0 ≈ jkr && djkr0 ≈ djkr && ykr0 ≈ ykr && dykr0 ≈ dykr && jkr_s0 ≈ jkr_s && djkr_s0 ≈ djkr_s
            end
        end
    end
end
