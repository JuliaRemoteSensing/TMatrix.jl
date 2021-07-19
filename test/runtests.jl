using TMatrix
using TMatrix.Wrapper
using Libdl
using Test

@testset "TMatrix.jl" begin
    tm = Libdl.dlopen("../shared/tmatrix.so")
    @testset "Calculate vig function" begin
        @testset "vig($nmax, $m, $x)" for (nmax, m, x) in [(nmax, m, x) for nmax in [5, 10, 100] for m in [0, 1, 2, 3] for x in [-0.9, -0.5, -0.1, 0.1, 0.5, 0.9]]
            @test begin
                dv10, dv20 = TMatrix.Wrapper.vig(tm, nmax, m, x)
                dv1, dv2 = TMatrix.vig(nmax, m, x)
                dv1 ≈ dv10
                dv2 ≈ dv20
            end
        end
    end

    @testset "Calculate vigampl function" begin
        @testset "vigampl($nmax, $m, $x)" for (nmax, m, x) in [(nmax, m, x) for nmax in [5, 10, 100] for m in [0, 1, 2, 3] for x in [-1.0, -0.5, -0.1, 0.1, 0.5, 1.0]]
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
            (rev, m, a_to_c, ngauss) for rev in [0.5, 1.0, 2.0] for m in [1.0, 0.9 + 0.001im, 1.1 - 0.001im] for a_to_c in [0.5, 1.0, 2.0] for ngauss in [4, 20, 1000]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r=rev,
                    shape=TMatrix.SHAPE_SPHEROID,
                    refractive_index=m,
                    axis_ratio=a_to_c,
                )
                r2, drr = TMatrix.Wrapper.rsp1(tm, ngauss, rev, a_to_c)
                r, dr = TMatrix.calc_r(scatterer, ngauss)

                r.^2 ≈ r2 && dr ./ r ≈ drr
            end
        end

        @testset "for cylinders with rev = $rev, m = $m, d_to_h = $d_to_h and ngauss = $ngauss" for (
            rev,
            m,
            d_to_h,
            ngauss,
        ) in [
            (rev, m, d_to_h, ngauss) for rev in [0.5, 1.0, 2.0] for m in [1.0, 0.9 + 0.001im, 1.1 - 0.001im] for d_to_h in [0.5, 1.0, 2.0] for ngauss in [4, 20, 1000]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r=rev,
                    shape=TMatrix.SHAPE_CYLINDER,
                    refractive_index=m,
                    axis_ratio=d_to_h,
                )
                r2, drr = TMatrix.Wrapper.rsp3(tm, ngauss, rev, d_to_h)
                r, dr = TMatrix.calc_r(scatterer, ngauss)

                r.^2 ≈ r2 && dr ./ r ≈ drr
            end
        end

        @testset "for Chebyshev particles with rev = $rev, m = $m, ε = $ε, ngauss = $ngauss and ncheb = $ncheb" for (
            rev,
            m,
            ε,
            ngauss,
            ncheb,
        ) in [
            (rev, m, ε, ngauss, ncheb) for rev in [0.5, 1.0, 2.0] for m in [1.0, 0.9 + 0.001im, 1.1 - 0.001im] for ε in [0.5, 1.0, 2.0] for ngauss in [4, 20, 1000] for ncheb in [2, 3, 4, 10]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r=rev,
                    shape=TMatrix.SHAPE_CHEBYSHEV,
                    refractive_index=m,
                    axis_ratio=ε,
                    n=ncheb,
                )
                r2, drr = TMatrix.Wrapper.rsp2(tm, ngauss, rev, ε, ncheb)
                r, dr = TMatrix.calc_r(scatterer, ngauss)

                r.^2 ≈ r2 && dr ./ r ≈ drr
            end
        end
    end

    @testset "Calculate rev from rea" begin
                @testset "for spheroids with rea = $rea and a_to_c = $a_to_c" for (rea, a_to_c) in [
            (rea, a_to_c) for rea in [0.5, 1.0, 2.0] for a_to_c in [0.5, 1.000001, 2.0]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r=rea,
                    shape=TMatrix.SHAPE_SPHEROID,
                    radius_type=TMatrix.RADIUS_EQUAL_AREA,
                    refractive_index=1.0 + 0.0im,
                    axis_ratio=a_to_c,
                )
                rev = scatterer.rev
                ratio = TMatrix.Wrapper.sarea(tm, a_to_c)
                rev ≈ rea * ratio
            end
        end

        @testset "for cylinders with rea = $rea and d_to_h = $d_to_h" for (rea, d_to_h) in [
            (rea, d_to_h) for rea in [0.5, 1.0, 2.0] for d_to_h in [0.5, 1.0, 2.0]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r=rea,
                    shape=TMatrix.SHAPE_CYLINDER,
                    radius_type=TMatrix.RADIUS_EQUAL_AREA,
                    refractive_index=1.0 + 0.0im,
                    axis_ratio=d_to_h,
                )
                rev = scatterer.rev
                ratio = TMatrix.Wrapper.sareac(tm, d_to_h)
                rev ≈ rea * ratio
            end
        end

        @testset "for Chebyshev particles with rea = $rea, ε = $ε and ncheb = $ncheb" for (rea, ε, ncheb) in [
            (rea, ε, ncheb) for rea in [0.5, 1.0, 2.0] for ε in [0.5, 1.0, 2.0] for ncheb in [2, 3, 4, 10]
        ]
            @test begin
                scatterer = TMatrix.Scatterer(
                    r=rea,
        shape=TMatrix.SHAPE_CHEBYSHEV,
        radius_type=TMatrix.RADIUS_EQUAL_AREA,
                    refractive_index=1.0 + 0.0im,
                    axis_ratio=ε,
                    n=ncheb,
                )
                rev = scatterer.rev
                ratio = TMatrix.Wrapper.surfch(tm, ncheb, ε)
                rev ≈ rea * ratio
            end
        end
    end
end
