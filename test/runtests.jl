using TMatrix
using TMatrix.Wrapper
using Libdl
using Test

@testset "TMatrix.jl" begin
    tm = Libdl.dlopen("../shared/tmatrix.so")

    @testset "Calculate r(θ) and dr/dθ" begin
        @testset "for spheroids" begin
            for rev in [0.5, 1.0, 2.0] 
                for m in [1.0, 0.9 + 0.001im, 1.1 - 0.001im]
                    for a_to_c in [0.5, 1.0, 2.0]
                        for ngauss in [4, 20, 1000]
                            @test begin
                                scatterer = TMatrix.Scatterer(r=rev, shape=TMatrix.SHAPE_SPHEROID, refractive_index=m, axis_ratio=a_to_c)
                                r2, drr = TMatrix.Wrapper.rsp1(tm, ngauss, rev, a_to_c)
                                r, dr = TMatrix.calc_r(scatterer, ngauss)

                                r.^2 ≈ r2 && dr ./ r ≈ drr
                            end
                        end
                    end
                end
            end
        end

        @testset "for cylinder" begin
            for rev in [0.5, 1.0, 2.0] 
                for m in [1.0, 0.9 + 0.001im, 1.1 - 0.001im]
                    for d_to_h in [0.5, 1.0, 2.0]
                        for ngauss in [4, 20, 1000]
                            @test begin
                                scatterer = TMatrix.Scatterer(r=rev, shape=TMatrix.SHAPE_CYLINDER, refractive_index=m, axis_ratio=d_to_h)
                                r2, drr = TMatrix.Wrapper.rsp3(tm, ngauss, rev, d_to_h)
                                r, dr = TMatrix.calc_r(scatterer, ngauss)

                                r.^2 ≈ r2 && dr ./ r ≈ drr
                            end
                        end
                    end
                end
            end
        end
        
        @testset "for Chebyshev particles" begin
            for rev in [0.5, 1.0, 2.0] 
                for m in [1.0, 0.9 + 0.001im, 1.1 - 0.001im]
                    for ε in [0.5, 1.0, 2.0]
                        for ngauss in [4, 20, 1000]
                            for ncheb in [2, 3, 4, 10]
                                @test begin
                                    scatterer = TMatrix.Scatterer(r=rev, shape=TMatrix.SHAPE_CHEBYSHEV, refractive_index=m, axis_ratio=ε, n=ncheb)
                                    r2, drr = TMatrix.Wrapper.rsp2(tm, ngauss, rev, ε, ncheb)
                                    r, dr = TMatrix.calc_r(scatterer, ngauss)

                                    r.^2 ≈ r2 && dr ./ r ≈ drr
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end