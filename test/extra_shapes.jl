@testset "Extra shapes" begin
    λ = 1
    m = 1311 // 1000

    @testset "Bicone" begin
        @testset "r = $r, h = $h" for (r, h) in [(1, 1), (1, 5 // 4), (1, 4 // 5)]
            @test begin
                scatterer = TMatrix.Bicone(r = r, h = h, m = m, λ = λ)
                TT = TMatrix.calc_tmatrix!(scatterer)
                _, _, ω = TMatrix.cross_section(TT, scatterer.λ)
                abs(ω - 1) < 1e-5
            end
        end
    end

    @testset "Capsule" begin
        @testset "r = $r, h = $h" for (r, h) in [(1, 2), (1 // 2, 3 // 2), (1 // 2, 2 // 3)]
            @test begin
                scatterer = TMatrix.Capsule(r = r, h = h, m = m, λ = λ)
                TT = TMatrix.calc_tmatrix!(scatterer)
                _, _, ω = TMatrix.cross_section(TT, scatterer.λ)
                abs(ω - 1) < 1e-4
            end
        end
    end
end
