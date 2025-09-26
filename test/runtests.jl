
using Test
using ASSET


@testset "ASSET" begin
    let T = Float32, dims = (3, 4, 2)
        @testset "CalibratedData" begin
            d = rand(T, dims)
            
            D = CalibratedData(d, d, d, d)
            @test typeof(D) <: CalibratedData{T,length(dims)}
            @test D.d === d
            
            @test size(D) == dims
            @test axes(D) == axes(d)
            @test eltype(D) == T
            @test_nowarn show(IOBuffer(), D)
        end

        @testset "AbstractBkg" begin

            struct MockBkg{R} <: AbstractBkg
                params::R
                regul::Function
            end
            get_bkg(B::MockBkg) = ones(typeof(B.params), dims) .* B.params
            
            B = MockBkg(one(T), x -> x)
            @test typeof(B) <: MockBkg{T}
            @test get_bkg(B) == ones(T, dims)
            @test_throws MethodError B + ones(dims)
            @test_throws MethodError B - ones(dims)
        end
    end
end
