using RayCast.Bresenham: cast_heading, cast
using StaticArrays
using Distributions
using Random
using RayCast: rangebearing2point
#using RoboLib.Util: @gridstride, rangebearing2point
#using CuArrays
#using CUDAnative


struct TestFunctor end
@inline (t::TestFunctor)(x,y) = false

struct TestFunctorArr{A} a::A end
@inline (m::TestFunctorArr)(x,y) = m.a[x, y]

@inline testfn(x, y) = false

@inline testfnarr(x, y, arr) = false

struct TestData{T, I, FTOR, FTORA, FN, FNA, RNG}
    x0y0::Vector{SVector{2, T}}
    x1y1::Vector{Vector{SVector{2, T}}}
    x0y0_i::Vector{SVector{2, I}}
    x1y1_i::Vector{Vector{SVector{2, I}}}

    theta::Vector{T}
    range::T

    ftor::FTOR
    ftorarr::FTORA
    fn::FN
    fnarr::FNA

    rng::RNG
    xyhit::Vector{Vector{SVector{2, I}}}

    function TestData{T, I}(nrays, npoints, range, dim, seed) where {T,I}
        rng = Random.MersenneTwister(1234)
        d = Uniform(range + 1, dim - range)
        theta = LinRange(0, 2*pi, nrays)

        x0y0 = [SVector{2, T}(rand(rng, d), rand(rng, d)) for _ in 1:npoints]

        x0y0_i = similar(x0y0, SVector{2, I})
        x1y1 = [similar(x0y0, nrays) for _ in 1:npoints]
        x1y1_i = [similar(x0y0_i, nrays) for _ in 1:npoints]

        for i in eachindex(x0y0)
            x0, y0 = x0y0[i]
            x0y0_i[i] = SVector{2, I}(round(I, x0), round(I, y0))
            for j in eachindex(theta)
                th = theta[j]
                x1, y1 = rangebearing2point(x0, y0, th, range)

                x1y1[i][j] = SVector{2, T}(x1, y1)
                x1y1_i[i][j] = SVector{2, I}(round(I, x1), round(I, y1))
            end
        end
        ftor = TestFunctor()
        fn = testfn

        arr = zeros(Bool, dim, dim)
        ftorarr = TestFunctorArr(arr)
        fnarr = let arr = arr
            (x, y)->arr[x, y]
        end

        xyhit = deepcopy(x1y1_i)

        new{T, I, typeof(ftor), typeof(ftorarr), typeof(fn), typeof(fnarr), typeof(rng)}(x0y0, x1y1, x0y0_i, x1y1_i, theta, range, ftor, ftorarr, fn, fnarr, rng, xyhit)
    end
end

