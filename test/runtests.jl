using Test

using RayCast.Bresenham: cast_heading, cast
using RayCast: rangebearing2point
using IterTools
using StaticArrays
using RecursiveArrayTools
using Distributions
using Random

include("util.jl")

TD = TestData{Float64, Int}(136, 50, 50, 200, 1234)

function cast_alloc_test!(buf, x0y0, x1y1, occ)
    alloc = @allocated _cast_alloc_test!(buf, x0y0, x1y1, occ)
    return alloc == 0
end

function _cast_alloc_test!(buf, x0y0, x1y1, occ)
    @inbounds for i in eachindex(x0y0)
        x0, y0 = x0y0[i]
        @inbounds for j in eachindex(x1y1)
            x1, y1 = x1y1[i][j]
            buf[i][j] = cast(x0, y0, x1, y1, occ)
        end
    end
end

function cast_validate_line!(x0y0, x1y1, occ)
    pass = true
    for i in eachindex(x0y0)
        x0, y0 = x0y0[i]
        for j in eachindex(x1y1)
            x1, y1 = x1y1[i][j]
            pass &= _cast_validate_line!(x0, y0, x1, y1, occ)
        end
    end
    pass
end

function _cast_validate_line!(x0::Integer, y0::Integer, x1::Integer, y1::Integer, occ)
    _cast_validate_line!(x0, y0, x1, y1, occ, isequal)
end

# if dealing with floats than we may be off by at most 0.5
# due to rounding errors
function _cast_validate_line!(x0, y0, x1, y1, occ)
    _cast_validate_line!(x0, y0, x1, y1, occ, (x,y)->abs(x-y) <= 0.5)
end

function _cast_validate_line!(x0, y0, x1, y1, occ, comparator)
    # build line eq
    dx = x1 - x0
    dy = y1 - y0
    function errfunc(x,y)
        # Ax + By + C where C==0 b/c no y intercept
        # pos iff below line
        dy*x - dx*y
    end

    # test that the algorithm correctly follows the line
    bounded = true
    # maxerror follows directly from Bresenham line algorithm
    maxerror = max(abs(dx), abs(dy))
    function test(x,y)
        err = errfunc(x,y)
        bounded &= abs(err) <= maxerror
        return false
    end

    # computed endpoints
    xt, yt = cast(x0, y0, x1, y1, occ)

    if !comparator(x1, xt) || !comparator(y1, yt)
        return false
    elseif  !bounded
        return false
    end
    return true
end

## TODO: clean up and refactor for shared tests/benchmarking
#using CuArrays
#using CUDAnative
#using RoboLib.Util: @gridstride
#function cast_cuda_test(buf, x0y0, x1y1, occ)
#    npoints, nrays = length(buf), length(buf[1])
#    nqueries = npoints * nrays
#
#    bufdtype = eltype(buf[1])
#    querydtype = eltype(x0y0)
#
#    # the following are all 1 x nqueries
#    x1y1 = reshape(VectorOfArray(x1y1), nqueries)
#    buf = reshape(VectorOfArray(buf), nqueries)
#    x0y0rep = Vector{querydtype}(undef, nqueries)
#    for i in 1:npoints, j in 1:nrays
#        offset = (i-1) * nrays
#        x0y0rep[offset+j] = x0y0[i]
#    end
#
#    bufd = cu(buf)
#    x0y0repd = cu(buf)
#    x1y1d = cu(x1y1)
#
#    numblocks = ceil(Int, length(bufd))
#    @cuda threads=256 blocks=numblocks cudacast!(bufd, x0y0repd, x1y1d, occ)
#    return bufd == x1y1d
#
#end
#
#function cudacast!(buf, x0y0, x1y1, occ)
#    @gridstride length(buf) i begin
#        x0, y0 = x0y0[i]
#        x1, y1 = x1y1[i]
#        buf[i] = cast(x0, y0, x1, y1, occ)
#    end
#    nothing
#end
#
@testset "Bresenham" begin
    @inferred cast(1,1,10,10,TD.fn)
    @inferred cast(1,1,10,10,TD.ftor)
    @inferred cast(1.0,1,10,10,TD.fn)
    @inferred cast(1.0,1,10,10,TD.ftor)
    @inferred cast(Int, 1.0,1,10,10,TD.fn)
    @inferred cast(Int, 1.0,1,10,10,TD.ftor)

    @inferred cast_heading(1,1,pi/2,10, TD.fn)
    @inferred cast_heading(1,1,pi/2,10, TD.ftor)
    @inferred cast_heading(Int,1,1,pi/2,10, TD.fn)
    @inferred cast_heading(Int,1,1,pi/2,10, TD.ftor)
    @inferred cast_heading(Int,1.0,1,pi/2,10, TD.fn)
    @inferred cast_heading(Int,1.0,1,pi/2,10, TD.ftor)

    @test_broken cast_alloc_test!(TD.xyhit, TD.x0y0_i, TD.x1y1_i, TD.fn)
    @test cast_alloc_test!(TD.xyhit, TD.x0y0_i, TD.x1y1_i, TD.ftor)
    @test cast_validate_line!(TD.x0y0, TD.x1y1, TD.ftor)
    #@test cast_cuda_test(TD.xyhit, TD.x0y0_i, TD.x1y1_i, TD.ftor)
end




