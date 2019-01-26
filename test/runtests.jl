using Test

using RayCast.Bresenham: cast_heading, cast
using StaticArrays

@inline function occupied(x,y)::Bool
    return false
end

function alloc_test()
    n, r = 360, 100
    thetas = collect(LinRange(0, 2*pi, n))
    sc = sincos.(thetas)
    xy = [SVector(trunc(Int, c*r), trunc(Int, s*r)) for (s, c) in sc]
    buf = [SVector(0,0) for _ in 1:n]
    alloc = @allocated _alloc_test(xy, buf)
    return iszero(alloc)
end

function _alloc_test(xy, buf)
    for i in 1:length(xy)
        x, y = xy[i]
        buf[i] = cast(1, 1, x, y, occupied)
    end
end

function validate_line()
    r = 10
    for (i, theta) in enumerate(LinRange(0,2*pi, 36))
        s, c = sincos(theta)

        # target endpoints
        x = round(Int, r * c)
        y = round(Int, r * s)

        # build line eq
        dx = x
        dy = y
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
        xn, yn = cast(0,0,x,y, test)

        if x != xn || y != yn
            return false
        elseif  !bounded
            return false
        end
    end
    return true
end

@testset "Bresenham" begin

    @test validate_line()
    @test alloc_test()
end




