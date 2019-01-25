using RayCast.Bresenham: cast_heading, cast
using BenchmarkTools
using StaticArrays
using Profile
using InteractiveUtils

@inline function test2(x,y)::Bool
    return false
end
function footest(n, r)
           thetas = collect(LinRange(0, 2*pi, n))
           sc = sincos.(thetas)
           xy = [SVector(trunc(Int, c*r), trunc(Int, s*r)) for (s, c) in sc]
           barbuf = Vector{SVector{2, Int}}(undef, n)
           foobuf = Vector{SVector{2, Int}}(undef, n)
           #@btime bar($thetas, $r, $barbuf)
           #Profile.clear_malloc_data()
           @btime foo($xy, $foobuf)
           return barbuf, foobuf
end

function validate_ranges()
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
            # so that line steps until the endpoint
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


function foo(xy, buf)
    for (i, (x, y)) in enumerate(xy)
        xhit, yhit = cast(Int, 1, 1, x, y, test2)
        buf[i] = SVector(xhit, yhit)
    end
    buf
end

function bar(thetas, r, buf)
    for (i, t) in enumerate(thetas)
        x, y = cast_heading(Int, 1, 1, t, r, (x,y)->false)
        buf[i] = SVector(x, y)
    end
    buf
end