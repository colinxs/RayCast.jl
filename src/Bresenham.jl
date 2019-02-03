module Bresenham

export cast, cast_heading

using StaticArrays
using InteractiveUtils
using RoboLib.Util: rangebearing2point

@inline function cast_heading(I::Type{<:Integer}, x0::Real, y0::Real, heading::Real, maxrange::Real, occupied)
    x1, y1 = rangebearing2point(x0, y0, heading, maxrange)
    cast(I, x0, y0, x1, y1, occupied)
end

@inline function cast_heading(x0::Real, y0::Real, heading::Real, maxrange::Real, occupied)
    cast_heading(Int, x0, y0, heading, maxrange, occupied)
end

@inline cast(x0, y0, x1, y1, occupied) = cast(Int, x0, y0, x1, y1, occupied)

@inline function cast(I::Type{<:Integer}, x0, y0, x1, y1, occupied)
    cast(round(I, x0), round(I, y0), round(I, x1), round(I, y1), occupied)
end

function cast(x0::I, y0::I, x1::I, y1::I, occupied) where I <: Integer
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)

    xstep = ifelse(x0 < x1, one(I), -one(I))
    ystep = ifelse(y0 < y1, one(I), -one(I))

    if dy > dx
       #InteractiveUtils.@code_warntype _bline_high(x0, y0, x1, y1, dx, dy, xstep, ystep, occupied)
       _bline_high(x0, y0, x1, y1, dx, dy, xstep, ystep, occupied)
    else
       _bline_low(x0, y0, x1, y1, dx, dy, xstep, ystep, occupied)
    end
end

_helper(x::Integer, y::Integer, f)::Bool = f(x,y)::Bool


function _bline_high(x0::I, y0::I, x1::I, y1::I, dx::I, dy::I, xstep::I, ystep::I, occupied) where I <: Integer
    err = I(2)*dx - dy
    x, xn, y = x0, x0, y0
    for yn=y0:ystep:y1
        x = xn
        y = yn
        if occupied(x, y) break end
        if (err > zero(I))
            xn += xstep
            err -= I(2)*dy
        end
        err += I(2)*dx
    end
    return SVector{2, I}(x, y)
end

function _bline_low(x0::I, y0::I, x1::I, y1::I, dx::I, dy::I, xstep::I, ystep::I, occupied) where I <: Integer
    err = I(2)*dy - dx
    x, y, yn = x0, y0, y0
    for xn=x0:xstep:x1
        x = xn
        y = yn
        if occupied(x, y) break end
        if (err > zero(I))
            yn += ystep
            err -= I(2)*dx
        end
        err += I(2)*dy
    end
    return SVector{2, I}(x,y)
end

end # module