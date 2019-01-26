module Bresenham

export cast, cast_heading

using StaticArrays
using RoboLib.Util: rangebearing2point

@inline function cast_heading(I::Type{<:Integer}, x0::Real, y0::Real, heading::Real, maxrange::Real)
    x1, y1 = rangebearing2point(x0, y0, heading, maxrange)
    cast(I, x0, y0, x1, y1)
end

@inline function cast_heading(x0::Real, y0::Real, heading::Real, maxrange::Real)
    cast_heading(Int, x0, y0, heading, maxrange)
end

@inline cast(x0, y0, x1, y1) = cast(Int, x0, y0, x1, y1)

@inline function cast(I::Type{<:Integer}, x0, y0, x1, y1)
    cast(round(I, x0), round(I, y0), round(I, x1), round(I, y1))
end

@inline function cast(x0::I, y0::I, x1::I, y1::I, occupied)::SVector{2, I} where I <: Integer
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)

    xstep = ifelse(x0 < x1, one(I), -one(I))
    ystep = ifelse(y0 < y1, one(I), -one(I))

    if dy > dx
       _bline_high(x0, y0, x1, y1, dx, dy, xstep, ystep, occupied)
    else
       _bline_low(x0, y0, x1, y1, dx, dy, xstep, ystep, occupied)
    end
end

function _bline_high(x0::I, y0::I, x1::I, y1::I, dx::I, dy::I, xstep::I, ystep::I, occupied::Function)::SVector{2, I} where I <: Integer
    err = I(2)*dx - dy
    x, xn, y = x0, x0, y0
    for yn=y0:ystep:y1
        x = xn
        y = yn
        if occupied(x, y) return SVector(x, y) end
        if (err > zero(I))
            xn += xstep
            err -= I(2)*dy
        end
        err += I(2)*dx
    end
    return SVector(x, y)
end

function _bline_low(x0::I, y0::I, x1::I, y1::I, dx::I, dy::I, xstep::I, ystep::I, occupied::Function)::SVector{2, I} where I <: Integer
    err = I(2)*dy - dx
    x, y, yn = x0, y0, y0
    for xn=x0:xstep:x1
        x = xn
        y = yn
        if occupied(x, y) return SVector(x, y) end
        if (err > zero(I))
            yn += ystep
            err -= I(2)*dx
        end
        err += I(2)*dy
    end
    return SVector(x,y)
end

end # module