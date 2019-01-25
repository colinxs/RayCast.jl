module Bresenham

export calc_ranges!

using Distances
using UnsafeArrays
using StaticArrays
using RoboLib.Util: rangebearing2point

#TODO(cxs): decide typing or nah

@inline function cast_heading(I::Type{<:Integer}, x0::Real, y0::Real, heading::Real, maxrange::Real, occupied::Function)
    x1, y1 = rangebearing2point(x0, y0, heading, maxrange)
    cast(I, x0, y0, x1, y1, occupied::Function)
end

@inline function cast_heading(x0::Real, y0::Real, heading::Real, maxrange::Real, occupied::Function)
    cast_heading(Int, x0, y0, heading, maxrange, occupied)
end

@inline cast(x0, y0, x1, y1, maxx, maxy, maxrange, occupied::Function) = cast(Int, x0, y0, x1, y1, occupied::Function)

@inline function cast(I::Type{<:Integer}, x0, y0, x1, y1, occupied::Function)
    cast(round(I, x0), round(I, y0), round(I, x1), round(I, y1), occupied)
end

function cast(x0::I, y0::I, x1::I, y1::I, occupied::Function) where I <: Integer
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)

    xstep = ifelse(x0 < x1, one(I), -one(I))
    ystep = ifelse(y0 < y1, one(I), -one(I))

    #if dy > dx
    #   return _bline_high(x0, y0, x1, y1, dx, dy, xstep, ystep, occupied)
    #else
    #   return _bline_low(x0, y0, x1, y1, dx, dy, xstep, ystep, occupied)
    #end
    steep = dy > dx
    #x0, y0 = ifelse(steep, (y0, x0), (x0, y0))
    #x1, y1 = ifelse(steep, (y1, x1), (x1, y1))
    #xstep, ystep = ifelse(steep, (ystep, xstep), (xstep, ystep))
    #dx, dy = ifelse(steep, (dy, dx), (dx, dy))
    x0, y0, x1, y1, xstep, ystep, dx, dy = ifelse(steep, (y0, x0, y1, x1, ystep, xstep, dy, dx), (x0, y0, x1, y1, xstep, ystep, dy, dx))
    #testoccupied(x,y) = ifelse(steep, occupied(y,x), occupied(x,y))
    #test = ifelse(steep, (x,y)->occupied(y,x), (x,y)->occupied(x,y))
    #if steep
    #    testoccupied(x,y) = occupied(y,x)
    #else
    #    testoccupied(x,y) = occupied(x,y)
    #end
    #@println(testoccupied)

    x, y = _bline_low(x0, y0, x1, y1, dx, dy, xstep, ystep, occupied)
    #x, y = ifelse(steep, (y, x), (x, y))
    return x, y
end

#function _bline_high(x0::I, y0::I, x1::I, y1::I, dx::I, dy::I, xstep::I, ystep::I, occupied::F) where I <: Integer where F <: Function
#    err = convert(I, 2)*dx - dy
#    x, y = x0, y0
#    #y1 = max(one(I), min(y1, maxy))
#    for yi=y0:ystep:y1
#        y = yi
#        #if x < one(I) return x + one(I), y
#        #elseif x > maxx return x - one(I), y
#        #elseif (x-x0)^2 + (y-y0)^2 > maxrange2 return x, y
#        #elseif occupied(x, y) return x, y
#        #end
#        if (err > zero(I))
#            x += xstep
#            err -= convert(I, 2)*dy
#        end
#        err += convert(I, 2)*dx
#    end
#    return x, y
#end
#function _bline_high2(x0::I, y0::I, x1::I, y1::I, dx::I, dy::I, xstep::I, ystep::I, occupied::F) where I <: Integer where F <: Function
#    err = convert(I, 2)*dx - dy # 2dy -dx
#    x, y = x0, y0
#    #y1 = max(one(I), min(y1, maxy))
#    for yi=y0:ystep:y1 #xi=x0:xstep:x1
#        y = yi # x = xi
#        #if x < one(I) return x + one(I), y # x, y+ 1
#        #elseif x > maxx return x - one(I), y
#        #elseif (x-x0)^2 + (y-y0)^2 > maxrange2 return x, y
#        #elseif occupied(x, y) return x, y
#        #end
#        if (err > zero(I))
#            x += xstep # y+=ystep
#            err -= convert(I, 2)*dy # err -dx
#        end
#        err += convert(I, 2)*dx
#    end
#    return x, y
#end

function _bline_low(x0::I, y0::I, x1::I, y1::I, dx::I, dy::I, xstep::I, ystep::I, occupied::F) where I <: Integer where F <: Function
    err = convert(I, 2)*dy - dx
    #x0+=xstep
    x, y, yn = x0, y0, y0
    #println("start")
    #println(err)
    #x1 = max(one(I), min(x1, maxx))
    for xn=x0:xstep:x1
        x = xn
        y = yn
        #println("x, ",x," y, ", y, " err, ", err)
        #if y < one(I) return x, y + one(I)
        #elseif y > maxy return x, y - one(I)
        #elseif (x-x0)^2 + (y-y0)^2 > maxrange2 return x, y
        if occupied(x, y) return x, y end
        #end
        if (err > zero(I))
            yn += ystep
            err -= convert(I, 2)*dx
        end
        err += convert(I, 2)*dy
    end
    return x, y
end

end # module