module Bresenham

export calc_ranges!

using Distances
using UnsafeArrays
using StaticArrays

@inline binary_occupied(x::Integer, y::Integer, grid::Matrix{Bool})::Bool = @inbounds grid[x, y]

abstract type RayData end

# TODO: separate params from data with functor object

struct BHamData{I<:Integer} <: RayData
    xy::Vector{SVector{2, I}}
    maxx::I
    maxy::I
    maxrange::Float64
    occupied::Function

    function BHamData{I}(nrays::Integer, maxx, maxy, maxrange, occupied) where I<:Integer
        return new{I}(Vector{SVector{2, I}}(undef, nrays), maxx, maxy, maxrange, occupied)
    end
end
BHamData(nrays::Integer, maxx, maxy, maxrange, occupied) = BHamData{Int}(nrays, maxx, maxy, maxrange, occupied)

@inline function calc_hit_heading(I::Type{<:Integer}, x0::Real, y0::Real, heading::Real, maxrange::Real, maxx::Real, maxy::Real, occupied::Function)
    x1, y1 = range2point(x0, y0, heading, maxrange)
    calc_hit(round(I, x0), round(I, y0), round(I, x1), round(I, y1), round(I, maxx), round(I, maxy), occupied)
end

@inline function calc_hit_heading(x0::Real, y0::Real, heading::Real, maxrange::Real, maxx::Real, maxy::Real, occupied::Function)
    calc_hit_heading(Int, x0, y0, heading, maxrange, maxx, maxy, occupied)
end

@inline function range2point(x0::Real, y0::Real, heading::Real, max_range::Real)
    @fastmath x1 = x0 + max_range * cos(heading)
    @fastmath y1 = y0 + max_range * sin(heading)
    return x1, y1
end

function bline_high(x0::I, y0::I, x1::I, y1::I, dx::I, dy::I, xstep::I, ystep::I, maxx::I, maxy::I, occupied::F) where I <: Integer where F <: Function
    err = convert(I, 2)*dx - dy
    x, y = x0, y0
    y1 = max(one(I), min(y1, maxy))
    for yi=y0:ystep:y1
        y = yi
        if x < one(I) return x + one(I), y
        elseif x > maxx return x - one(I), y
        elseif occupied(x, y) return x, y
        end
        if (err > zero(I))
            x += xstep
            err -= convert(I, 2)*dy
        end
        err += convert(I, 2)*dx
    end
    return x, y
end


function bline_low(x0::I, y0::I, x1::I, y1::I, dx::I, dy::I, xstep::I, ystep::I, maxx::I, maxy::I, occupied::F) where I <: Integer where F <: Function
    err = convert(I, 2)*dy - dx
    x, y = x0, y0
    x1 = max(one(I), min(x1, maxx))
    for xi=x0:xstep:x1
        x = xi
        if y < one(I) return x, y + one(I)
        elseif y > maxy return x, y - one(I)
        elseif occupied(x, y) return x, y
        end
        if (err > zero(I))
            y += ystep
            err -= convert(I, 2)*dx
        end
        err += convert(I, 2)*dy
    end
    return x, y
end


function calc_hit(x0::I, y0::I, x1::I, y1::I, maxx::I, maxy::I, occupied::F)::Tuple{I, I} where I <: Integer where F <: Function
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)

    xstep = ifelse(x0 < x1, one(I), -one(I))
    ystep = ifelse(y0 < y1, one(I), -one(I))

    if dy > dx
        x1, y1 = bline_high(x0, y0, x1, y1, dx, dy, xstep, ystep, maxx, maxy, occupied)
    else
        x1, y1 = bline_low(x0, y0, x1, y1, dx, dy, xstep, ystep, maxx, maxy, occupied)
    end
    return x1, y1
end


#function calc_ranges!(data::BHamData, x0y0x1y1::Array{<:Integer, 2})
#    maxy, maxx = size(data.grid)
#    #nthreads = Threads.nthreads()
#    #threadrange = Distributed.splitrange(size(x0y0x1y1, 2), nthreads)
#    #Threads.@threads for i=1:nthreads
#        #@inbounds for col=threadrange[i]
#        @inbounds for col=1:size(x0y0x1y1, 2)
#            x0, y0, x1, y1 = uview(x0y0x1y1, :, col)
#            data.xy[:, col] .= calc_hit(x0, y0, x1, y1, maxx, maxy, data.occupied, data.grid)
#        end
#    #end
#end
#
function calc_ranges!(data::BHamData, x0::Real, y0::Real, headings::Vector{<:Real})
    #@inbounds Threads.@threads for col=1:length(headings)
    for col=1:length(headings)
        x, y = calc_hit_heading(x0, y0, headings[col], data.maxrange, data.maxx, data.maxy, data.occupied)
        data.xy[col] = SVector(x, y)
    end
end


#function calc_ranges!(data::BHamData, x0y0heading::Array{<:Real, 2})
#    maxy, maxx = size(data.grid)
#    #nthreads = Threads.nthreads()
#    #threadrange = Distributed.splitrange(size(x0y0x1y1, 2), nthreads)
#    #Threads.@threads for i=1:nthreads
#        #@inbounds for col=threadrange[i]
#        @inbounds for col=1:size(x0y0heading, 2)
#            x0, y0, heading = uview(x0y0heading, :, col)
#            data.xy[:, col] .= calc_hit_heading(x0, y0, heading, data.max_range, maxx, maxy, data.occupied, data.grid)
#        end
#    #end
#end

end # module