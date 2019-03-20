module RayCast

# TODO copied from robolib...
@inline function rangebearing2point(x0::Real, y0::Real, bearing::Real, range::Real)
    s, c = sincos(bearing)
    x1 = x0 + range * c
    y1 = y0 + range * s
    return x1, y1
end

#using Logging, LoggingExtras
#const LOGGER = DemuxLogger(FileLogger(joinpath(@__DIR__, "logs/log.txt")))
#global_logger(LOGGER)

include("Bresenham.jl")

end # module
