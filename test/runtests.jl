using Test
using RayCast.Bresenham: cast_heading

@testset "RayCast" begin

#Todo(cxs): move Test package to extras

@test begin
    thetas = LinRange(0, 2*pi, 100)
    sidemax = ceil(sqrt(100^2/2))
    rmax = sqrt(2*sidemax^2)
    for theta in thetas
        x, y = cast_heading(0,0,theta,100, (x,y)->false)
        r = sqrt(x^2 + y^2)
        if r > rmax return false end
    end
    return true
end

end




