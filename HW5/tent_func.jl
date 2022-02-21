
"""
    Tent function, the basis function
"""
function UnitTentFunc(x, y)
    
    if x < 0 && y < 0 && x > -1 && y > -1   # bottom left
        return ((x + 1)*(y + 1))/((0 + 1)*(0 + 1))
    end
    if x > 0 && y < 0 && x < 1 && y > -1 # bottm right  
        return (x - 1)*(y + 1)/((0 - 1)*(0 - (-1)))
    end
    if x > 0 && y > 0 && x < 1 && y < 1  # Top right
        return ((x - 1)*(y - 1))/((0 - 1)*(0 - 1))
    end
    if x < 0 && y > 0 && x > - 1 && y < 1 # Top Left 
        return ((x + 1)*(y - 1))/((0 + 1)*(0 - 1))
    end
return 0 end

using Plots
plotlyjs()
"""

"""
function RunThis()
    xs = collect(LinRange(-1, 1, 100))
    as = collect(LinRange(-1, 1, 100))
    x_grid = [x for x = xs for y = as]
    a_grid = [y for x = xs for y = as]
    plot(
        x_grid, 
        a_grid, 
        UnitTentFunc.(x_grid, a_grid), 
        linetype=:path3d) |> display
return end
RunThis();