
"""
    Make the laplacian for the m x m grid, together with the RHS 
    over the domain of [0, 1] x [0, 1]. 
    Parameter: 
        m, the number of interior grid points, 0 and m + 1 are the boundary grid point. 
        if m is given then the grid width would be 1/(m + 1)

"""
function Make2D5PointsLaplacianMatrix(
    m::Int64;
    f=nothing,       # keyword argument
    u_left=nothing, 
    u_right=nothing, 
    u_top=nothing, 
    u_bottom=nothing
)
    @assert m > 2 "m should be larger than 2. "
    zeroFunc(x, y) = 0
    
    u_left = isnothing(u_left) ? (x) -> 0 : u_left
    u_right = isnothing(u_right) ? (x) -> 0 : u_right
    u_top = isnothing(u_top) ? (x) -> 0 : u_top
    u_bottom = isnothing(u_bottom) ? (x) -> 0 : u_bottom
    f = isnothing(f) ? zeroFunc : f

    Tl(i, j) = i + (j - 1)*m
    # make the matrix
    h = 1/(m + 1)              
    rhsMod = zeros(m^2)
    RowIdx = Vector{Float64}()
    ColIdx = Vector{Float64}()
    Vals   = Vector{Float64}()
    
    function AddCoeff(TupleIdx1, TupleIdx2, Val)
        i1, j1 = TupleIdx1
        i2, j2 = TupleIdx2
        # interior
        if i2 >= 1 && i2 <= m && j2 >= 1 && j2 <= m
            push!(RowIdx, Tl(i1, j1))
            push!(ColIdx, Tl(i2,j2))
            push!(Vals, Val)
            return
        end
        # boundary. 
        rhsMod[Tl(i1, j1)] = f(i1*h, j1*h)  
        if i2 == 0 || i2 == m + 1 # left or right B.C
            if i2 == 0 # left
                rhsMod[Tl(i1, j1)] -= u_left(j1*h)
            else      # right
                rhsMod[Tl(i1, j1)] -= u_right(j1*h)
            end
        elseif j2 == 0 || j2 == m + 1 # top or bottom
            if j2 == 0 # bottom 
                rhsMod[Tl(i1, j1)] -= u_bottom(i1*h)
            else       # top
                rhsMod[Tl(i1, j1)] -= u_top(i1*h)
            end
        end
    return end
    
    # implement stencils
    for i in 1:m, j in 1:m
        AddCoeff((i, j), (i, j), 4/h^2)
        AddCoeff((i, j), (i - 1, j), 1/h^2)
        AddCoeff((i, j), (i, j - 1), 1/h^2)
        AddCoeff((i, j), (i + 1, j), 1/h^2)
        AddCoeff((i, j), (i, j + 1), 1/h^2)
    end

return sparse(RowIdx, ColIdx, Vals), rhsMod end

using SparseArrays
f = (x, y) -> 1
M, b = Make2D5PointsLaplacianMatrix(3, 
    f=nothing, u_right=(x)->1, u_top=nothing
)
