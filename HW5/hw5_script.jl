
"""
    Make the laplacian for the m x m grid, together with the RHS 
    over the domain of [0, 1] x [0, 1]. 
    ### Parameter: 
        m: the number of interior grid points, 0 and m + 1 are the boundary grid point. 
        if m is given then the grid width would be 1/(m + 1)

"""
function Make2D5PointsLaplacianMatrixAndRHS(
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
    
    function AddCoeff(TupleIdx1, TupleIdx2, coef)
        i1, j1 = TupleIdx1
        i2, j2 = TupleIdx2
        # interior
        if i2 >= 1 && i2 <= m && j2 >= 1 && j2 <= m
            push!(RowIdx, Tl(i1, j1))
            push!(ColIdx, Tl(i2,j2))
            push!(Vals, coef)
        # boundary. 
        elseif i2 == 0 || i2 == m + 1 # left or right B.C
            if i2 == 0 # left
                rhsMod[Tl(i1, j1)] -= u_left(j2*h)*coef
            else      # right
                rhsMod[Tl(i1, j1)] -= u_right(j2*h)*coef
            end
        elseif j2 == 0 || j2 == m + 1 # top or bottom
            if j2 == 0 # bottom 
                rhsMod[Tl(i1, j1)] -= u_bottom(i2*h)*coef
            else       # top
                rhsMod[Tl(i1, j1)] -= u_top(i2*h)*coef
            end
        else 
            error("this should not happen, I don't expect this")
        end

    return end
    
    # implement stencils
    for i in 1:m, j in 1:m
        AddCoeff((i, j), (i, j), -4/h^2)
        AddCoeff((i, j), (i - 1, j), 1/h^2)
        AddCoeff((i, j), (i, j - 1), 1/h^2)
        AddCoeff((i, j), (i + 1, j), 1/h^2)
        AddCoeff((i, j), (i, j + 1), 1/h^2)
        rhsMod[Tl(i, j)] += f(i*h, j*h)
    end

return sparse(RowIdx, ColIdx, Vals), rhsMod end


function Make2D9PointLaplacianAndRHS(
    m::Int64;
    f=nothing, # keyword argument
    u_left=nothing, 
    u_right=nothing, 
    u_top=nothing, 
    u_bottom=nothing
)
    error("I haven't implement this part yet. ")
return end


"""
    Return the interior points of a partitioning of the grid points
    in a vector of tuples of x, y if the parameter u is not given. 
    if it's given then it will just put these tuples of values into the 
    given u function, values are returned as a vector of natural ordering. 
    parameter: 
        m:: An integer larger than 2. 
"""
function GridPointNaturallyOrdered(m::Int64, u=nothing)
    v = Vector()
    h = 1/(1 + m)
    if isnothing(u)
        for i in 1:m, j in 1:m
            push!(v, (i*h, j*h))
        end
        return v
    end

    for i in 1:m, j in 1:m
        push!(v, u(i*h, j*h))
    end
return v end

"""
    convert the naturally ordered vector into a grid, including the 
    boundary conditions. 
    
    * Returns a matrix that is literally the function values over the grid points. 
    the (i,j) element is (j*h, i*h). 
    * Plotting of the matrix will recover the function's plot over [0, 1] x [0, 1]
"""
function NaturallyOrderedToGrid(
    m, v;
    bc_left=nothing,
    bc_right=nothing,
    bc_top=nothing, 
    bc_bottom=nothing
)
    g = zeros(m, m)
    for i in 1:m, j in 1:m
        g[j, i] = v[(j - 1)*m + i]
    end
    if isnothing(bc_left) && isnothing(bc_right) && isnothing(bc_top) && isnothing(bc_bottom)
       return g 
    end
    @error("Not yet implemented.")
    gWithBC = zeros(m + 2, m + 2)
    gWithBC[2:end - 1, 2:end - 1] = g 
    
return gWithBC end

"""
    Performs conjugate gradient to solve the sparse linear effectively. The implementation of this method 
    is faster than writting out an FFT based Poisson solver. 
"""
function SimpleConjugateGradient(A::AbstractMatrix, b::AbstractVector; epsilon=1e-10)


return end


### Hands on testing part ----------------------------------------------------------------------------------------------

using SparseArrays, Plots

"""
    Basic test. 
"""
function BasicTests()
    f = (x, y) -> 1
    M, b = Make2D5PointsLaplacianMatrixAndRHS(3,
        f = nothing, u_right=(x)->1, u_top=(x) -> 1, u_bottom =(x)-> 1, u_left=(x) ->1
    )
    M |> display 
    b |> display
    NaturallyOrderedToGrid(m, GridPointNaturallyOrdered(m, (x, y)-> x + y)) |> heatmap |> display
return end


function Problem1()
    # super fine grid point solution come first. 
    M = 2^10
    f(x, y) = x^2 + y^2
    u_bc(x) = 1
    A, b = Make2D5PointsLaplacianMatrixAndRHS(
        M, 
        f=f,
        u_left=u_bc, u_right=u_bc, u_top=u_bc, u_bottom=u_bc
    )

    @info "The system we are solving is: "
    A |> display
    b |> display
    uVeryFine = NaturallyOrderedToGrid(M, A\b)
    uVeryFine |> heatmap |> display

    # estimating the error using super fine grid as a reference
    Errors = Vector()       # L2 error over the grid
    gridWidth = Vector()    
    for m in 2 .^ collect(2:8)
        h = 1/(m + 1)
        A, b = Make2D5PointsLaplacianMatrixAndRHS(
            m, 
            f=f,
            u_left=u_bc, u_right=u_bc, u_top=u_bc, u_bottom=u_bc
        )
        uCoarse = NaturallyOrderedToGrid(m, A\b)
        skip = convert(Int64, M/m |> floor)
        println("skip is: $(skip)")
        push!(Errors, norm(uVeryFine[1:skip:end, 1:skip:end] - uCoarse))
        push!(gridWidth, h)
    end
    Errors |> display
    gridWidth |> display
    fig = plot(gridWidth .|> log2, Errors .|> log2)
    plot!(fig, gridWidth .|> log2, gridWidth.^2 .|> log2)
    fig |> display
return end

Problem1()

