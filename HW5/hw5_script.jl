"""
    The cofficient of finite difference of 2D 5 points laplacian without
    the h^2 multiplier and assume equally spaced gridpoints. 
    
"""
function Laplacian5Stencil()
    stencil = Vector()
    push!(stencil, (0, 0, -4))
    push!(stencil, (-1, 0, 1))
    push!(stencil, (0, -1, 1))
    push!(stencil, (1, 0, 1))
    push!(stencil, (0, 1, 1))

return stencil end

"""
    The cofficient of finite difference of 2D 9 points laplacian without
    the h^2 multiplier and assume equally spaced gridpoints. 
"""
function Laplacian9Stencil()
    stencil = Vector()
    push!(stencil, (-1, -1, 1/6))
    push!(stencil, (-1, 1, 1/6))
    push!(stencil, (1, -1, 1/6))
    push!(stencil, (1, 1, 1/6))
    push!(stencil, (0, -1, 4/6))
    push!(stencil, (0, 1, 4/6))
    push!(stencil, (1, 0, 4/6))
    push!(stencil, (-1, 0, 4/6))
    push!(stencil, (0, 0, -20/6))

return stencil end

"""
    Make a Linear system representin the Poisson Problem. 
    Grid point: 
        Equally space on both x and y direction on the domain of [0, 1] x [0, 1]
    Boundary Conditions: 
        4 pices of boundary conditions are keyword paremter, should be passed in
        as a scalar function. 
    Keyworld Parameters: 
        * u_left, u_right, u_top, u_bottom are the function modeling the bouncary
        conditions. 
        * f is the righthand side function. 
        * stencil: It's a stencil generative function. 
    
    
"""
function MakeLaplacianSystem(
    m::Int64;
    f=nothing,       # keyword argument
    u_left=nothing, 
    u_right=nothing, 
    u_top=nothing, 
    u_bottom=nothing, 
    stencil::Function=Laplacian5Stencil
)
    @assert m > 2 "m should be larger than 2. "
    zeroFunc(x, y) = 0
    returnsMatrixOnly = isnothing(u_top) &&
        isnothing(u_bottom) && 
        isnothing(u_left) && 
        isnothing(u_right) && 
        isnothing(f) 
    
    # if every keyword argument are nohting, then we just need the 
    # laplacian matrix 
        
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

    # Construct using Stencil 
    stencil = stencil()
    for i in 1:m, j in 1:m
        for (k, l, v) in stencil
            AddCoeff((i, j), (i + k, j + l), v/h^2)
        end
        rhsMod[Tl(i, j)] += f(i*h, j*h)
    end

    if returnsMatrixOnly
        return sparse(RowIdx, ColIdx, Vals)
    end
    
return sparse(RowIdx, ColIdx, Vals), rhsMod end


"""
    Return the interior points of a partitioning of the grid points
    in a vector of tuples of x, y if the parameter u is not given. 
    if it's given then it will just put these tuples of values into the 
    given u function, values are returned as a vector of natural ordering. 

    parameter: 
        m:: An integer larger than 2. 
"""
function Grid2Vec(m::Int64, u=nothing)
    v = Vector{Float64}()
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
    the (i,j) element is (j*h, i*h), so it's the 2d quadrant rotated to the right 
    by 180

    * Plotting of the matrix using heatmap will recover the function's plot 
    over [0, 1] x [0, 1] like, exactly as what you expected because 
    julia plots it upside down. 

"""
function Vec2Grid(
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
    error("Not Yet implemented yet. ")
return end


### Hands on testing part ----------------------------------------------------------------------------------------------

using SparseArrays, Plots

"""
    Basic test. 
"""
function BasicTests()
    f = (x, y) -> x^2 + y^2
    m = 3
    M, b = MakeLaplacianSystem(
        3,
        f=nothing, u_right=(x)->1, u_top=(x) -> 1, u_bottom =(x)-> 1, u_left=(x) ->1
    )
    @info "Carrying out some basic testing: "
    M |> display 
    b |> display
    Vec2Grid(m, Grid2Vec(m, (x, y)-> x + y)) |> heatmap |> display

    @info "Testing whether the 9 points Laplacian works"
    m = 63
    A, b = MakeLaplacianSystem(
        m, f=f,
        stencil=Laplacian9Stencil
    )
    println("M matrix: ")
    Vec2Grid(
        m, 
        A\b
    ) |> heatmap |> display

return end

BasicTests()

function Problem1()
    # super fine grid point solution come first. 
    M = 2^10 - 1
    f(x, y) = x^2 + y^2
    u_bc(x) = 1
    A, b = MakeLaplacianSystem(
        M,
        f=f,
        u_left=u_bc, u_right=u_bc, u_top=u_bc, u_bottom=u_bc
    )

    @info "The system we are solving is: "
    A |> display
    b |> display
    uVeryFine = Vec2Grid(M, A\b)
    uVeryFine |> heatmap |> display

    # estimating the error using super fine grid as a reference
    Errors = Vector()       # L2 error over the grid
    gridWidth = Vector()    
    for m in 2 .^ (collect(2:8)) .- 1
        h = 1/(m + 1)
        A, b = MakeLaplacianSystem(
            m, 
            f=f,
            u_left=u_bc, u_right=u_bc, u_top=u_bc, u_bottom=u_bc
        )
        uCoarse = Vec2Grid(m, A\b)
        skip = convert(Int64, (M + 1)/(m + 1) |> floor)
        push!(Errors, h*norm(uVeryFine[skip:skip:end, skip:skip:end] - uCoarse))
        push!(gridWidth, h)
    end
    # Printing things out: 
    @info "These are a print out for the error vectors and gird width. "
    Errors |> display
    gridWidth |> display
    fig = plot(
        gridWidth .|> log2,
        Errors .|> log2, 
        title="Log2 vs Log2 Error", 
        label="5 p stencil", 
        legend=:bottomright
    )
    plot!(fig, gridWidth .|> log2, gridWidth.^2 .|> log2, label="reference h^2")
    xlabel!(fig, "log2(h)")
    ylabel!(fig, "log2(E)")
    fig |> display
    savefig(fig, "p1_fig.png")
    # Estimating rate of convergence. 

return end

# Problem1()

"""

"""
function Problem2()
    function SolveWithDeferredCorrection(M::Int64)
        # super fine grid point solution come first. 
        M = 2^10 - 1
        f(x, y) = x^2 + y^2
        u_bc(x) = 1
        A, b = MakeLaplacianSystem(
            M,
            f=f,
            u_left=u_bc, u_right=u_bc, u_top=u_bc, u_bottom=u_bc
        )
        # deferred correlations
        LaplacianF(x, y) = 4
        b += (1/(M + 1)^2)*Grid2Vec(M, LaplacianF)/12
    return Vec2Grid(M, A\b) end
    uVeryFine = SolveWithDeferredCorrection(2^10 - 1)
    @info "uVeryFine 9 points stencils solution looks like: "
    uVeryFine |> display

end

Problem2()