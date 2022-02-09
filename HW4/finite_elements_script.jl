function MidPointRuleSingleInterval(f::Function, a, b)
return (b - a)*f((a/2 + b/2)) end



function TrapzRuleSingleInterval(f::Function, a, b)
return (b - a)*(f(a)/2 + f(b)/2) end



function TwoPointsGaussQuadratureSingleInterval(f::Function, a, b)
    m = a/2 + b/2
    l = m - (1/sqrt(3))*(m - a)
    r = m + (1/sqrt(3))*(b - m)
return 0.5*(b - a)*(f(l) + f(r)) end


"""
    A specific steady state BVP finite element problem. 
    -(r(x)u'(x))' = f(x), x in [0, 1], u(0) = alpha, u(1) = beta
    Using:  
        - Any user defined mesh grid in 1D. 
        - Piecewiese Linear Basis function only. 

    * Field Explanations
        - 
"""
mutable struct FiniteElementProblem{T<:Number}
    r::Function         # heat conductivity
    alpha::T            # left boundary condition, dirichlet
    beta::T             # right boundary condition, dirichlet
    f::Function         # RHS function
    integral_rule:: Function
                        # an implementation of taking integral. 

    x_grids::AbstractVector{T}
                        # Number of grid inner grid points for the problem, 
                        # includes the boundary. 
    max_h::T               # Maximum meshgrid spacing. 

    A::AbstractMatrix{T}# Discrete Linear Operator for F.E. 
    rhs::Vector{T}                 # The rhs function INCLUDING the boundary conditions. 
    c::Vector{T}                   # Weights on Basis Functions. 


    function FiniteElementProblem(
        r::Function, 
        f::Function, 
        alpha, 
        beta, 
        x_grids::AbstractVector,
        integral_rule::Function=MidPointRuleSingleInterval
    )
        this = new{AbstractFloat}()
        this.r = r
        this.f = f
        this.integral_rule = integral_rule
        this.alpha = alpha
        this.beta = beta
        this.x_grids = x_grids
        this.max_h = maximum(x_grids[2:end] - x_grids[1:end - 1])
        
        @assert minimum(x_grids[2:end] - x_grids[1:end - 1]) > 0 "Meshgrid Error."
        
        MakeDiscreteOperator!(this)    
        MakeRHS!(this)

    return this end

end

"""
    Re-estrablish the discretized operator in the field 
    of the struct using the current parameters in the struct. 
"""
function MakeDiscreteOperator!(this::FiniteElementProblem)
    # make diagonals
    diagonal = Vector()
    subdiagonal = Vector()
    for (xl, xm, xr) in zip(
        this.x_grids[1:end - 2], 
        this.x_grids[2:end - 1],
        this.x_grids[3:end]
    )
        push!(
            diagonal, 
            this.integral_rule(this.r, xl, xm)/(xm - xl)^2 +
            this.integral_rule(this.r, xm, xr)/(xr - xm)^2
        )
    end

    # make sub-diagonals
    for (xl, xr) in zip(
        this.x_grids[2:end - 2], 
        this.x_grids[3:end - 1]
    )
        push!(
            subdiagonal, 
            -this.integral_rule(this.r, xl, xr)/(xr - xl)^2
        )

    end
    this.A = SymTridiagonal(diagonal, subdiagonal)
return end


"""
    establish the new RHS for the finite element equation using the 
    current parameters defined in the struct. 
"""
function MakeRHS!(this::FiniteElementProblem)
    rhsVec = Vector()
    func(x, z) = this.f(x)*(x - z)
    for (xl, xm, xr) in zip(
        this.x_grids[1:end - 2], 
        this.x_grids[2:end - 1],
        this.x_grids[3:end]
    )
        push!(
            rhsVec, 
            this.integral_rule((x) -> func(x, xl), xl, xm)/(xm - xl) - 
            this.integral_rule((x) -> func(x, xr), xm, xr)/(xr - xm)
        )

    end
    this.rhs = rhsVec
return end

"""
    Given a series of points to query the solution value at, 
    we solve this the system. 
"""
function Solve(this::FiniteElementProblem)
    if !isdefined(this, :c)
        this.c = this.A\this.rhs
    end
return copy(vcat(this.alpha, this.c, this.beta)) end




## Using this code to solve some problems in HW. 

using LinearAlgebra, Plots, Logging

function BasicTests()
    xgrid = LinRange(0, 1, 11)
    r(x) = (1 + x^2)
    alpha, beta = 0, 0
    delta = 0.1
    f(x) = 2*(3x^2 - x + 1)

    P = FiniteElementProblem(r, f, alpha, beta, xgrid)
    fig = plot(xgrid, Solve(P), label="F.E Soln")
    plot!(fig, LinRange(0, 1, 1000), (x)-> x*(1 - x), label="u(x)") |> display
    
return P end

P = BasicTests()

function Problem1PartC()
    gridPoints = 0:10
    Errors = Vector{AbstractFloat}()
    GridSize = Vector{AbstractFloat}()
    for m in .^(collect(gridPoints)) .+ 1
        xgrid = LinRange(0, 1, m)
        
    end
return end
