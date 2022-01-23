# For HW 2. Problem 1 
# similate the non-homogenous heat conductivity problem. 

using SparseArrays

"""
    (κ(x)u'(x))' = f(x), with u(0) = 1, u'(1) = 0 
    1d BVP.
    Parameter m: 
        It defines the grids. 
        Grids: x0, x1, ... xm, x(m + 1)
        h: 1/(m + 1)
        u0 = 0, which simplifies the system to a (m + 1)×(m + 1) system. 
        So, m is the number of interior points for the domain of the function. 
    Parameter c: 
        A function for scalar that is the thermal conductivity on the rod. 
    Returns: 
        The matrix A to solve and an additional vector for modifying the 
        RHS, which is the boundary conditions.
"""
function FiniteDiffMatrix(m::Int64, c::Function)
    h = 1/(m + 1)
    dict = Dict{Tuple{Int64, Int64}, Float64}()
    # first row
    
    dict[1, 1] = (1/h^2)*(c(3*h/2) - c(h/2))
    dict[1, 2] = (1/h^2)*c(3*h/2)
    for i in 2:m
        dict[i, i - 1] = c((i - 1/2)*h)/h^2
        dict[i, i] = - c((i + 1/2)*h)/h^2 - c((i - 1/2)*h)/h^2
        dict[i, i + 1] = c((i + 1/2)*h)/h^2
    end
    # Last row
    dict[m + 1, m] = -4/h^2
    dict[m + 1, m - 1] = 1/h^2
    dict[m + 1, m + 1] = 3/h^2 
    # coordinate format. 
    CooFormatx = Vector{Float64}()
    CooFormaty = Vector{Float64}()
    CooFormatval = Vector{Float64}()
    for (K, V) in dict
        push!(CooFormatx, K[1])
        push!(CooFormaty, K[2])
        push!(CooFormatval, V)
    end
    # RHS fixes by boundaries conditions. 
    b = zeros(m + 1)
    b[1] = c(h/2)/h^2
    # sparse matrix. 
    return sparse(CooFormatx, CooFormaty, CooFormatval), b
end


"""
    Solve the system: (κ(x)u'(x))' = f(x), with u(0) = 1, u'(1) = 0 
    for c(x) = (1 - x)^2. 
"""
function SolveFor(f::Function, m::Int64)
    c(x) = (1 + x)^2
    A, b = FiniteDiffMatrix(m, c)
    display(A)
    display(b)
    RHS = f.(LinRange(0, 1, m + 2)[2:end]) + b
    u = A\RHS
    return vcat(1, u)
end


# Basic testing
using Plots
m = 63
f = (x) -> 2*(3x^2 - 2x + 1)
u = SolveFor(f, m)
plot(LinRange(0, 1, m + 2), u)

# Verify the Errors rate
Errors = Vector{Float64}()
Hs = Vector{Float64}()
for m in 2.0 .^ collect(4:12) .- 1
    h = 1/(m + 1)
    u = SolveFor(f, m)
    push!(Hs, norm())
end