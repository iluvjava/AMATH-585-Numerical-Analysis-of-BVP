"""
    (c(x)u'(x))' = f(x), with u(0) = 1, u'(1) = 0 
    1d BVP.
    Parameter m: 
        It defines the grids. 
        Grids: x0, x1, ... xm, x(m + 1)
        h: 1/(m + 1)
        u0 = 0, which simplifies the system to a (m + 1)x(m + 1) system. 
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
    dict[1, 1] = -(c(3*h/2) + c(h/2))/h^2
    dict[1, 2] = c(3*h/2)/h^2
    for i in 2:m
        dict[i, i - 1] = c((i - 1/2)*h)/h^2
        dict[i, i]     = -(c((i + 1/2)*h) + c((i - 1/2)*h))/h^2
        dict[i, i + 1] = c((i + 1/2)*h)/h^2
    end
    
    # Last row
    dict[m + 1, m]     = (c((m + 3/2)*h) + c((m + 1/2)*h))/h^2
    dict[m + 1, m + 1] = -(c((m + 3/2)*h) + c((m + 1/2)*h))/h^2
    
    # coordinate format. 
    CooFormatx = Vector{Float64}()
    CooFormaty = Vector{Float64}()
    CooFormatval = Vector{Float64}()
    for (K, V) in dict
        push!(CooFormatx, K[1])
        push!(CooFormaty, K[2])
        push!(CooFormatval, V)
    end
    
    # boundary modifications
    b = zeros(m + 1)
    b[1] =  - c(h/2)/h^2
    
    # sparse matrix. 
return sparse(CooFormatx, CooFormaty, CooFormatval), b end


"""
    Solve the system: (c(x)u'(x))' = f(x), with u(0) = 1, u'(1) = 0 
    for c(x) = (1 + x)^2. 
    
    Parameters m: 
        The number of interior points for the 1d Grid points, It will inlcude the boundary into the results as well. 
        Therefore, the resutls returned would be m + 2 many points. 
"""
function SolveFor(f::Function, m::Int64)
    c(x) = (1 + x^2)
    h = 1/(m + 1)
    A, b = FiniteDiffMatrix(m, c)
    y = f.(LinRange(0, 1, m + 2)[2:end])
    RHS = y + b
    #display(A)
    # display(RHS)
    u = A\RHS
return LinRange(0, 1, m + 2) |> collect, vcat(1, u) end 

# Coarse Grid Accelerated by Richardson, Running and plotting and stuff
using SparseArrays, LinearAlgebra, Plots

function Problem3Part1()
    f(x) = 2(3x^2 - 2x + 1)   # RHS function. 
    uhat = (x) -> (1 - x)^2
    # h = 0.1 with coarser grid 
    m = 9
    hCoarse = 0.1
    tCoarse, uCoarse = SolveFor(f, m)
    m = 19
    hFine = 0.05
    tFiner, uFiner = SolveFor(f, m)
    richardExtrapolate = (4*uFiner[1:2:end] - uCoarse)/3
    # Plot the results out
    fig = plot(tCoarse, uCoarse - uhat.(tCoarse))
    plot!(fig, tCoarse, uFiner[1:2:end] - uhat.(tCoarse))
    plot!(fig, tCoarse, richardExtrapolate - uhat.(tCoarse))
    display(fig)

    # print the results out. 
    println("Measuring the L2 Norm over the corase grid points for these 3 methods. ")
    print("L2 Norm of the error using only the Coarse Grid is: ")
    sqrt(hCoarse)*norm(uCoarse - uhat.(tCoarse)) |> println
    print("L2 Norm of the error using the fine grid and measured over the coarse grid is: ")
    sqrt(hCoarse)*norm(uFiner[1:2:end] - uhat.(tCoarse)) |> println
    print("L2 Norm of the Error for Richardson on the Coarse Grid")
    sqrt(hCoarse)*norm(richardExtrapolate - uhat.(tCoarse)) |> println
    
return end

Problem3Part1()

