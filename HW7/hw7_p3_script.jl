using LinearAlgebra


# ---- 1D Problem  Generation --------------------------------------------------
"""
    Make a finite difference method for the following system: 
        u''(x) = f(x) u(0) = a, u(1) = b. 
    Dirichilet boundary conditions
    Parameter: 
        n: grid point including the boundary point. The number of interval is 
        n - 1, and the interval size is h = 1/(n - 1). 
"""
function MakeFiniteDiffProblem(
    n;
    f::T1=nothing, 
    a::T2=0, 
    b::T2=0
) where {T1<:Union{Function, Nothing}, T2<:Number}
    @assert n >= 3 "The number of grid point including"*
        "the boundary conditions has to be at least 3"
    if isnothing(f)
        f = (x) -> 0
    end
    h = 1/(n - 1)
    gridPoints = LinRange(0, 1, n) |>  collect
    SymM = SymTridiagonal(fill(-2.0, n - 2), fill(1.0, n - 3))
    SymM ./= h^2
    RHS = f.(gridPoints[2:end - 1])
    RHS[1] -= a/h^2
    RHS[end] -= b/h^2 
return SymM, RHS end

function MakeFiniteDiffMatrix(n)
    A, _ = MakeFiniteDiffProblem(n)
return A end

# Iterative Solvers for linear system ------------------------------------------
function WeightedJacobiUpdate(
    A, b,
    x0;
    w=2/3
)
    D = reshape((diagm(A)./w).^(-1), length(x0), 1)
return D.*b + (I - D.*A)*x0 end


"""
    Performs the Guass Sediel Iterations. 
"""
function GSUpdate(
    A,
    b,
    x0
)
    U = triu(A, 1)
    L = tril(A)
return L\(b - U*x0) end

function diagm(T::SymTridiagonal)
return T.dv end


# Multi-Grid related methods ---------------------------------------------------

function UpScaleWith(v)
    vFine = zeros(2*length(v) + 1)
    vFine[2:2:end - 1] = v
    vFine[1] = v[1]/2
    vFine[3:2:end - 2] = 0.5*(v[1:end - 1] + v[2: end])
    vFine[end] = v[end]/2
return vFine end

"""
    Give it an initial guess, it will modify it by performing 
    one step of 2-Grid Multi-Grid method and returns the 
    residual vector. 
"""
function MultiGridOnce!(
    A,
    A_coarse,
    rhs,
    x0; 
    solver::Function = WeightedJacobiUpdate,
    damp_itr=1
)
    for _ in 1:damp_itr
        x0 .= solver(A, rhs, x0)
    end
    r = rhs - A*x0
    
    # rCoarse = r[2:2:end - 1]
    rCoarse = view(r, 2:2:length(r) - 1)
    eCoarse = A_coarse\rCoarse
    eFine = UpScaleWith(eCoarse)
    x0 .+= eFine
    

return rhs - A*x0 end

function f(x)
    phi(x) = 20*pi*x^3
    Dphi(x) = 60*pi*x^2
    DDphi(x) = 129*pi*x
    a = 0.5
return -20 + a*DDphi(x)*cos(phi(x)) - a*(Dphi(x)^2)*sin(phi(x)) end


function f(x::AbstractVecOrMat)
return f.(x) end


function u(x)
    a = 0.5
    phi(x) = 20*pi*x^3
return 1 + 12*x - 10*x^2 + a*sin(phi(x)) end


function u(x::AbstractVecOrMat)
return u.(x) end


### Using the above function to solve problems ---------------------------------

using Plots

function PerformMultiGrid(n = 257, solver=WeightedJacobiUpdate)
    xgridWithBC = LinRange(0, 1, n) |> collect
    lbc = 0
    rbc = 0
    A, b = MakeFiniteDiffProblem(n, f=f, a=lbc, b=rbc)
    Ac = MakeFiniteDiffMatrix(div((n - 1), 2) + 1)
    x0 = zeros(size(b))
    TotalItr = 0
    for Itr in 1:500
        r = MultiGridOnce!(A, Ac, b, x0, solver=solver) 
        relRes = norm(r)/norm(b)
        println(relRes)
        if relRes < 1e-4
            TotalItr = Itr
            println("total Itr: $(Itr)")
            break
        end
    end

return TotalItr end

function MultiGridConvergence()
    WJacobiItr = Vector{Int64}()
    GSItr =  Vector{Int64}()
    GridPointsCount = 2 .^ collect(4:14)
    for n in GridPointsCount
        n = n + 1
        push!(WJacobiItr, PerformMultiGrid(n, WeightedJacobiUpdate))
        push!(GSItr, PerformMultiGrid(n, GSUpdate))
    end
    fig = plot(
        GridPointsCount.|>log2, WJacobiItr, 
        label="Jacobi MG Iter", 
        marker=:x, 
        xlabel="log2(n)", ylabel="Iterations", 
        legend=:topleft
    )
    plot!(
        fig, 
        GridPointsCount.|>log2, GSItr, 
        label="GS MG Iter", 
        marker=:+
    )
    fig |> display
    savefig(fig, "p3_mg_itr.png")

end

MultiGridConvergence();