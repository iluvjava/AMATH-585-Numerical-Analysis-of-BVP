
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


"""
    Unpreconditioned Conjugate Gradient Algorithm. 
"""
function CGS(
    A,
    b;
    x0=nothing, 
    eps=1e-10, 
    maxitr=20,
    eps2=Inf
)
    if isnothing(x0)
        x0 = similar(b)
    end
    r = b - A*x0
    p = r
    v = Vector{typeof(x0)}()
    push!(v, x0)
    Ap = similar(b)
    for _ in 1:maxitr
        Ap .= A*p
        a = dot(r, r)/dot(p, Ap)
        if a < 0
            @warn ("Matrix is not symmetric positive definite for CGS. a=$(a)")
        end
        push!(v, v[end] + a*p)

        b = 1/dot(r, r)
        r = r - a*Ap
        
        if norm(r, Inf) < eps && norm(v[end] - v[end - 1], Inf) < eps2
            println("norm small as: $(norm(r, Inf))")
            break
        end
        b *= dot(r,r)
        p = r + b*p
    end
    
return hcat(v...) end


"""
    Performs the Guass Sediel Iterations. 
"""
function GS(
    A,
    b;
    x0=nothing, 
    eps=1e-10, 
    maxitr=20,
    eps2=Inf
)
    U = triu(A, 1)
    L = tril(A)
    if isnothing(x0)
        x0 = zeros(size(b))
    end
    v = Vector{typeof(x0)}()
    push!(v, x0)
    for _ in 1:maxitr
        push!(v, L\(b - U*v[end]))
        res = norm(b - A*v[end], Inf)
        if norm(v[end] - v[end - 1], Inf) < eps2 && 
            res < eps
            break
        end
    end
return hcat(v...) end

### Using the code above to answer the questiosn we have -----------------------

using LinearAlgebra, Plots

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


function Problem2GS(m=257)
    lbd = 1; rbd = 3
    A, b = MakeFiniteDiffProblem(m, f=f, a=lbd, b=rbd)
    # A |> display
    # b |> display
    GridPointsBd = (LinRange(0, 1, m) |> collect)
    IntGridPoints = GridPointsBd[2:end - 1]
    solns = GS(A, b, x0=1 .+ 2*IntGridPoints, maxitr=20)

    # Plotting the solutions along with the correct solutions. 
    
    fig = plot(
        IntGridPoints, 
        u(IntGridPoints), 
        legend=:topleft, 
        label="u(x)",
        line=(:dot, 2), 
        xlabel="x",
        ylabel="solution values", 
        title="GS solution Convergence"
    )
    for J in 1:5:size(solns, 2)
        plot!(
            fig, GridPointsBd, vcat(lbd, solns[:, J], rbd),
            label="u_$(J-1)(x)"
        )
    end
    fig |> display
    savefig(fig, "p2_gs_solns.png")
    
    fig = plot(
        legend=:topleft, 
        label="u(x)",
        line=(:dot, 2), 
        xlabel="x",
        ylabel="errors", 
        title="GS Error Oscillation Patterns"
    )
    @info "Problem 2 GS Errors: "
    for J in 1:5:size(solns, 2)
        err = vcat(lbd, solns[:, J], rbd) - u(GridPointsBd)
        plot!(
            fig, GridPointsBd, 
            err,
            label="E_$(J-1)(x)"
        )
        err = norm(err, Inf)/sqrt(m - 1)
        println("err: $(err)")
    end
    fig |> display
    savefig(fig, "p2_gs_errors.png")
end

function Problem2CG(m=257)
    lbd = 1; rbd = 3
    A, b = MakeFiniteDiffProblem(m, f=f, a=lbd, b=rbd)
    # A |> display
    # b |> display
    GridPointsBd = (LinRange(0, 1, m) |> collect)
    IntGridPoints = GridPointsBd[2:end - 1]
    solns = CGS(A, b, x0=1 .+ 2*IntGridPoints, maxitr=20)

    # Plotting the solutions along with the correct solutions. 
    
    fig = plot(
        IntGridPoints, 
        u(IntGridPoints), 
        legend=:topleft, 
        label="u(x)",
        line=(:dot, 2), 
        xlabel="x",
        ylabel="solution values", 
        title="CG solution Convergence"
    )
    for J in 1:5:size(solns, 2)
        plot!(
            fig, GridPointsBd, vcat(lbd, solns[:, J], rbd),
            label="u_$(J-1)(x)"
        )
    end
    fig |> display
    savefig(fig, "p2_cg_solns.png")
    
    fig = plot(
        legend=:topleft, 
        label="u(x)",
        line=(:dot, 2), 
        xlabel="x",
        ylabel="errors", 
        title="CG Error Oscillation Patterns"
    )
    @info "Problem 2 CG Errors: "
    for J in 1:5:size(solns, 2)
        err = vcat(lbd, solns[:, J], rbd) - u(GridPointsBd)
        plot!(
            fig, GridPointsBd, 
            err,
            label="E_$(J-1)(x)"
        )
        err = norm(err, Inf)/sqrt(m - 1)
        println("err: $(err)")
    end
    fig |> display
    savefig(fig, "p2_cg_errors.png")
return end

Problem2CG();
Problem2GS();