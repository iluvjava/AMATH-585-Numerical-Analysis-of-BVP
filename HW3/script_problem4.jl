
mutable struct NonlinearProblem
    # primary paramters 
    m           # M + 1 is the number of grid points including the boundary 
                # conditions. 
    epsilon     # Problem parameters.

    a           # Left Boundary
    b           # right boundary
    alpha       # boundary condition when t = a. 
    beta        # boundary condition when t = b. 
    
    # second dary parameters
    h           # time domain division 

    
    # running parameters
    u     # the kth guess 

    function NonlinearProblem(m, a, b, alpha, beta, epsilon, u0)
        this = new()
        this.m = m
        this.epsilon = epsilon
        this.a = a
        this.b = b
        this.alpha = alpha
        this.beta = beta
        this.h = (b - a)/(m + 1)
        this.u = u0
    return this end

end


"""
    Get the jacobian matrix given a trajector over the time interval. 
"""
function EvalJacobiAt(this::NonlinearProblem, u) 
    @assert length(u) == this.m "The given trajectory must agree with"*
    " the problem parameters, more precisely, the grid discretization."*
    "We expect shape ($(this.m),) for function: EvalJacobiAt,"*
    " but we got: $(size(u)). We only expect values at interior grid points."
    m = this.m; h = this.h; eps = this.epsilon
    u0 = this.alpha; uLast = this.beta
    row = Vector{Float64}()
    col = Vector{Float64}()
    values = Vector{Float64}()
    # construct the matrix row by row. 
    # First row
    push!(row, 1);push!(col, 1); push!(values, 
        (-2eps)/h^2 + (u[2] - u0)/(2h) - 1
    )
    push!(row, 1);push!(col, 2); push!(values, 
        eps/h^2 + u[1]/(2h)
    )
    # Interior row in the middle.
    for Jdx in 2: m - 1
        nonDiag = 
        push!(row, Jdx); push!(col, Jdx - 1); push!(
            values, 
            eps/h^2 - u[Jdx]/(2h)
        )
        push!(row, Jdx); push!(col, Jdx + 1); push!(
            values, 
            eps/h^2 + u[Jdx]/(2h)
        )
        push!(row, Jdx); push!(col, Jdx); push!(
            values, 
            (-2eps)/h^2 + (u[Jdx + 1] - u[Jdx - 1])/(2h) - 1
        )
    end
    # Last row 
    push!(row, m); push!(col, m - 1); push!(values, 
        eps/h^2 - u[end]/(2h)
    )
    push!(row, m); push!(col, m); push!(values, 
        (-2eps)/h^2 + (uLast - u[end])/(2h) - 1
    )
    
    return sparse(row, col, values) 
end


function EvalJacobiAt(this::NonlinearProblem)
    return EvalJacobiAt(this, this.u) 
end


"""
    Get the value of the non-linear finite diff given the trajactory of the 
    pendulum. 
"""
function EvalAt(this::NonlinearProblem, u)
    h = this.h
    m = this.m
    epsilon = this.epsilon
    v = zeros(this.m)
    u0 = this.alpha
    uLast = this.beta
    
    # Boundary for t = 0, alpha. 
    v[1] = epsilon*(
            (u0 - 2u[1] + u[2])/h^2
        ) + u[1]*(
            (u[2] - u0)/(2h) - 1
        )
    
    # Interior 
    for i in 2:m - 1
        v[i] = epsilon*(
            (u[i - 1] - 2u[i] + u[i + 1])/h^2
        ) + u[i]*(
            (u[i + 1] - u[i - 1])/(2h) - 1
        )
    end
    
    # Boundary for t = T, beta. 
    v[end] = epsilon*(
        (u[m - 1] - 2u[m] + uLast)/h^2
    ) + u[m]*(
        (uLast - u[m - 1])/(2h) - 1
    )

return v end



"""
Evalute the pendulum problem at the most recent guesses from the newton's 
iteration. 
"""
function EvalAt(this::NonlinearProblem)
    return EvalAt(this::NonlinearProblem, this.u) 
end



"""
    () operator overrided as performing one step of the newton's method 
    using the current solution, which is a running variable. 
"""
function (this::NonlinearProblem)()
    newGuess = this.u - EvalJacobiAt(this)\EvalAt(this)
    this.u = newGuess
    return copy(newGuess)
end


### Using the code we have to solve the system with a given set of parameters. 
using LinearAlgebra, SparseArrays, Plots
pyplot() # change backend

function PerformNewtonIteration(m, a, b, beta, alpha, epsilon)
    # Suggested Initial Guess: 
    w0 = .5(a - b + beta - alpha)
    xbar = .5(a + b - alpha - beta)
    x = LinRange(a, b, m + 2)
    u0 = @. x - xbar + w0*tanh(w0*(x - xbar)/(epsilon))
    # Set up the problem and solve
    P = NonlinearProblem(m, a, b, alpha, beta, epsilon, u0[2:end - 1])
    previousGuess = P.u
    for _ in 1:30
        newGuess = P()
        if norm(previousGuess - newGuess, Inf) < 1e-10
            return x, vcat(alpha, newGuess, beta)
        end
        previousGuess = newGuess
    end
    
return error("Newton's Method Failed to converge. ") end


function Problem4Solve(m=127)
    a = 0
    b = 1
    alpha = -1 
    beta = 1.5
    epsilon = 0.01
    fig = plot(
        title="u'' + u'(u - 1) = 0, u(a) = alpha, u(b) = beta",
        legend=:bottom,
        size=(800, 500)
    )
    for (m, marker) in zip([159, 79, 39, 19], 
            [:vline, :hline, :+, :x]
        )
        x, u = PerformNewtonIteration(m, a, b, beta, alpha, epsilon)
        plot!(
            fig, x, u, 
            label="h=$(1/(m + 1))", 
            markershape=marker, markersize=5
        )
    end
    xlabel!(fig, "x")
    ylabel!(fig, "u")
    display(fig)
    savefig(fig, "problem4_plot.png")
return end

Problem4Solve()
