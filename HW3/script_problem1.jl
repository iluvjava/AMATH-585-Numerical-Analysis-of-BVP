"""
    A Struct for storing all the relevant paramters for the 
    nont linear pendulum problems
"""
mutable struct PendulumProblem
    # primary paramters 
    m           # M + 1 is the number of grid points including the boundary conditions. 
    alpha       # boundary condition when t = 0. 
    beta        # boundary condition when t = T. 
    # second dary parameters
    h           # time domain division 

    # running parameters
    theta_k     # the kth guess 

    function PendulumProblem(m, alpha, beta, T, theta0)
        this = new()
        this.m = m
        this.alpha = alpha
        this.beta = beta
        this.h = T/(m + 1)
        this.theta_k = theta0
        return this 
    end

end


"""
    Get the jacobian matrix given a trajector over the time interval. 
"""
function EvalJacobiAt(this::PendulumProblem, theta) 
    @assert length(theta) == this.m "The given trajectory must agree with the problem parameters, more precisely, the grid discretization."*
    "We expect shape ($(this.m),) for function: EvalJacobiAt, but we got: $(size(theta)). "
    diagonals = zeros(this.m)
    @. diagonals = -2/this.h^2 + cos(theta)
    subdiagonals = fill(1/this.h^2, this.m)
    return SymTridiagonal(diagonals, subdiagonals)
end

function EvalJacobiAt(this::PendulumProblem)
    return EvalJacobiAt(this, this.theta_k)
end


"""
    Get the value of the non-linear finite diff given the trajactory of the 
    pendulum. 
"""
function EvalAt(this::PendulumProblem, theta)
    h = this.h
    v = zeros(this.m)
    
    # Boundary for t = 0, alpha. 
    v[1] = (this.alpha - 2*theta[1] + theta[2])/h^2 + sin(theta[1])
    
    # Interior 
    for i in 2:this.m - 1
        v[i] = (theta[i - 1] - 2*theta[i] + theta[i + 1])/h^2 + sin(theta[i])
    end
    
    # Boundary for t = T, beta. 
    v[end] = (theta[end - 1] - 2*theta[end] + this.beta)/h^2 + sin(theta[end])
    return v
end



"""
Evalute the pendulum problem at the most recent guesses from the newton's iteration. 
"""
function EvalAt(this::PendulumProblem)
    return EvalAt(this::PendulumProblem, this.theta_k)
end



"""
    Performs one step of newton iterations based on the initial conditions given
    for the construction of the problem instance. 
"""
function (this::PendulumProblem)()
    theta_k = this.theta_k - EvalJacobiAt(this, this.theta_k)\EvalAt(this, this.theta_k)
    this.theta_k = theta_k
    return copy(this.theta_k)
end


"""
    Solve this given all the parameters needed for the Newton's Method for 
    Iterations. 
"""
function SolveSystemFor(m, alpha, beta, T, initial_guess)
    t = collect(LinRange(0, T, m + 2)[2:end - 1])
    theta0 = nothing
    if isa(initial_guess, Function)
        theta0 = initial_guess.(t)
    else
        theta0 = initial_guess
    end
    problem = PendulumProblem(m, alpha, beta, T, theta0)
    trajectories = Vector{Vector}()
    push!(trajectories, theta0)
    
    for i in 1:200
        push!(trajectories, problem())
        if norm(trajectories[end] - trajectories[end - 1], Inf) < 1e-8
            break
        end
    end
    for v in trajectories
        push!(v, beta)
        pushfirst!(v, alpha)
    end
    return vcat(0, t, T), trajectories
end


### ----------------------------------------------------------------------------

using LinearAlgebra, Plots

function Problem1A()
    m = 64
    alpha = 0.7
    beta = 0.7
    T = 2*pi 
    t = collect(LinRange(0, T, m + 2)[2:end - 1])
    theta0 = @. sin(2t) + 0.7
    problem = PendulumProblem(m, alpha, beta, T, theta0)
    trajectories = Vector()
    push!(trajectories, theta0)
    for i in 1:40
        push!(trajectories, problem())
        if norm(trajectories[end] - trajectories[end - 1], Inf) < 1e-8
            break
        end
    end

    fig = plot(t, trajectories[1], label="initial guess", title="Problem 1 (a)")
    plot!(t, trajectories[end], label = "conveged")
    xlabel!(fig, "t")
    ylabel!(fig, "theta")
    savefig(fig, "problem1(a).png")
    return 
end

Problem1A();

function Problem1B()
    m = 127
    T = 10
    alpha = 0
    beta = 0
    f(t) = 2.8sin((2pi*t)/20)
    t, trajectories = SolveSystemFor(m, alpha, beta, T, f)

    # Boopstrap Boundary Conditions. 
    for epsilon in LinRange(0, 0.7, 10)
        t, trajectories = SolveSystemFor(m, epsilon, epsilon, T, trajectories[end][2:end-1])
    end
    
    # Bootstrap time interval.
    for T in LinRange(10, 20, 20)
        t, trajectories = SolveSystemFor(m, 0.7, 0.7, T, trajectories[end][2:end-1])
    end

    fig = plot(t, trajectories[end], ylim=(0, 1.1pi), xlim=(t[1], t[end]), title="converged solution") 
    display(fig)
    xlabel!(fig, "t")
    ylabel!(fig, "theta")
    savefig(fig, "problem1(b).png")
return end

Problem1B()