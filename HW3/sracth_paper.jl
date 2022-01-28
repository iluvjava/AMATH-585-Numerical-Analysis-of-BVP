### Check the eigenvalues of the matrix for problem 2. 

function P2Make(m)
    h = 1/(m + 1)
    row = Vector{Int64}()
    col = Vector{Int64}()
    val = Vector{Float64}()
    
    # The first row
    push!(row, 1); push!(col, 1); push!(val, 2/h^2 + 1 + h^2)
    push!(row, 1); push!(col, 2); push!(val, -1/h^2)
    
    # the interior points 
    for i in 2: m - 1
        push!(row, i); push!(col, i); push!(val, 2/h^2 + 1 + (i*h)^2)
        push!(row, i); push!(col, i - 1); push!(val, -1/h^2)
        push!(row, i); push!(col, i + 1); push!(val, -1/h^2)
    end
    
    # the last row
    push!(row, m); push!(col, m); push!(val, 2/h^2 + 1 + (m*h)^2)
    push!(row, m); push!(col, m-1); push!(val, -1/h^2)
    
    return sparse(row, col, val)
end


using SparseArrays, LinearAlgebra
# Test it for different values of h. 
for m = 3: 300
    A = Matrix(P2Make(m))
    h = 1/(m + 1)
    upperbound = 4/h^2 + 1 + (m - 1)^2*h^2
    lowerbound = 1 + (2h)^2
    println("eigmin: $(eigmin(A)), lowerbound: $(lowerbound)")
    println("eigmax: $(eigmax(A)), upper bound $(upperbound)")
end

