# Problem 1 and Problem 2 for the HW. 

# To run the script, you need to add these packages using pkg to your Julia. 
using Plots, Latexify, Logging, DataFrames


"""
    Second order second differential for fxn at x with h. 
"""
function FiniteDiff2nd2ndAt(fxn, h, x)
    return (fxn(x + h) + fxn(x - h) - 2fxn(x))/h^2
end



"""
    one layer of richardson extrapolations 
"""
function Richardson(fxn, h, x, i=1)
    if i == 1
        return FiniteDiff2nd2ndAt(fxn, h, x)
    end
    if i == 2
        return (4*Richardson(fxn, h/2, x, 1) - Richardson(fxn, h, x, 1))/3
    end
    if i == 3
        u = Dict()
        for h_ = [h, h/2, h/4, 0, -h/4, -h/2, -h]
            u[x + h_] = fxn(x + h_)
        end
        phi1(h) = h^(-2)*(u[x + h] - 2u[x] + u[x - h])
        phi2(h) = (4phi1(h/2) - phi1(h))/3
        phi3(h) = (16phi2(h/2) - phi2(h))/15
        return phi3(h)
    end
    error("not yet implemented")
end


function Problem1()
    x = pi/6
    Errors = Vector();
    Computed = Vector();
    hs = 10.0 .^ collect(-1:-1:-16)
    for h in hs
        FiniteDiff = FiniteDiff2nd2ndAt(sin, h, x)
        push!(Computed, FiniteDiff)
        push!(Errors, abs(-sin(x) - FiniteDiff))
    end
    @info "The error for problem 1 is: "
    fig = plot(log10.(hs), log10.(Errors))
    savefig(fig, "problem1_log_log.png")
    @info "The table is: "
    df = DataFrame(h=hs, u=Computed, error=Errors)
    latex = latexify(df, env=:mdtable)
    display(latex)
    latex = latexify(df, env=:tabular, latex=false)
    @info "The table in LaTeX is: "
    println(latex)
end


function Problem2()
    Results = Vector()
    push!(Results, "Computed Results")
    Schemes = Vector()
    push!(Schemes, "Schemes")
    Errors = Vector()
    push!(Errors, "Errors")
    h = 0.2
    push!(Results, Richardson(sin, h, x))
    push!(Schemes, "\$\\varphi_1\$ with h = $h")
    push!(Results, Richardson(sin, h/2, x))
    push!(Schemes, "\$\\varphi_1\$ with h = $(h/2)")
    push!(Results, Richardson(sin, h/4, x))
    push!(Schemes, "\$\\varphi_1\$ with h = $(h/4)")
    push!(Results, Richardson(sin, h, x, 2))
    push!(Schemes, "\$\\varphi_2\$ with h = $(h)")
    push!(Results, Richardson(sin, h/2, x, 2))
    push!(Schemes, "\$\\varphi_2\$ with h = $(h/2)")
    push!(Results, Richardson(sin, h, x, 3))
    push!(Schemes, "\$\\varphi_3\$ with h = $(h)")
    append!(Errors, abs.(Results[2:end] .+ .5))
    TableToPrint = hcat(Schemes, Results, Errors) 
    @info "The table for problem 2 is: "
    latexify(TableToPrint, env=:tabular, latex=false) |> print
    
end

Problem1()
Problem2()




