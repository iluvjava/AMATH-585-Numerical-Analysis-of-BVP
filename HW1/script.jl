# Problem 1 and Problem 2 for the HW. 

using Plots, Latexify, Logging, DataFrames





"""
    Second order second differential for fxn at x with h. 
"""
function FiniteDiff2nd2ndAt(fxn, h, x)
    return (fxn(x + h) + fxn(x - h) - 2fxn(x))/h^2
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
    latex = latexify(df, env=:tabular)
    println(latex)
end


Problem1()



