function [soln, RelativeResErr] = StationaryIterative(A, b, epsilon, arg3, w, x0, maxitr)
%%% Function performs stational iterations 
%%% INTPUT: 
%%%     A: 
%%%         The sparse matrix for performing the vector operation.
%%%     b: 
%%%         The b vector for the RHS of the system. 
%%%     arg3: 
%%%         The type of method that we are using for the system. By default
%%%         it uses the Jacobi iteration if this is not set. 
%%%
%%% OUTPUT: 
%%%     soln: The solution in the end. 
%%%     RelativeResErr: List of ||Ax - b||/||b|| during the iteration. 
 

    if ~exist("epsilon", "var") || isempty(epsilon)
        epsilon = 1e-8;
    end
    if ~exist("arg3", "var") || isempty(arg3)
        arg3 = "jb";
    end
    if ~exist("w", "var") || isempty(w)
        w = 1.5;
    end
    if ~exist("x0", "var") || isempty(x0)
        x0 = zeros(size(b));
    end
    if ~exist("maxitr", "var") || isempty(maxitr)
        maxitr = max(2*size(b, 1), 1000);
    end

    L = tril(A, -1); U = triu(A, 1); d = diag(A); 
    D = diag(d);
    
    if ~any(["jb", "gs", "sor"] == arg3)
        error("arg3 must be one of the following: jb, gs, sor. ")
    else  
        maxEig = abs(eigs((L + U)./d, 1, 'largestabs'));
        if arg3 == "jb"
            disp("Stationary Iterative Method: jb")
            if maxEig > 1 || isnan(maxEig)
                disp("Jacobi might not converge, max eig of the matrix is:");
                disp(num2str(maxEig));
                disp("Here is the plot of the matrix")
                figure;
                imagesc((L + U)./d); colorbar;
            end
        else
            if arg3 == "gs"
                disp("Staionary Iterative method: GS");
                w = 1;
            else  % sor
                disp("Stationary Iterative method: SOR");
                w = 2/(1 + sqrt(1 - maxEig^2));
                disp("sor determined relaxtion factor is: ");
                disp(num2str(w));
            end
        end
    end
    
    RelativeResErr = zeros(1, maxitr);
    RelativeResErr(1) = norm(b - A*x0)/norm(b);
    
    
    for Itr = 1: maxitr
        if arg3 == "jb"
            x0 = (b - (L + U)*x0)./d;
        else  % gs or sor
            x0 = (D + w*L)\(w*b - (w*U + (w - 1)*D)*x0);
        end
        
        RelativeResErr(Itr + 1) = norm(b - A*x0)/norm(b);
        
        if RelativeResErr(Itr + 1) < epsilon
            break;
        end

        if isnan(RelativeResErr(Itr + 1)) || isinf(RelativeResErr(Itr + 1))
            disp("Error during stationary iteration hit nan or inf."); 
            break;
        end

    end

    RelativeResErr = RelativeResErr(1:Itr + 1);
    soln = x0;

end

