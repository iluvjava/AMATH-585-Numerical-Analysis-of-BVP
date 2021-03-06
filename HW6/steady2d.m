close all; clear all; clc; 

%% Stationary Iterative Methods Convergence. 
[A, b] = MakeTestProblem(8);
[~, errsJB] = StationaryIterative(A, b);
[~, errsGS] = StationaryIterative(A, b, [], "gs");
[~, errsSOR] = StationaryIterative(A, b, [], "sor");
ErrsJBLogged = log10(errsJB);
ErrsGSLogged = log10(errsGS);
ErrsSORLogged = log10(errsSOR);
figure;
plot(ErrsJBLogged,"-+");
hold on;
plot(ErrsGSLogged, "-o");
plot(ErrsSORLogged, "-.");
legend(["JB", "GS","SOR"]); 
xlabel("Iteration");
ylabel("Log10 of Relative Residual");
title("Stationary Iterative Method Convergence");
saveas(gcf, "stationary_methods.png");
% Stationary Iterative Method: Convergence estimate wrt h, the grid size. 

%% CG with and without and with Preconditioning. 
[A, b] = MakeTestProblem(20);
[~, ErrsCG] = PerformCG(A, b);
[~, ErrsPCG] = PerformCG(A, b, 1); 
figure; 
plot(log10(ErrsCG)); 
hold on;
plot(log10(ErrsPCG));
xlabel("Iteration Count");
ylabel("Log10 Relative Residual"); 
legend(["cg", "pcg+ichol"]);
saveas(gcf, "pcg_vs_cgs_itr.png");

%% Convergence rate wrt to h the grid size for all method. 
GridDivisions = 5:30;
ItrJB = [];
ItrGS = [];
ItrSOR = [];
ItrCGS = [];
ItrPCG = [];
for n = GridDivisions
    ItrJB(end + 1) = GetIterationCountForMethod(@StationaryIterative, n);
    ItrGS(end + 1) = GetIterationCountForMethod( ...
        @(A, b) StationaryIterative(A, b, [], "gs"), n ...
        );
    ItrSOR(end + 1) = GetIterationCountForMethod( ...
        @(A, b) StationaryIterative(A, b, [], "sor"), n ...
        );
    ItrCGS(end + 1) = GetIterationCountForMethod(@PerformCG, n);
    ItrPCG(end + 1) = GetIterationCountForMethod( ...
        @(A, b) PerformCG(A, b, 1), n ...
    );
end

%% Plotting it out. 
close all;
figure;
plot(GridDivisions, ItrJB);hold on
plot(GridDivisions, ItrGS);
plot(GridDivisions, ItrSOR);
plot(GridDivisions, ItrCGS);
plot(GridDivisions, ItrPCG);
legend(["JB", "GS", "SOR", "CGS", "PCG"], "location", "northwest");
xlabel("Number of Grids Partition on one Dimension");
ylabel("Iterations of Methods");
title("Iteration Count vs Grid Division")
saveas(gcf, "h_vs_methods_itr.png");

figure; 
loglog(GridDivisions, ItrJB, "-x");hold on
loglog(GridDivisions, ItrGS, "-o");
loglog(GridDivisions, ItrSOR, "-.");
legend(["JB", "GS", "SOR"], "location", "northwest");
xlabel("Log of Iteration count");
ylabel("Log of Numer of Grid Partition on one Dimension");
title("The Stationary Methods"); 
saveas(gcf, "h_vs_stationary_methods.png");


figure; 
loglog(GridDivisions, ItrPCG, '-x'); hold on ;
loglog(GridDivisions, ItrSOR, '-o');
legend(["pcg", "sor"]);
xlabel("Log of Iteration count");
ylabel("Log of Numer of Grid Partition on one Dimension");
title("The PGC and SOR")
saveas(gcf, "h_vs_pcg_sor.png");



%% Numerically Compute it. 
% Iteration vs number of grid division on one dimension
ConvergenceRateSOR = LogLogSlopeEstimate(GridDivisions, ItrSOR);
disp(ConvergenceRateSOR);
ConvergenceRateGS = LogLogSlopeEstimate(GridDivisions, ItrGS);
disp(ConvergenceRateGS);
ConvergenceRatePCG = LogLogSlopeEstimate(GridDivisions, ItrPCG);
disp(ConvergenceRatePCG);


function [soln, RelativeErrs] = PerformCG(A, b, precon)
    if nargin < 3
        precon = 0; 
    end
    tol = 1e-8;
    if precon == 0
        [soln, FLAG, ~, ~, RelativeErrs] = cgs(A, b, tol, length(b)*10);
        disp("The Flag for cgs is: ");
        disp(num2str(FLAG));
    else
        L = ichol(A);
        M = L*L';
        [soln, FLAG, ~, ~, RelativeErrs] = pcg(A, b, tol, length(b)*10, M);
        disp("The Flag for pcg is: ");
        disp(num2str(FLAG));
    end
end

function TotalItr = GetIterationCountForMethod(fxn, n)
    [A, b] = MakeTestProblem(n); 
    [~, Errs] = fxn(A, b);
    TotalItr = length(Errs); 
end

function slope = LogLogSlopeEstimate(x, y)
    slope = mean( ...
            (log(y(1:end - 1)) - log(y(2: end))) ...
            /(log(x(1: end - 1)) - log(x(2: end))) ...
        );
end