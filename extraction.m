% Script to extract rigourous lower and upper bounds on epsilon from the dual and
% primal solutions of the SDP. 
 
% Adaptation of a method developed in:
% J. Bavaresco, M. Murao, M. T. Quintino "Strict Hierarchy between Parallel, Sequential, and
% Indefinite-Causal-Order Strategies for Channel Discrimination", arXiv:2011.08300

% This script computes bounds for the following function:
% id 5865
% truth table: "0001011011101001" 
% polynomial representation: x1 + x3x4 + x2x3 + x2x4 + x2x3x4

n = 4; % Number of input bits of the Boolean function considered
T = 2; % Number of queries
dim_H = n + 1;
dO = dim_H^T; % Norm of the process matrix

symbolic = true; % To avoid numerical imprecision we use symbolic representation.
func = boolean_function_from_table([0 0 0 1 0 1 1 0 1 1 1 0 1 0 0 1]); % Truth table for this function
oracles = oracles_map(n, T); % Get all the query oracles Ox for each bit string x

% Loads the numerical solutions obtained from solving the primal and dual SDPs
load("primal_5865.mat")
load("dual_5865.mat")

% Extraction of exact solutions of the primal for general supermaps and FO-supermaps. 
% This gives a lower bound on the objective function 1-epsilon, and hence an upper bound on epsilon.
p_fo = [];
p_gen = [];
W_gen_primal_extract = extract_primal(W_Gen, T, n, 5, symbolic); % Supermap class 5 is general supermaps
W_fo_primal_extract = extract_primal(W_FO, T, n, 2, symbolic); % Supermap class 2 is FO-supermaps

% Compute the analytic value of the lower bound.
for x = dec2bin(0:2^n-1)' - '0'
    im = func(x');
    Ox = oracles(num2str(x'));
    p_fo = [p_fo, trace(W_fo_primal_extract{im+1}*transpose(Ox))];
    p_gen = [p_gen, trace(W_gen_primal_extract{im+1}*transpose(Ox))];
end

% Extraction of exact solution of the dual for general supermaps and FO-supermaps. 
% It gives an upper bound on the objective function 1-epsilon, and hence a lower bound on epsilon.
[W_gen_dual_extract, lambda_gen_extract] = extract_dual(W_gen_dual, lambda_gen, n, func, T, 5, symbolic); % Supermap class 5 is general supermaps
[W_fo_dual_extract, lambda_fo_extract] = extract_dual(W_fo_dual, lambda_fo, n, func, T, 2, symbolic); % Supermap class 2 is FO-supermaps

% Compute the values of the two bounds on epsilon.
primal_gen = 1-min(p_gen);
dual_gen = - (trace(W_gen_dual_extract)/dO - sum(lambda_gen_extract{1}) - sum(lambda_gen_extract{2}));

primal_fo = 1-min(p_fo);
dual_fo = - (trace(W_fo_dual_extract)/dO - sum(lambda_fo_extract{1}) - sum(lambda_fo_extract{2}));

bounds_gen = [dual_gen, primal_gen]
bounds_fo = [dual_fo, primal_fo]
