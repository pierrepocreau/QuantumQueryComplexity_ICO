% Script to extract rigorous lower and upper bounds from the dual and
% primal solutions of the SDP. We adapt a method developped in:
% Strict Hierarchy between Parallel, Sequential, and
% Indefinite-Causal-Order Strategies for Channel Discrimination, Bavaresco,
% Murao, Quintino.

% This script computes bounds for the following function:
% id 5865
% truth table: "0001011011101001" 
% polynomial representation: x1 + x3x4 + x2x3 + x2x4 + x2x3x4
bits = 4; % Number of bits of the Boolean function considered
T=2;
symbolic = true; % To avoid numerical imprecision we use the symbolic representation.

func = boolean([0 0 0 1 0 1 1 0 1 1 1 0 1 0 0 1]);
oracles = oracles_map(bits, T);

% Loads the solutions obtained from solving the primal and dual sdp
load("primal_5865.mat")
load("dual_5865.mat")

% Truncation of the primal solutions. It gives a lower bound on the
% objectif 1-epsilon, so an upper bound on epsilon.
p_fo = [];
p_gen = [];
W_gen_trunc = trunc_lower(W_Gen, T, 5, symbolic);
W_fo_trunc = trunc_lower(W_FO, T, 2, symbolic);

% Compute the value of the lower bound.
for x = dec2bin(0:2^bits-1)' - '0'
    im = func(x');
    Ox = oracles(num2str(x'));
    p_fo = [p_fo, trace(W_fo_trunc{im+1}*transpose(Ox))];
    p_gen = [p_gen, trace(W_gen_trunc{im+1}*transpose(Ox))];
end

% Truncation of the dual solution. It gives an upper bound on the objectif
% 1-epsilon, so a lower bound on epsilon.
[W_gen_dual_trunc, lambda_gen_trunc] = trunc_upper(W_gen_dual, lambda_gen, bits, func, T, 5, symbolic);
[W_fo_dual_trunc, lambda_fo_trunc] = trunc_upper(W_fo_dual, lambda_fo, bits, func, T, 2, symbolic);

% Compute the value of the two bounds on epsilon.
primal_gen = 1-min(p_gen);
dual_gen = - (trace(W_gen_dual_trunc)/25 - sum(lambda_gen{1}) - sum(lambda_gen{2}));

primal_fo = 1-min(p_fo);
dual_fo = - (trace(W_fo_dual)/25 - sum(lambda_fo{1}) - sum(lambda_fo{2}));

bounds_gen = [dual_gen, primal_gen]
bounds_fo = [dual_fo, primal_fo]
