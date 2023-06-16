%%
bits = 4;
dim_H = bits + 1; % We don't consider qubits, but more general system. Furthermore, |0...0> will not query the Oracle.
T = 2;
symmetries = 0;
symbolic = false;

%Function id 5865: "0001011011101001" x1 + x3x4 + x2x3 + x2x4 + x2x3x4
%(simplest NPN equivalent function with largest gap)
func = boolean([0 0 0 1 0 1 1 0 1 1 1 0 1 0 0 1]);
%func = boolean([0 1 1 0 0 1 1 0 1 0 0 1 1 1 0 0]), %other f  with gap
oracles = oracles_map(dim_H, bits, T);
%%
load("primal_5865.mat")
load("dual_5865.mat")
%load("primal_f.mat")
%load("dual_f.mat")


%%
% Truncation of primal lower bound
p_fo = [];
p_gen = [];
W_gen_trunc = trunc_lower(W_Gen, T, 5, symbolic);
W_fo_trunc = trunc_lower(W_FO, T, 2, symbolic);
for x = dec2bin(0:2^bits-1)' - '0'
    im = func(x');
    Ox = oracles(num2str(x'));
    p_fo = [p_fo, trace(W_fo_trunc{im+1}*transpose(Ox))];
    p_gen = [p_gen, trace(W_gen_trunc{im+1}*transpose(Ox))];
end

%%
[W_gen_dual_trunc, lambda_gen_trunc] = trunc_upper(W_gen_dual, lambda_gen, bits, func, T, 5, symbolic);
[W_fo_dual_trunc, lambda_fo_trunc] = trunc_upper(W_fo_dual, lambda_fo, bits, func, T, 2, symbolic);

%%
%Bound on epsilon on both cases:
primal_gen = 1-min(p_gen);
dual_gen = - (trace(W_gen_dual_trunc)/25 - sum(lambda_gen{1}) - sum(lambda_gen{2}));

primal_fo = 1-min(p_fo);
dual_fo = - (trace(W_fo_dual)/25 - sum(lambda_fo{1}) - sum(lambda_fo{2}));

%Bounds
bounds_gen = [dual_gen, primal_gen]
bounds_fo = [dual_fo, primal_fo]
