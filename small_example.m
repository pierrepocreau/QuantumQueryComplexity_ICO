%%
bits = 2;
dim_H = bits + 1; % We don't consider qubits, but more general system. Furthermore, |0...0> will not query the Oracle.
T = 1;
symmetries = 0;
symbolic = false;

func = boolean([0 0 1 0]);
oracles = oracles_map(dim_H, bits, T);
%% PRIMAL
settings = sdpsettings('showprogress',1,'savesolverinput',1,'savesolveroutput',1,'dualize',0,'solver','scs','scs.eps',1e-6,'scs.eps_abs',1e-6,'scs.eps_rel', 0, 'scs.max_iters',50000,'dimacs',1);

W = gen_variables(dim_H, T, symmetries);
constr_GEN = constraints(W, T, 5);

%Constraint and objective to maximise the probability of sucess for the
%worst possible input x.
oracles = oracles_map(dim_H, bits, T);
epsilon = sdpvar(1,1);
obj = epsilon;

if isa(W{1}, 'replab.CommutantVar')
    W{1} = W{1}.fullMatrix();
    W{2} = W{2}.fullMatrix();
end

constr = [];
for x = dec2bin(0:2^bits-1)' - '0'
    im = func(x');
    Ox = oracles(num2str(x'));
    constr = [constr, trace(W{im+1}*transpose(Ox)) >= epsilon];
end
constr_GEN = [constr_GEN, constr];

%Optimisation
optout_gen = optimize(constr_GEN, -obj, settings);
gen_primal = value(obj);
W_Gen{1} = value(W{1});
W_Gen{2} = value(W{2});
%% DUAL
%Generate variables
dim = dim_H^(2*T);
W = sdpvar(dim,dim,'symmetric');
lambda{1} = sdpvar(3,1);
lambda{2} = sdpvar(1,1);

%Process constraints
dim = size(W,1);
dim_H = exp(log(dim)/(2*T));

d = dim_H * ones(1,2*T);
d = [1 d 1];
A{1}{1} = 1;
A{1}{2} = [];
for i = 1:T
    A{i+1}{1} = 2*i;
    A{i+1}{2} = 2*i + 1;
end
A{T+2}{1} =  2*T+2;
A{T+2}{2} = [];

constr_GEN_dual = in_valid_dual_cone(W, d, A); 
constr_lambda = [sum(lambda{1}) + sum(lambda{2}) >= 1, lambda{1} >= 0, lambda{2} >= 0];

oracles = oracles_map(dim_H, bits, T);

constr_0 = 0;
constr_1 = 0;
c0 = 1;
c1 = 1;
for x = dec2bin(0:2^bits-1)' - '0'
    im = func(x');
    Ox = oracles(num2str(x'));
    if im == 0
        constr_0 = constr_0 + lambda{im+1}(c0) * Ox;
        c0 = c0 + 1;
    else
        constr_1 = constr_1 + lambda{im+1}(c1) * Ox; 
        c1 = c1 + 1;
    end
end
constr0 =  W - constr_0 >= 0;
constr1 =  W - constr_1 >= 0;

constr_GEN_dual = [constr_GEN_dual, constr_lambda, constr1, constr0];

obj = trace(W)/dim_H^T;

optout_gen = optimize(constr_GEN_dual, obj, settings);
gen_dual = value(obj);
W_gen_dual = value(W);
lambda_gen{1} = value(lambda{1});
lambda_gen{2} = value(lambda{2});
%% Lower bound
W_gen_lower = trunc_lower(W_Gen, T, 5, symbolic);
p_gen = [];
for x = dec2bin(0:2^bits-1)' - '0'
    im = func(x');
    Ox = oracles(num2str(x'));
    p_gen = [p_gen, trace(W_gen_lower{im+1}*transpose(Ox))];
end
%%
[W_gen_dual_trunc, lambda_gen_trunc] = trunc_upper(W_gen_dual, lambda_gen, bits, func, T, 5, symbolic);
bounds_gen = [min(p_gen), trace(W_gen_dual_trunc)/dim_H^T];
