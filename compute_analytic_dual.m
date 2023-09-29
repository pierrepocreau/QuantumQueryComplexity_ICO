% Script to compute the dual SDP formulated in the paper,
% for the Boolean function of id 5865
% of truth table: "0001011011101001" 
% and polynomial representation: x1 + x3x4 + x2x3 + x2x4 + x2x3x4

bits = 4; % Number of bits of the Boolean function considered
dim_H = bits + 1; % The oracles take a register of dimension 5
T=2;
settings = sdpsettings('showprogress',1,'savesolverinput',1,'savesolveroutput',1,'dualize',0,'solver','scs','scs.eps',1e-6,'scs.eps_abs',1e-6,'scs.eps_rel', 0, 'scs.max_iters',50000,'dimacs',1);

f = boolean_function_from_table([0 0 0 1 0 1 1 0 1 1 1 0 1 0 0 1]); % Loads the function

% Generate the dual variables
dim = dim_H^(2*T);
W = sdpvar(dim,dim,'symmetric');
lambda{1} = sdpvar(8,1);
lambda{2} = sdpvar(8,1);

% Generate the dual cone contraints and the constraints on the lambdas.
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

constr_QCFO_dual = in_QCFO_dual_cone(W, d, A);
constr_GEN_dual = in_valid_dual_cone(W, d, A); 
constr_lambda = [sum(lambda{1}) + sum(lambda{2}) <= 1, lambda{1} >= 0, lambda{2} >= 0];


% Generate the constraints linked to the oracles
oracles = oracles_map(bits, T);

constr_0 = 0;
constr_1 = 0;
c0 = 1;
c1 = 1;
for x = dec2bin(0:2^bits-1)' - '0'
    im = f(x');
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
constr_QCFO_dual = [constr_QCFO_dual, constr_lambda, constr1, constr0];

% Dual objective
obj = 1 + trace(W)/dim_H^2 - sum(lambda{1}) - sum(lambda{2}); 

% Optimisation
optout_gen = optimize(constr_GEN_dual, obj, settings);
W_gen_dual = value(W);
lambda_gen{1} = value(lambda{1});
lambda_gen{2} = value(lambda{2});
% Extract the value of epsilon from the objective (1 - epsilon)
gen_dual = - (trace(value(W_gen_dual))/dim_H^2 - sum(lambda_gen{1}) - sum(lambda_gen{2}));

optout_fo = optimize(constr_QCFO_dual, obj, settings);
W_fo_dual = value(W);
lambda_fo{1} = value(lambda{1});
lambda_fo{2} = value(lambda{2});
% Extract the value of epsilon from the objective (1 - epsilon)
fo_dual = - (trace(value(W_fo_dual))/dim_H^2 - sum(lambda_fo{1}) - sum(lambda_fo{2}));

[fo_dual, gen_dual]