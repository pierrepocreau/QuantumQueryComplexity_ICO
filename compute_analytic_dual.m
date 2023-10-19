% Script to compute the dual SDP formulated in the paper,
% for the Boolean function of id 5865 with truth table: "0001011011101001" 
% and polynomial representation: x1 + x3x4 + x2x3 + x2x4 + x2x3x4

n = 4; % Number of bits of the Boolean function considered
dim_H = n + 1; % The oracles take a register of dimension 5
T = 2;
% Yalmip settings
settings = sdpsettings('solver','scs','scs.eps',1e-6,'scs.eps_abs',1e-6,'scs.eps_rel', 0, 'scs.max_iters',50000,'dimacs',1);
% settings = sdpsettings('solver','mosek','dimacs',1);

f = boolean_function_from_table([0 0 0 1 0 1 1 0 1 1 1 0 1 0 0 1]); % Loads the function

% Generate the dual variables and their constraints
dim = dim_H^(2*T);
W = sdpvar(dim,dim,'symmetric');
lambda{1} = sdpvar(8,1);
lambda{2} = sdpvar(8,1);

% Specifying the scenario: dimensions, and which space belongs to which operation
d = dim_H * ones(1,2*T);

spaces{1}{1} = [];
for i = 1:T
    spaces{i+1}{1} = 2*i - 1;
    spaces{i+1}{2} = 2*i;
end
spaces{T+2}{1} = [];

%% Building the constraints and objective
Wproj_QCFO = project_onto_dual_QCFO_superops(W, d, spaces);
constr_QCFO_dual = W == Wproj_QCFO;

Wproj_gen = project_onto_dual_valid_superops(W, d, spaces);
constr_GEN_dual = W == Wproj_gen;

constr_lambda = [sum(lambda{1}) + sum(lambda{2}) <= 1, lambda{1} >= 0, lambda{2} >= 0];

% Generate the constraints linked to the oracles
oracles = oracles_map(n, T);

operator0 = 0;
operator1 = 0;
num0s = 1;
num1s = 1;
for x = dec2bin(0:2^n-1)' - '0'
    im = f(x');
    Ox = oracles(num2str(x'));
    if im == 0
        operator0 = operator0 + lambda{im+1}(num0s) * Ox;
        num0s = num0s + 1;
    else
        operator1 = operator1 + lambda{im+1}(num1s) * Ox; 
        num1s = num1s + 1;
    end
end
constr0 =  W - operator0 >= 0;
constr1 =  W - operator1 >= 0;

constr_GEN_dual = [constr_GEN_dual, constr_lambda, constr1, constr0];
constr_QCFO_dual = [constr_QCFO_dual, constr_lambda, constr1, constr0];

% Dual objective
obj = 1 + trace(W)/dim_H^T - sum(lambda{1}) - sum(lambda{2}); 

%% Optimisation
optout_gen = optimize(constr_GEN_dual, obj, settings);
W_gen_dual = value(W);
lambda_gen{1} = value(lambda{1});
lambda_gen{2} = value(lambda{2});
% Extract the value of epsilon from the objective (1 - epsilon)
eps_Gen = - (trace(value(W_gen_dual))/dim_H^T - sum(lambda_gen{1}) - sum(lambda_gen{2}));

optout_fo = optimize(constr_QCFO_dual, obj, settings);
W_fo_dual = value(W);
lambda_fo{1} = value(lambda{1});
lambda_fo{2} = value(lambda{2});
% Extract the value of epsilon from the objective (1 - epsilon)
eps_FO = - (trace(value(W_fo_dual))/dim_H^T - sum(lambda_fo{1}) - sum(lambda_fo{2}));

[eps_FO, eps_Gen]