% Script to solve the primal SDP formulated in the paper,
% for the Boolean function of id 5865 with truth table: "0001011011101001" 
% and polynomial representation: x1 + x3x4 + x2x3 + x2x4 + x2x3x4

n = 4; % Number of bits of the Boolean function considered
dim_H = n + 1; % The oracles take a register of dimension 5
T = 2;
% Yalmip settings
settings = sdpsettings('solver','scs','scs.eps_abs',1e-6,'scs.eps_rel', 0, 'scs.max_iters',50000,'dimacs',1);
% settings = sdpsettings('solver','mosek','dimacs',1);

f = boolean_function_from_table([0 0 0 1 0 1 1 0 1 1 1 0 1 0 0 1]); % Loads the function

% Generate the process matrix variables. Can assume this is real since the oracles are real.
dim = dim_H^(2*T);
W{1} = sdpvar(dim,dim,'symmetric');
W{2} = sdpvar(dim,dim,'symmetric');

% Specifying the scenario: dimensions, and which space belongs to which operation
d = dim_H * ones(1,2*T);

spaces{1}{1} = [];
for i = 1:T
    spaces{i+1}{1} = 2*i - 1;
    spaces{i+1}{2} = 2*i;
end
spaces{T+2}{1} = [];

constr_QCFO = is_QCFO(W,d, spaces);
constr_GEN = is_valid_superop(W,d, spaces);

% Generate the constraints that all x must be computed with an error
% smaller than epsilon.
oracles = oracles_map(n, T);
epsilon = sdpvar(1,1);

constr = [];
for x = dec2bin(0:2^n-1)' - '0'
    im = f(x');
    Ox = oracles(num2str(x'));
    constr = [constr, trace(W{im+1}*transpose(Ox)) >= 1 - epsilon];
end
constr_FO = [constr_QCFO, constr];
constr_GEN = [constr_GEN, constr];

% Optimisation
optout_GEN = optimize(constr_GEN, -(1-epsilon), settings); % The optimisation is by default a minimisation.
eps_GEN = value(epsilon); 

optout_FO = optimize(constr_FO, -(1-epsilon), settings);
eps_FO = value(epsilon);

[eps_GEN, eps_FO]