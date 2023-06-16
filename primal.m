% Script to compute the primal SDP formulated in the paper,
% for the Boolean function of id 5865
% of truth table: "0001011011101001" 
% and polynomial representation: x1 + x3x4 + x2x3 + x2x4 + x2x3x4

bits = 4; % Number of bits of the Boolean function considered
dim_H = bits + 1; % The oracles takes a register of dimension 5
T=2;
settings = sdpsettings('showprogress',1,'savesolverinput',1,'savesolveroutput',1,'dualize',0,'solver','scs','scs.eps',1e-6,'scs.eps_abs',1e-6,'scs.eps_rel', 0, 'scs.max_iters',50000,'dimacs',1);

f = boolean([0 0 0 1 0 1 1 0 1 1 1 0 1 0 0 1]); % Loads the function

% Generates the process matrix variable and the constraints for both the FO
% supermap and the general supermap
W = gen_variables(dim_H, T, 0);
constr_QCFO = constraints(W, T, 2);
constr_GEN = constraints(W, T, 5);

% Generate the constraints that all x must be computed with an error
% smaller than epsilon.
oracles = oracles_map(dim_H, bits, T);
epsilon = sdpvar(1,1);

constr = [];
for x = dec2bin(0:2^bits-1)' - '0'
    im = f(x');
    Ox = oracles(num2str(x'));
    constr = [constr, trace(W{im+1}*transpose(Ox)) >= 1 - epsilon];
end
constr_FO = [constr_QCFO, constr];
constr_GEN = [constr_GEN, constr];

% Optimisation
optout_gen = optimize(constr_GEN, -(1-epsilon), settings); % The optimisation is by default a minimisation.
gen_primal = value(epsilon); 

optout_fo = optimize(constr_FO, -(1-epsilon), settings);
fo_primal = value(epsilon);

[fo_primal, gen_primal]



