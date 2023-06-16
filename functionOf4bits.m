% This script makes an exhaustive search (on 4-bit Boolean functions) for an advantage of general
% supermaps over fixed order ones.
% 
% The file "representative_n4.mat" contains the truth tables of all the 222 NPN representatives of 4-bit
% boolean functions.
%
% For each function, we compute the minimum error with which it can be
% computed using 2 queries, with a fixed order quantum supermap, and with a
% general supermap.

load('representative_n4.mat');
bits = 4;
dim_H = bits + 1; % For 4-bits, the oracle are of dimension 5.
T=2;

settings = sdpsettings('showprogress',1,'savesolverinput',1,'savesolveroutput',1,'dualize',0,'solver','scs','scs.eps',1e-6,'scs.eps_abs',1e-6,'scs.eps_rel', 0, 'scs.max_iters',50000,'dimacs',1);

bin_rep = dec2bin(representative_n4) - '0';
bin_rep = bin_rep(1:3,:); % Selecting a small sample to test.

% Stores the minimum errors reached with two queries for each function.
EFOs = [];
EGens = [];

for tt = bin_rep'
    func = boolean([0 tt']); % The NPN representatives are specified on truth tables of 15 bits, we have to pad to have 2^4 = 16 bits.
    
    % Generate the process matrix of correct dimensions
    W = gen_variables(dim_H, T, 0); % We don't use symmetries.
    
    % Generate the constraints for fixed order and general supermaps.
    constr_QCFO = constraints(W, T, 2);
    constr_GEN = constraints(W, T, 5);
    
    % Generate the constraints that all x must be computed with an error
    % smaller than epsilon.
    oracles = oracles_map(dim_H, bits, T);
    epsilon = sdpvar(1,1);

    constr = [];
    for x = dec2bin(0:2^bits-1)' - '0'
        im = func(x');
        Ox = oracles(num2str(x'));
        constr = [constr, trace(W{im+1}*transpose(Ox)) >= 1-epsilon];
    end

    constr_FO = [constr_QCFO, constr];
    constr_GEN = [constr_GEN, constr];

    % Optimisation
    optout_gen = optimize(constr_GEN, -(1 - epsilon), settings);
    EGen = value(obj);

    optout_fo = optimize(constr_FO, -(1 - epsilon), settings);
    EFO = value(obj);

    % Storing the results
    [EFO, EGen]
    EFOs = [EFOs, EFO];
    EGens = [EGEns, EGen];
end
