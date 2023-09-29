% This script makes an exhaustive search (on 4-bit Boolean functions) for an advantage of general
% supermaps over fixed order ones.
% 
% The file "representative_n4.mat" contains the truth tables of all the 222 NPN representatives of 4-bit
% Boolean functions.
%
% For each function we compute the minimum error with which it can be computed using T=2 queries, 
% both with fixed-order quantum supermaps and with general supermaps.

load('representative_n4.mat');
n = 4; % Number of bits
dim_H = n + 1; % For 4-bits, the oracles have dimension 5.
T = 2; % Number of queries

settings = sdpsettings('solver','scs','scs.eps_abs',1e-6,'scs.eps_rel', 0, 'scs.max_iters',50000,'dimacs',1);

% Truth tables for all the functions
bin_rep = dec2bin(representative_n4) - '0';

% Stores the minimum errors reached with T queries for each function.
epsilons_FO = [];
epsilons_Gen = [];

for tt = bin_rep'
    func = boolean([0 tt']); % The NPN representatives are specified on truth tables of 15 bits, we have to pad to have 2^4 = 16 bits.
    id = bin2dec(num2str(tt'));
    
    % Generate the process matrix variables
    dim = dim_H^(2*T);
    W{1} = sdpvar(dim,dim,'symmetric');
    W{2} = sdpvar(dim,dim,'symmetric');
    
    % Specifying the scenario: dimensions, and which space belongs to which operation
    d = dim_H * ones(1,2*T);
    spaces = {{[]}, {1, 2}, {3, 4}, {[],}};

    % Constraints for the two types of supermaps
    constr_FO = is_QCFO(W,d,spaces);
    constr_GEN = is_valid_superop(W,d,spaces);
    
    % Generate the constraints that all x must be computed with an error smaller than epsilon.
    oracles = oracles_map(n, T);
    epsilon = sdpvar(1,1);

    constr = [];
    for x = dec2bin(0:2^n-1)' - '0'
        im = func(x');
        Ox = oracles(num2str(x'));
        constr = [constr, trace(W{im+1}*transpose(Ox)) >= 1-epsilon];
    end

    constr_FO = [constr_FO, constr];
    constr_GEN = [constr_GEN, constr];

    % Optimisation
    optout_gen = optimize(constr_GEN, -(1 - epsilon), settings);
    eps_Gen = value(epsilon);

    optout_fo = optimize(constr_FO, -(1 - epsilon), settings);
    eps_FO = value(epsilon);

    % Storing the results
    disp(['Results for function id ', num2str(id), ': (', num2str(eps_FO), ', ', num2str(eps_Gen), ')']);
    epsilons_FO = [epsilons_FO, eps_FO];
    epsilons_Gen = [epsilons_Gen, eps_Gen];
end
