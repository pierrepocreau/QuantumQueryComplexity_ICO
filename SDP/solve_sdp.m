function [obj,W,optout] = solve_sdp(W, constr, T, func, bits, oracles, settings)
    
    % If we use symmetries
    if isa(W{1}, 'replab.CommutantVar')
        W{1} = W{1}.fullMatrix();
        W{2} = W{2}.fullMatrix();
    end

    %--------------------
    %Oracles computations
    %--------------------
    dim = size(W{1},1); % Total dimension of the process matrix.

    % Computes an operator aggregating the Oracles for the queries x, such
    % that f(x) = 0 (resp. f(x) = 1).

    inst_0 = zeros(dim); 
    inst_1 = zeros(dim); 

    inputs = dec2bin(0:2^bits-1) - '0'; %Generate all possile inputs
    for x = inputs'
        Ox = oracles(num2str(x'));
        if func(x') == 0
            inst_0 = inst_0 + Ox;
        else
            inst_1 = inst_1 + Ox;
        end
    end

    ops(:,:,1) = inst_0;
    ops(:,:,2) = inst_1;

    %-------------------------
    %SDP options and objective
    %-------------------------
    
    % Maximising the probability in worst case.
    epsilon = sdpvar(1,1);

    if isa(W{1}, 'replab.CommutantVar')
        W{1} = W{1}.fullMatrix();
        W{2} = W{2}.fullMatrix();
    end

    for x = dec2bin(0:2^bits-1)' - '0'
        im = func(x');
        Ox = oracles(num2str(x'));
        constr = [constr, trace(W{im+1}*transpose(Ox)) >= epsilon];
    end
    
    % Maximising the probability in the mean case.
    %obj = trace(W{1}*transpose(ops(:,:,1))) + trace(W{2}*transpose(ops(:,:,2)));
    %obj = real(obj)/2^bits; %Maximize the mean probability of computing f(x) for all x.

    obj = epsilon;
    optout = optimize(constr, -obj, settings);
end

