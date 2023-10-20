function oracles = oracles_map(n_bits, T)
    % Generate a dictionary giving, for each bitstring x, the
    % Choi matrix of the corresponding query oracle O_x.
    % The query oracle is tensored with itself T times, corresponding to 
    % the number of queries used in the computation.
    
    oracles = containers.Map;
    inputs = dec2bin(0:2^n_bits-1) - '0'; % Generate all possible inputs x
    
    for x = inputs'
        Ox = pure_CJ(query_oracle_x(x));
        Ox = Ox*Ox';
    
        Ox = Tensor(Ox,T); % Tensor T times with itself
    
        oracles(num2str(x')) = Ox;
    end
end

function oracle = query_oracle_x(x)
% Generate the oracle for a bit string x
    n = length(x); % Number of bits
    oracle = zeros(n);
    oracle(1, 1) = 1;
    for i = 1:n
        oracle(i+1, i+1) = (-1)^x(i);
    end
end