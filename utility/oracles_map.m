function oracles = oracles_map(n_bits, T)
% Generate a dictionary mapping for each bit vector x, its
% corresponding query oracle in its Choi form.
% The query oracles are tensored T times, corresponding to the number
% of query allowed for the computation (up to T=3).

oracles = containers.Map;
inputs = dec2bin(0:2^n_bits-1) - '0'; %Generate all possible inputs

for x = inputs'
    Ox = pure_CJ(query_oracle_x(x));
    Ox = Ox*Ox';

    if T == 2
        Ox = kron(Ox, Ox);
    end
    if T == 3
        Ox = kron(Ox, kron(Ox, Ox));
    end

    oracles(num2str(x')) = Ox;
end
end

function oracle = query_oracle_x(x)
% Generate the oracle of a bit string x.
    bits = length(x);
    oracle = zeros(bits);
    oracle(1, 1) = 1;
    for i = 1:bits
        oracle(i+1, i+1) = (-1)^x(i);
    end
end