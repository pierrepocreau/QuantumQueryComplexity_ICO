function oracles = oracles_map(dim_H, n_bits, T)
    % Generate a dictionnary mapping for each bit vector x, its
    % corresponding query oracle in its Choi form.
    % The query oracles are tensored T times, corresponding to the number
    % of query allowed for the computation.

    oracles = containers.Map;
    inputs = dec2bin(0:2^n_bits-1) - '0'; %Generate all possile inputs

    for x = inputs'
        Ox = to_choi(query_oracle_x(x), dim_H);
        if T == 2
            Ox = kron(Ox, Ox);
        end
        if T == 3
            Ox = kron(Ox, kron(Ox, Ox));
        end

        oracles(num2str(x')) = Ox;
    end