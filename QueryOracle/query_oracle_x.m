function handle = query_oracle_x(x)
% Creates a handle to query the input x.
% x is a bit vector.
% The handle function takes as input a vector i encoding a quantum
% register. Is the quantum register as n qubits then i is a 2**n bit
% vector.
    handle = @(i) query_x(x, i);
end

function i = query_x(x, i)
% A call to the oracle_x will modify the vector i.
% Depending on the size of x, the vector i might encode for too many index.
% i = [1, 0, ...] encode the quantum register |00...> which doesn't
% query x (this is necessary for reversibility of the query oracle).

%This function should be call with a one-hot encoding of the quantum
%register.
    assert(sum(i) == 1);

    bits = length(x);
    index = find(i, 1, 'first');

    if index ~=1 && index - 1 <= bits
        i = (-1) ^ x(index-1) * i;
    end
end


