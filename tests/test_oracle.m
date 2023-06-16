%x4 x3 x2 x1
input = [0 0 1 1];
oracle = query_oracle_x(input);

%Query oracle takes as input a vector with a one at the position to the bit
%we wish to query. Note that we need a way of not querying the oracle, here
%it's for the register [1 0 0... 0]
assert(isequal(oracle([1 0 0 0 0]),   [1 0 0 0 0]));
assert(isequal(oracle([0 1 0 0 0]), - [0 1 0 0 0]));
assert(isequal(oracle([0 0 1 0 0]), - [0 0 1 0 0]));
assert(isequal(oracle([0 0 0 1 0]),   [0 0 0 1 0]));
assert(isequal(oracle([0 0 0 0 1]),   [0 0 0 0 1]));

input = [1 1 0 1];
oracle = query_oracle_x(input);

assert(isequal(oracle([1 0 0 0 0]),   [1 0 0 0 0]));
assert(isequal(oracle([0 1 0 0 0]), - [0 1 0 0 0]));
assert(isequal(oracle([0 0 1 0 0]),   [0 0 1 0 0]));
assert(isequal(oracle([0 0 0 1 0]), - [0 0 0 1 0]));
assert(isequal(oracle([0 0 0 0 1]), - [0 0 0 0 1]))