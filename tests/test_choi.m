%Compares two different implementations of the choi.
%x4 x3 x2 x1
input = [0 0 1 1];
oracle = query_oracle_x(input);
U = [1  0  0  0  0;
     0 -1  0  0  0;
     0  0 -1  0  0;
     0  0  0  1  0;
     0  0  0  0  1];

assert(isequal(to_choi(oracle, 5), pure_CJ(U)*pure_CJ(U)'));


input = [1 0 0 1];
oracle = query_oracle_x(input);
U = [1  0  0  0  0;
     0 -1  0  0  0;
     0  0  1  0  0;
     0  0  0  1  0;
     0  0  0  0  -1];

assert(isequal(to_choi(oracle, 5), pure_CJ(U)*pure_CJ(U)'));
