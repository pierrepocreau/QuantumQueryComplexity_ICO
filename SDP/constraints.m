function constr = constraints(W,T,supermapClass)
    % If we use symmetries
    if isa(W{1}, 'replab.CommutantVar')
        W{1} = W{1}.fullMatrix();
        W{2} = W{2}.fullMatrix();
    end

%------------------------
%Causal order constraints
%------------------------
%Tools for constraints calculation
k = T;
dim = size(W{1},1);
dim_H = exp(log(dim)/(2*T));

%   Parties
% {{[1]}, {[2, 3]}, {[4, 5]}, {[6]}}

d = dim_H * ones(1,2*k);
d = [1 d 1];
assert(dim == prod(d));
A{1}{1} = 1;
A{1}{2} = [];
for i = 1:k
    A{i+1}{1} = 2*i;
    A{i+1}{2} = 2*i + 1;
end
A{k+2}{1} =  2*k+2;
A{k+2}{2} = [];

constr = [];
%constr = [W{1} >= 0, W{2} >= 0]; %Might be more efficient when we use the replabCommutant var.

switch supermapClass
    case 1 % Parallel
        constr = [constr, is_QCPAR(W,d, A)];
    case 2 % QC-FO
        constr = [constr, is_QCFO(W,d, A)];
    case 3 % QC-CC
        constr = [constr, is_QCCC(W,d, A)];
    case 4 % QC-QC
        constr = [constr, is_QCQC(W,d, A)];
    otherwise % General supermaps
        constr = [constr, is_valid_superop(W,d, A)];
end
end