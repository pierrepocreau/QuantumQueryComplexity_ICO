bits = 3;
dim_H = bits + 1; % We don't consider qubits, but more general system. Furthermore, |0...0> will not query the Oracle.
T=2;
symmetries=0;

% x1x2x4 + x1 + x3 + x4
%x4 x3 x2 x1 = 0000, 0001, 0010, 0011, 0100, ...
func = boolean([0 1 0 1 1 0 1 0 1 0 1 1 0 1 0 0]);
ambainis = boolean([0 1 1 1 0 1 0 0 0 0 1 0 1 1 1 0]);
amb2 = boolean([0 1 1 1 1 1 1 0]);
f=amb2;
%Pour les and et equality, on arrive bien au bornes montrées par Ambainis
%dans Optimal one-shot quantum algo...

%Pour and3 bits avec T=2, on a bien le même résultat que Montanaro "On exact quantum
%query complexity

settings = sdpsettings('showprogress',1,'savesolverinput',1,'savesolveroutput',1,'dualize',0,'solver','scs','scs.eps',1e-6,'scs.eps_abs',1e-6,'scs.eps_rel', 0, 'scs.max_iters',50000,'dimacs',1);

dim = 4*dim_H^(2*T);
W = sdpvar(dim,dim,'symmetric');
supermapClass = 2;

%------------------------
%Causal order constraints
%------------------------

%   Parties
% {{[1]}, {[2, 3]}, {[4, 5]}, {[6]}}

d = dim_H * ones(1,2*T);
d = [2 d 2];
assert(dim == prod(d));
A{1}{1} = 1;
A{1}{2} = [];
for i = 1:T
    A{i+1}{1} = 2*i;
    A{i+1}{2} = 2*i + 1;
end
A{T+2}{1} =  2*T+2;
A{T+2}{2} = [];

constr = [];

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

%Constraint and objective to maximise the probability of sucess for the
%worst possible input x.
oracles = oracles_map(dim_H, bits, T);
epsilon = sdpvar(1,1);
obj = epsilon;

for x = dec2bin(0:2^bits-1)' - '0'
    im = f(x');
    Ox = oracles(num2str(x'));
    if im == 0
        constr = [constr, trace(PartialTrace(Tensor([1 0; 0 0], eye(size(Ox, 1)), eye(2))*W*Tensor(eye(2),transpose(Ox), eye(2)),1,[dim/2, 2]) * [1 0; 0 0]) >= 0.999];
    else
        constr = [constr, trace(PartialTrace(Tensor([1 0; 0 0], eye(size(Ox, 1)), eye(2))*W*Tensor(eye(2),transpose(Ox), eye(2)),1,[dim/2, 2]) * [1 0; 0 0]) >= 0.999];
    end
    if im == 0
        constr = [constr, trace(PartialTrace(Tensor([0 0; 0 1], eye(size(Ox, 1)), eye(2))*W*Tensor(eye(2),transpose(Ox), eye(2)),1,[dim/2, 2]) * [1 0; 0 0]) >= epsilon];
    else
        constr = [constr, trace(PartialTrace(Tensor([0 0; 0 1], eye(size(Ox, 1)), eye(2))*W*Tensor(eye(2),transpose(Ox), eye(2)),1,[dim/2, 2]) * [0 0; 0 1]) >= epsilon];
    end
end


%Optimisation
optout_gen = optimize(constr, -obj, settings);
obj = value(obj);
W_ = value(W);

%%

for x = dec2bin(0:2^bits-1)' - '0'
    im = f(x');
    Ox = oracles(num2str(x'));
    output = PartialTrace(W_*Tensor(transpose(Ox), eye(2)),1,[dim/2, 2])*([1; 1]);
    output_bis = Tensor(output, [0; 1]);
    CZ = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1];
    X = [0 1; 1 0];
    output_bis = Tensor(eye(2), X) * CZ * output_bis
end