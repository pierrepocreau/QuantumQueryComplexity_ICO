bits = 4;
dim_H = bits + 1; % We don't consider qubits, but more general system. Furthermore, |0...0> will not query the Oracle.
T=2;
symmetries=0;

%Function id 6630: "0001100111100110": x3x4 + x2 + x2x4 + x2x3 + x2x3x4 + x1
%Should be in the same NPN class as: 
%x1x2x4 + x1 + x3 + x4 : 0110011010011100 (the initial function we looked at)
%f = boolean([0 1 1 0 0 1 1 0 1 0 0 1 1 1 0 0]);

%Function id 5865: "0001011011101001" x1 + x3x4 + x2x3 + x2x4 + x2x3x4
%(simplest NPN equivalent function with largest gap)
f = boolean([0 0 0 1 0 1 1 0 1 1 1 0 1 0 0 1]);


dim = 4*dim_H^(2*T); %times 4 for input and ouput dimensions
W = sdpvar(dim,dim,'symmetric');

%   Parties
% {{[1]}, {[2, 3]}, {[4, 5]}, {[6]}}
k = T;
d = dim_H * ones(1,2*k);
d = [2 d 2];
assert(dim == prod(d));
A{1}{1} = 1;
A{1}{2} = [];
for i = 1:k
    A{i+1}{1} = 2*i;
    A{i+1}{2} = 2*i + 1;
end
A{k+2}{1} =  2*k+2;
A{k+2}{2} = [];

constrGen = is_valid_superop(W,d, A);
%%

%Constraint and objective to maximise the probability of sucess for the
%worst possible input x.
oracles = oracles_map(dim_H, bits, T);
epsilon = sdpvar(1,1);

constr = [];
for x = dec2bin(0:2^bits-1)' - '0'
    im = f(x');
    Ox = oracles(num2str(x'));
    mes = [1 0; 0 0] * im + (1-im)*[0 0; 0 1];
    output0 = PartialTrace(W * kron(kron([1 0; 0 0], transpose(Ox)), eye(2)), 1, [dim/2, 2]);
    output1 = PartialTrace(W * kron(kron([1 0; 0 0], transpose(Ox)), eye(2)), 1, [dim/2, 2]);   
    constr = [constr, trace(output0*[1 0; 0 0]) >= 1]; % Constraint that when input is 0 we do nothing
    constr = [constr, trace(output1*mes) >= 1 - epsilon]; % When input is 1 we want to do the computation
    constr = [constr, output1(2, 1)^2 - output1(1, 1) + output1(1, 1)^2]; % Pure constraint, but it is non-convex ! take (a,b) and (c,d) which verifies the constraint and then look at the middle point...
end
constrGEN = [constrGen, constr];
%%
%Optimisation
settings = sdpsettings('showprogress',1,'savesolverinput',1,'savesolveroutput',1,'dualize',0,'solver','Gurobi');

optout_gen = optimize(constrGEN, -(1-epsilon), settings) %it's a minimisation
gen_primal = value(epsilon); 
W_Gen = value(W);




