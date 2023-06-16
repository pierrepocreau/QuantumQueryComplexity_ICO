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

settings = sdpsettings('showprogress',1,'savesolverinput',1,'savesolveroutput',1,'dualize',0,'solver','scs','scs.eps',1e-6,'scs.eps_abs',1e-6,'scs.eps_rel', 0, 'scs.max_iters',50000,'dimacs',1);

W = gen_variables(dim_H, T, symmetries);

constr_QCFO = constraints(W, T, 2);
constr_GEN = constraints(W, T, 5);

%Constraint and objective to maximise the probability of sucess for the
%worst possible input x.
oracles = oracles_map(dim_H, bits, T);
epsilon = sdpvar(1,1);

if isa(W{1}, 'replab.CommutantVar')
    W{1} = W{1}.fullMatrix();
    W{2} = W{2}.fullMatrix();
end

constr = [];
for x = dec2bin(0:2^bits-1)' - '0'
    im = f(x');
    Ox = oracles(num2str(x'));
    constr = [constr, trace(W{im+1}*transpose(Ox)) >= 1 - epsilon];
end
constr_FO = [constr_QCFO, constr];
constr_GEN = [constr_GEN, constr];

%Optimisation
optout_gen = optimize(constr_GEN, -(1-epsilon), settings); %it's a minimisation
gen_primal = value(epsilon); 
W_Gen{1} = value(W{1});
W_Gen{2} = value(W{2});

optout_fo = optimize(constr_FO, -(1-epsilon), settings);
fo_primal = value(epsilon);
W_FO{1} = value(W{1});
W_FO{2} = value(W{2});

%Results
[fo_primal, gen_primal]



