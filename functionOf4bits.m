bits = 4;
dim_H = bits + 1; % We don't consider qubits, but more general system. Furthermore, |0...0> will not query the Oracle.
T=2;
symmetries=0;

%Pour les and et equality, on arrive bien au bornes montrées par Ambainis
%dans Optimal one-shot quantum algo...

%Pour and3 bits avec T=2, on a bien le même résultat que Montanaro "On exact quantum
%query complexity
settings = sdpsettings('showprogress',1,'savesolverinput',1,'savesolveroutput',1,'dualize',0,'solver','scs','scs.eps',1e-6,'scs.eps_abs',1e-6,'scs.eps_rel', 0, 'scs.max_iters',50000,'dimacs',1);
%%
load('representative_n4.mat');
bin_rep = dec2bin(representative_n4) - '0';
bin_rep = bin_rep(1:3,:);

flags = [];
for id = bin_rep'
    func = boolean([0 id']); % padding as it gives truth table of 15 bits while we need 16 bits for n=4.
    flag = 1; %if flag = 1 then function is pair
    for x = dec2bin(0:2^bits-1)' - '0'
        im = func(x');
        im2 = func((1-x)');
        if im ~= 1-im2
            flag = 0;
        end
    end
    id
    flags = [flag, flags];
end


%%
load('representative_n4.mat');
bin_rep = dec2bin(representative_n4) - '0';
bin_rep = bin_rep(199:end,:);
FO = [];
GEN = [];
for id = bin_rep'
    func = boolean([0 id']); % padding as it gives truth table of 15 bits while we need 16 bits for n=4.

    W = gen_variables(dim_H, T, symmetries);

    constr_QCFO = constraints(W, T, 2);
    constr_GEN = constraints(W, T, 5);

    %Constraint and objective to maximise the probability of sucess for the
    %worst possible input x.
    oracles = oracles_map(dim_H, bits, T);
    epsilon = sdpvar(1,1);
    obj = epsilon;

if isa(W{1}, 'replab.CommutantVar')
    W{1} = W{1}.fullMatrix();
    W{2} = W{2}.fullMatrix();
end
    constr = [];
    for x = dec2bin(0:2^bits-1)' - '0'
        im = func(x');
        Ox = oracles(num2str(x'));
        constr = [constr, trace(W{im+1}*transpose(Ox)) >= epsilon];
    end
    constr_FO = [constr_QCFO, constr];
    constr_GEN = [constr_GEN, constr];

    %Optimisation
    optout_gen = optimize(constr_GEN, -obj, settings);
    gen_primal = value(obj);
    W_Gen{1} = value(W{1});
    W_Gen{2} = value(W{2});

    optout_fo = optimize(constr_FO, -obj, settings);
    fo_primal = value(obj);
    W_FO{1} = value(W{1});
    W_FO{2} = value(W{2});

    %Results
    [fo_primal, gen_primal]
    FO = [FO, fo_primal]
    GEN = [GEN, gen_primal]
end
