function [S_final,lambdas_trunc] = trunc_dual(S, lambdas, bits, func, T, supermapClass, symbolic)
% Truncation of the dual solution. It gives an upper bound on the objective
% (1-epsilon) and therefore a lower bound on epsilon.

    dim = size(S, 1);
    dim_H = exp(log(dim)/(2*T));

    d = dim_H * ones(1,2*T);
    d = [1 d 1];
    A{1}{1} = 1;
    A{1}{2} = [];
    for i = 1:T
        A{i+1}{1} = 2*i;
        A{i+1}{2} = 2*i + 1;
    end
    A{T+2}{1} =  2*T+2;
    A{T+2}{2} = [];
    
    %Clean the lambdas
    lambdas_trunc = cell(1,2);
    [N, D] = rat(lambdas{1});
    lambdas_trunc{1} = N./D;
    [N, D] = rat(lambdas{2});
    lambdas_trunc{2} = N./D;
    
    %Use a symbolic representation.
    if symbolic == true
        lambdas_trunc{1} = sym(lambdas_trunc{1});
        lambdas_trunc{2} = sym(lambdas_trunc{2});
    end
    
    %The lambdas should be lower than one, we distributed the gap if needed.
    gap = 1 - sum(lambdas_trunc{1}) - sum(lambdas_trunc{2});
    if gap < 0
        lambdas_trunc{1} = lambdas_trunc{1} + gap/16;
        lambdas_trunc{2} = lambdas_trunc{2} + gap/16;
    end

    %Cleaning of the dual S
    [N, D] = rat(S);
    S_truncated = N./D;
    
    if symbolic == true
        S_truncated = sym(S_truncated);
    end
    
    %The dual should be self-adjoint.
    S_autoAdjoint = (S_truncated + S_truncated') / 2;

    %Project onto the dual cone.
    switch supermapClass
        case 2 % QC-FO
            S_valid = project_onto_dual_QCFO_superops(S_autoAdjoint, d, A);
        otherwise % General supermaps
            S_valid = project_onto_dual_valid_superops(S_autoAdjoint, d, A);
    end
    
    %Compute the expression of the operator0 and operator1
    oracles = oracles_map(dim_H, bits, T);
    operator0 = 0;
    operator1 = 0;
    c0 = 1;
    c1 = 1;
    for x = dec2bin(0:2^bits-1)' - '0'
        im = func(x');
        Ox = oracles(num2str(x'));
        if im == 0
            operator0 = operator0 + lambdas_trunc{im+1}(c0) * Ox;
            c0 = c0 + 1;
        else
            operator1 = operator1 + lambdas_trunc{im+1}(c1) * Ox;
            c1 = c1 + 1;
        end
    end

   if symbolic == true
        operator0 = sym(operator0);
        operator1 = sym(operator1);
   end
   
   %We first make S_pos positive
   if symbolic == true
       S_double = double(S_valid);
       vp_min = min(eig(S_double));
       eps = vp_min/(vp_min - 1) + 10^-8;
   else
       vp_min = min(eig(S_valid));
       eps = vp_min/(vp_min - 1) + 10^-8;
   end
   S_pos = (1-eps)*S_valid + eps*eye(dim);

   %Make sure that the conditions W - operator0 >= 0 and W-operator1 >=
   %are respected. There could be a better way of doing this step.
   dS_pos = double(S_pos);
   coef = 0;
   step = 2*10^-5;
   
   while min(eig(dS_pos - operator0)) < 0
       dS_pos = dS_pos + step*operator0;
       coef = coef + step;
   end
   S_pos = S_pos + coef * operator0;

   dS_pos = double(S_pos);
   coef = 0;
   while min(eig(dS_pos - operator1)) < 0
       dS_pos = dS_pos + step*operator1;
       coef = coef + step;
   end
   S_final = S_pos + coef * operator1;

   %Check that all the constraints are verified
   [~, flagPos] = chol(S_final);
   [~, flagConstr0] = chol(S_final - operator0);
   [~, flagConstr1] = chol(S_final - operator1);
   assert(flagPos == 0);
   assert(flagConstr0 == 0);
   assert(flagConstr1 == 0);    

    assert(sum(lambdas_trunc{1}) + sum(lambdas_trunc{2}) -1 <= 0) % lambdas sum to 1
    switch supermapClass
        case 2 % QC-FO
            assert(0 == max(max(abs(S_final - project_onto_dual_QCFO_superops(S_final, d, A)))))
        otherwise % General supermaps
            assert(0 == max(max(abs(S_final - project_onto_dual_valid_superops(S_final, d, A)))))
    end
end

