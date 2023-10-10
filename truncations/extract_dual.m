function [S_final,lambdas_frac] = extract_dual(S, lambdas, n, func, T, supermapClass, symbolic)
% Extract an exact dual solution from a numerical approximation. 
% This function gives an upper bound on the objective (1-epsilon) and therefore a lower bound on epsilon.
% It follows Algorithm 1 as described in the paper 
% "Quantum Query Complexity of Boolean Functions under Indefinite Causal Order".

    dim_H = n+1;
    
    d = dim_H * ones(1,2*T);
    spaces{1}{1} = []; % trivial P
    for i = 1:T
        spaces{i+1}{1} = 2*i - 1;
        spaces{i+1}{2} = 2*i;
    end
    spaces{T+2}{1} =  []; % trivial F
    
    % Verify all lambdas are non-negative
    for i = 1:2
        lambdas{i} = max(lambdas{i},zeros(1,length(lambdas{i})));
    end

    % Rationalise the lambdas
    lambdas_frac = cell(1,2);
    [N, D] = rat(lambdas{1});
    lambdas_frac{1} = N./D;
    [N, D] = rat(lambdas{2});
    lambdas_frac{2} = N./D;
    
    % Use a symbolic representation
    if symbolic == true
        lambdas_frac{1} = sym(lambdas_frac{1});
        lambdas_frac{2} = sym(lambdas_frac{2});
    end
    
    % The lambdas should be smaller than one. We distribute the difference delta if needed.
    lambdas_final = cell(1,2);
    delta = 1 - sum(lambdas_frac{1}) - sum(lambdas_frac{2});
    if delta < 0
        % Determine how many lambdas we can adjust while keeping them positive.
        all_lambdas = sort([lambdas{1}, lambdas{2}],'descend');
        N = 2^n;
        while all_lambdas(N) < -delta/N
            N = N - 1;
        end

        % Adjust those lambdas
        for x = 1:length(lambdas_frac{1})
            lambdas1 = lambdas_frac{1};
            lambdas_final1 = zeros(1,length(lambdas1));
            if lambdas1(x) >= -delta/N
                lambdas_final1(x) = lambdas1(x) + delta/N;
            else
                lambdas_final1(x) = lambdas1(x);
            end
            lambdas_final{1} = lambdas_final1;
        end
        for x = 1:length(lambdas_frac{2})
            lambdas2 = lambdas_frac{2};
            lambdas_final2 = zeros(1,length(lambdas2));
            if lambdas2(x) >= -delta/N
                lambdas_final2(x) = lambdas2(x) + delta/N;
            else
                lambdas_final2(x) = lambdas2(x);
            end
            lambdas_final{2} = lambdas_final2;
        end
    else
        lambdas_final = lambdas_frac;
    end
    
    % Rationalise the dual S
    [N, D] = rat(S);
    S_frac = N./D;
    
    if symbolic == true
        S_frac = sym(S_frac);
    end
    
    % Ensure the dual is Hermitian
    S_Herm = (S_frac + S_frac')/2;
    
    % Project S onto the dual cone
    switch supermapClass
        case 2 % QC-FO
            S_valid = project_onto_dual_QCFO_superops(S_Herm, d, spaces);
        otherwise % General supermaps
            S_valid = project_onto_dual_valid_superops(S_Herm, d, spaces);
    end
    
    % Make sure S is positive semidefinite
    if symbolic == true
       S_double = double(S_valid);
       eig_min = min(eig(S_double));
       mu = eig_min/(eig_min - 1) + 10^-8; % add some extra cusion to ensure strict positivity given numeric eig_min
    else
       eig_min = min(eig(S_valid));
       mu = eig_min/(eig_min - 1) + 10^-8;
    end
    S_pos = (1-mu)*S_valid + mu*eye(size(S_valid));
    
    % Compute the expressions of the operator0 and operator1 (O^{[0]} and O^{[1}}) in the paper
    oracles = oracles_map(dim_H, n, T);
    operator0 = 0;
    operator1 = 0;
    num0s = 0;
    num1s = 0;
    for x = dec2bin(0:2^n-1)' - '0'
        im = func(x');
        Ox = oracles(num2str(x'));
        if im == 0
            num0s = num0s + 1;
            operator0 = operator0 + lambdas_final{im+1}(num0s) * Ox;
        else
            num1s = num1s + 1;
            operator1 = operator1 + lambdas_final{im+1}(num1s) * Ox;
        end
    end
    
    if symbolic == true
        operator0 = sym(operator0);
        operator1 = sym(operator1);
    end
    
    % Make sure that the conditions S - operator0 >= 0 and S - operator1 >= 0 are respected. 
    step = 2*10^-5;
    
    dS_pos = double(S_pos);
    eta0 = 0;
    while min(eig(dS_pos - operator0)) < 0
       dS_pos = dS_pos + step*operator0;
       eta0 = eta0 + step;
    end
    S_pos = S_pos + eta0*perator0;
    
    dS_pos = double(S_pos);
    eta1 = 0;
    while min(eig(dS_pos - operator1)) < 0
       dS_pos = dS_pos + step*operator1;
       eta1 = eta1 + step;
    end
    S_final = S_pos + eta1*operator1;
    
    % Check that all the constraints are verified
    [~, flagPos] = chol(S_final);
    [~, flagConstr0] = chol(S_final - operator0);
    [~, flagConstr1] = chol(S_final - operator1);
    assert(flagPos == 0);
    assert(flagConstr0 == 0);
    assert(flagConstr1 == 0);    
    assert(sum(lambdas_final{1}) + sum(lambdas_final{2}) -1 <= 0) % lambdas sum to 1
    switch supermapClass
        case 2 % FO-supermaps
            assert(0 == max(max(abs(S_final - project_onto_dual_QCFO_superops(S_final, d, spaces)))))
        otherwise % General supermaps
            assert(0 == max(max(abs(S_final - project_onto_dual_valid_superops(S_final, d, spaces)))))
    end
end

