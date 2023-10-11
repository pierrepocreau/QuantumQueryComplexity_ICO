function [S_final,lambdas_frac] = extract_dual(S, lambdas, n, func, T, supermapClass, symbolic)
% Extract an exact dual solution from a numerical approximation. 
% This function gives an upper bound on the objective (1-epsilon) and therefore a lower bound on epsilon.
% It follows Algorithm 1 as described in the paper:
% A. A. Abbott, M. Mhalla, P. Pocreau, "Quantum Query Complexity of Boolean Functions under Indefinite Causal Order", arXiv:2307.10285 

    dim_H = n+1;
    
    d = dim_H * ones(1,2*T);
    dim = prod(d);
    spaces{1}{1} = []; % trivial P
    for i = 1:T
        spaces{i+1}{1} = 2*i - 1;
        spaces{i+1}{2} = 2*i;
    end
    spaces{T+2}{1} =  []; % trivial F
    
    % Verify all lambdas are non-negative
    for i = 1:2
        lambdas{i} = max(lambdas{i},zeros(length(lambdas{i}),1));
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
        all_lambdas = sort([lambdas{1}; lambdas{2}],'descend');
        N = 2^n;
        while all_lambdas(N) < -delta/N
            N = N - 1;
        end

        % Adjust the lambdas
        lambdas1 = lambdas_frac{1};
        lambdas_final1 = lambdas1;
        for x = 1:length(lambdas1)
            if lambdas1(x) >= -delta/N
                lambdas_final1(x) = lambdas1(x) + delta/N;
            end
        end
        lambdas_final{1} = lambdas_final1;

        lambdas2 = lambdas_frac{2};
        lambdas_final2 = lambdas2;
        for x = 1:length(lambdas2)
            if lambdas2(x) >= -delta/N
                lambdas_final2(x) = lambdas2(x) + delta/N;
            end
        end
        lambdas_final{2} = lambdas_final2;
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
    
    % Compute the expressions of the operator0 and operator1 (O^{[0]} and O^{[1}}) in the paper
    oracles = oracles_map(n, T);
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
    
    % First we make sure S >= 0 (this ensures the convergence of the following steps)
    if symbolic == true
        S_valid_double = double(S_valid);
        lambda_min = min(eig(S_valid_double));
    else
        lambda_min = min(eig(S_valid));
    end
    mu = lambda_min/(lambda_min-1) + 10^-9;
    S_pos = (1-mu)*S_valid + mu*eye(dim);

    % Make sure that the conditions S - operator0 >= 0 and S - operator1 >= 0 are respected. 
    step = 1*10^-5;

    if symbolic == true
        % Save the symbolic version of S, operatorX, then work numerically for now
        S_pos_sym = S_pos;
        S_pos = double(S_pos);
        operator0_sym = operator0;
        operator0 = double(operator0);
        operator1_sym = operator1;
        operator1 = double(operator1);
    end
    
    eta0 = 0;
    while min(eig(S_pos + (eta0-1)*operator0)) < 0
       eta0 = eta0 + step;
    end

    eta1 = 0;
    while min(eig(S_pos + (eta1-1)*operator1)) < 0
       eta1 = eta1 + step;
    end

    if symbolic == true
        S_final = S_pos_sym + sym(eto0)*operator0_sym + sym(eta1)*operator1_sym;
        operator0 = operator0_sym;
        operator1 = operator1_sym;
    else
        S_final = S_pos + eto0*operator0 + eta1*operator1;
    end
    
    % Check that all the constraints are verified
    [~, flagConstr0] = chol(S_final - operator0);
    [~, flagConstr1] = chol(S_final - operator1);
    assert(flagConstr0 == 0);
    assert(flagConstr1 == 0);    
    assert(sum(lambdas_final{1}) + sum(lambdas_final{2}) -1 <= 0) % lambdas sum to 1
    assert(all([lambdas_final{1}; lambdas_final{2}] >= 0)); % All lambdas non-negative
    switch supermapClass
        case 2 % FO-supermaps
            assert(0 == max(max(abs(S_final - project_onto_dual_QCFO_superops(S_final, d, spaces)))))
        otherwise % General supermaps
            assert(0 == max(max(abs(S_final - project_onto_dual_valid_superops(S_final, d, spaces)))))
    end
end

