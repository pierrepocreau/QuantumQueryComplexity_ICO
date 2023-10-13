function W_final = extract_primal(W, T, n, supermapClass, symbolic)
% Extract an exact primal solution from a numerical approximation. 
% This function gives a lower bound on the objective (1-epsilon) and therefore an upper bound on epsilon.
% It follows Algorithm 2 as described in the paper:
% A. A. Abbott, M. Mhalla, P. Pocreau, "Quantum Query Complexity of Boolean Functions under Indefinite Causal Order", arXiv:2307.10285 

    dim_H = n+1;

    d = dim_H * ones(1,2*T);
    d = [1 d 1]; % Trivial input and output spaces
    dim = prod(d);

    spaces{1}{1} = [];    
    spaces{1}{1} = 1;
    for i = 1:T
        spaces{i+1}{1} = 2*i;
        spaces{i+1}{2} = 2*i + 1;
    end
    spaces{T+2}{1} =  2*T+2;
    spaces{T+2}{2} = [];
    
    % Rationalise each element of the superinstrument W
    [N, D] = rat(W{1});
    W_frac{1} = N./D;
    [N, D] = rat(W{2});
    W_frac{2} = N./D;
    
    % Use a symbolic representation
    if symbolic == true
        W_frac{1} = sym(W_frac{1});
        W_frac{2} = sym(W_frac{2});
    end
    
    % Ensure that the superinstrument elements are Hermitian
    W_Herm{1} = (W_frac{1} + W_frac{1}') / 2;
    W_Herm{2} = (W_frac{2} + W_frac{2}') / 2;

    % Project W onto the space of valid QCFO or general supermaps
    switch supermapClass
        case 2 % QC-FO
            W_proj = project_onto_QCFOs(W_Herm{1} + W_Herm{2}, d, spaces);
        otherwise % General supermaps
            W_proj = project_onto_valid_superops(W_Herm{1} + W_Herm{2}, d, spaces);
    end
    
    % Ensure that W is in the valid space of supermaps by removing the orthogonal part 
    W_corr = W_proj - W_Herm{1} - W_Herm{2};
    W_valid{1} = W_Herm{1} + W_corr / 2;
    W_valid{2} = W_Herm{2} + W_corr / 2;
    
    % Ensure that the superinstrument elements are positive semi-definite
    % by mixing them with the identity.
    % This step can be performed analytically, however computing the
    % eigenvalues in symbolic mode is computationally inefficient.
    % Here we compute the eigenvalue numerically and shift them by a small constant.
    if symbolic == true
        W_double{1} = double(W_valid{1});
        W_double{2} = double(W_valid{2});
        lambda_min = min(min(eig(W_double{1})), min(eig(W_double{2})));
    else
        lambda_min = min(min(eig(W_valid{1})), min(eig(W_valid{2})));
    end
    mu = lambda_min/(lambda_min - 1) + 10^-8;

    if mu > 0
        W_pos{1} = (1-mu)*W_valid{1} + mu*eye(dim);
        W_pos{2} = (1-mu)*W_valid{2} + mu*eye(dim);
    end

    % Renormalise the process matrix
    norm = trace(W_pos{1} + W_pos{2});
    W_final{1} = dim_H^T/norm * W_pos{1};
    W_final{2} = dim_H^T/norm * W_pos{2};

    % Check that all the constraints are verified
    [~, flagPositivity1] = chol(W_final{1});
    [~, flagPositivity2] = chol(W_final{2});
    assert(flagPositivity1==0);
    assert(flagPositivity2==0);    
    switch supermapClass
        case 2 % QC-FO
            assert(0 == max(max(abs(W_final{1} + W_final{2} - project_onto_QCFOs(W_final{1} + W_final{2}, d, spaces)))));
        otherwise % General supermaps
            assert(0 ==  max(max(abs(W_final{1} + W_final{2} - project_onto_valid_superops(W_final{1} + W_final{2}, d, spaces)))));
    end

end