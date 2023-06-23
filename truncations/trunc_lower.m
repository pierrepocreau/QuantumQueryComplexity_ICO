function W_final = trunc_lower(W, T, supermapClass, symbolic)
%Modifies the process W obtained after SDP optimisation so that it 
%verifies the constraints exactly. It will provide a lower bound of the
%objectif we aim to maximise.

    dim = size(W{1}, 1);
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
    
    %This rationalization procedure cleans the matrix.
    [N, D] = rat(W{1});
    W_truncated{1} = N./D;
    [N, D] = rat(W{2});
    W_truncated{2} = N./D;
    
    %Use a symbolic representation.
    if symbolic == true
        W_truncated{1} = sym(W_truncated{1});
        W_truncated{2} = sym(W_truncated{2});
    end
    
    %Assure that the process is self-adjoint.
    W_autoAdj{1} = (W_truncated{1} + W_truncated{1}') / 2;
    W_autoAdj{2} = (W_truncated{2} + W_truncated{2}') / 2;

    %Assure that we have a valid process, by projecting it on the space.
    switch supermapClass
    case 2 % QC-FO
        W_valid = project_onto_QCFOs(W_autoAdj{1} + W_autoAdj{2}, d, A);
    otherwise % General supermaps
        W_valid = project_onto_valid_superops(W_autoAdj{1} + W_autoAdj{2}, d, A);
    end
    
    W_rest = W_valid - W_autoAdj{1} - W_autoAdj{2};
    W_ok{1} = W_autoAdj{1} + W_rest / 2;
    W_ok{2} = W_autoAdj{2} + W_rest / 2;
    
    %The matrices must be positive
    %Computing the eigenvalues in symbolic mode takes too much time, so we
    %computes them in numeric and shifts them a little bit.
    if symbolic == true
        W_double{1} = double(W_ok{1});
        W_double{2} = double(W_ok{2});
        vp_min = min(min(eig(W_double{1})), min(eig(W_double{2})));
        eps = vp_min/(vp_min - 1) + 10^-8;
    else
        vp_min = min(min(eig(W_ok{1})), min(eig(W_ok{2})));
        eps = vp_min/(vp_min - 1) + 10^-8;
    end

    if eps > 0
        W_ok{1} = (1-eps)*W_ok{1} + eps*eye(dim);
        W_ok{2} = (1-eps)*W_ok{2} + eps*eye(dim);
    end

    %Renormalize
    norm = trace(W_ok{1} + W_ok{2});
    W_final{1} = dim_H^T/norm * W_ok{1};
    W_final{2} = dim_H^T/norm * W_ok{2};

    % Check the process is indeed in the right space -> that projection is idempotent
    switch supermapClass
        case 2 % QC-FO
            assert(10^-8 >= max(max(abs(W_final{1} + W_final{2} - project_onto_QCFOs(W_final{1} + W_final{2}, d, A)))));
        otherwise % General supermaps
            assert(10^-8 >=  max(max(abs(W_final{1} + W_final{2} - project_onto_valid_superops(W_final{1} + W_final{2}, d, A)))));
    end
    [~, flag1] = chol(W_final{1});
    [~, flag2] = chol(W_final{2});
    assert(flag1==0);
    assert(flag2==0);
end