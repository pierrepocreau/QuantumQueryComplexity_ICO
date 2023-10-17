function W_proj = project_onto_dual_valid_superops(W, dims_raw, parties_raw)
    %% Setup and process the input

    % First put W in canonical ordering (this checks the input validity too)
    % The spaces P,AI,AO,...,F then correspond to dims 1,2,3,...,2*N+2
    if exist('parties_raw','var') && ~isempty(parties_raw)
        [W, dims, parties] = superop_to_canonical_ordering(W, dims_raw, parties_raw);
    else
        [W, dims, parties] = superop_to_canonical_ordering(W, dims_raw);
        parties_raw = parties;
    end

    N = length(parties) - 2;
    
    %% Implement projection separately and explicitly for each N
    switch N
        case 1
            % Not entirely sure it's the correct projector
            P = 1;
            AI = 2;
            AO = 3;
            F = 4;

            w_proj_AO = tr_replace(W,AO,dims);
                   
            w_proj_AIAO= tr_replace(w_proj_AO, AI, dims);
            
            W_proj = W - w_proj_AO + w_proj_AIAO;
        case 2
            %%
            %Assumes dim P and F are 1 ?
            P = 1;
            AI = 2;
            AO = 3;
            A = [AI, AO];
            BI = 4;
            BO = 5;
            B = [BI, BO];
            F = 6;

            w_proj_AO = tr_replace(W,AO,dims);
            w_proj_BO = tr_replace(W,BO,dims);
                    
            w_proj_AIAO= tr_replace(w_proj_AO, AI, dims);
            w_proj_BIBO= tr_replace(w_proj_BO, BI, dims);
            w_proj_AOBO = tr_replace(W, [AO, BO], dims);

            w_proj_AOBIBO = tr_replace(w_proj_AOBO, BI, dims);
            w_proj_AIAOBO = tr_replace(w_proj_AOBO, AI, dims);

            w_proj_AIAOBIBO = tr_replace(w_proj_AOBIBO, BO, dims);
            
            W_proj = W - w_proj_AO - w_proj_BO + w_proj_AIAO + w_proj_BIBO + w_proj_AOBO - w_proj_AOBIBO - w_proj_AIAOBO + w_proj_AIAOBIBO;
        otherwise
            error('Projector onto valid superoperators currently only implemented for N=2.');
    end

    % Put W back in its original ordering
    W_proj = superop_from_canonical_ordering(W_proj,dims_raw,parties_raw);

end
