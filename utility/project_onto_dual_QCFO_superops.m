function W_proj = project_onto_dual_QCFO_superops(W, dims_raw, parties_raw)
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

            W_proj = W - (tr_replace(W, BO, dims) - tr_replace(W, B, dims));
            W_proj = W_proj - (tr_replace(W_proj, [AO, B], dims) - tr_replace(W_proj, [A, B], dims));

        otherwise
            error('Projector onto valid superoperators currently only implemented up to N=2.');
    end

    % Put W back in its original ordering
    W_proj = superop_from_canonical_ordering(W_proj,dims_raw,parties_raw);
end