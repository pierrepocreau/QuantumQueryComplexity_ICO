function W_proj = project_onto_dual_valid_superops(W, dims_raw, parties_raw)
%project_onto_dual_valid_superops projects a superoperator onto the subspace
%of unnormalised process in the dual space of the valid superoperators.
%It corresponds to the space of unnormalised non-signalling channels.
%This projector can be written from the projector on the space of valid superoperators.

    W_proj = W - project_onto_valid_superops(W, dims_raw, parties_raw) + tr_replace(W, 1:length(dims_raw), dims_raw);
end
