function W_proj = project_onto_affine_dual_superops(W, dims, parties)
%project_onto_affine_dual_superops projects a superoperator onto the (unnormalised) 
% affine dual of valid superoperators.
% This projector can be written in terms of the projector onto the space of valid superops L^{Gen} 
% that defines the linear superoperator constraints (but doesn't require positive semidefiniteness)

    W_proj = W - project_onto_valid_superops(W, dims, parties) + tr_replace(W, 1:length(dims), dims);
end
