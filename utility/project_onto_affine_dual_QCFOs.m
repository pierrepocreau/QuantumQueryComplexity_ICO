function W_proj = project_onto_affine_dual_QCFOs(W, dims, parties)
%project_onto_affine_dual_QCFOs projects a superoperator onto the (unnormalised) 
% affine dual of valid QC-FOs.
% This projector can be written in terms of the projector onto the space of valid QC-FOs L^{FO} 
% that defines the linear QC-FO constraints (but doesn't require positive semidefiniteness)

    W_proj = W - project_onto_QCFOs(W, dims, parties) + tr_replace(W, 1:length(dims), dims);

end