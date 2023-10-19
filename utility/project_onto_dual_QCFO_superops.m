function W_proj = project_onto_dual_QCFO_superops(W, dims, parties)
%project_onto_dual_QCFO_superops projects a superoperator onto the subspace
%of unnormalised process in the dual space of the QCFOs.
%This projector can be written from the projector on the space of QCFOs.

    W_proj = W - project_onto_QCFOs(W, dims, parties) + tr_replace(W, 1:length(dims), dims);

end