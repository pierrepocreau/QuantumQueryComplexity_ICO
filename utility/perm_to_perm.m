function [ H_perm ] = perm_to_perm(permutation, n, T)
    % Representation of a permutation on T elements to a permutation on the 
    % n^(2*T) dimensions of our Hilbert space.
    dim = n^(2 * T);
    H_perm = 1:dim;

    for j = 1:dim
        init_vec = zeros(1,dim);
        init_vec(j) = 1;
        perm_vec = PermuteSystems(init_vec, permutation);

        %Permute the 1 in the vector representation of the ket vector.
        id_init = find(init_vec, 1);
        id_perm = find(perm_vec, 1);
        H_perm(id_init) = id_perm;
    end
end
