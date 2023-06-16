function W = gen_variables(dim_H,T,symmetries)
    if symmetries
        %Computes the representation of the permutation group on T elements
        %as a permutation of the Hilbert spaces.
        G = replab.S(T);
        elements = G.generators();
        images = {};
        for element = elements
            images{end + 1} = perm_to_perm(element{1}, dim_H, T);
        end

        dim = dim_H^(2*T);
        W{1} = replab.CommutantVar.fromPermutations(images, 'symmetric', 'real');
        W{2} = replab.CommutantVar.fromPermutations(images, 'symmetric', 'real');
        
    else
        dim = dim_H^(2*T);
        W{1} = sdpvar(dim,dim,'symmetric');
        W{2} = sdpvar(dim,dim,'symmetric');
    end
end

