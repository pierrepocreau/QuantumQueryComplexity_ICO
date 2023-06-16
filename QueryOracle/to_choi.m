function F = to_choi(func, input_dim)
% Creates the Choi matrix of the function func: H -> H, where dim(H) = n.
% The register is composed of log_2(n) qubits
    temp = zeros(input_dim^2, 1);

    % Calculation of f|1>>
    for k = 1:input_dim
        %This loop over the computational basis.
        vec = zeros(input_dim, 1);
        vec(k) = 1;
        q = func(vec);
        temp = temp + kron(vec, q);
    end
    %f|1>><<1|f'
    F = temp * temp';
end