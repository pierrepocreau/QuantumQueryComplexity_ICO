function handle = boolean(tt)
    % Return the function define by the truth table tt.
    handle = @(x) eval(tt, x);
end

function y = eval(tt, x)
    % Evaluate the function defined by its truth table tt, on input x.
    if size(output, 2) ~=  2^size(x, 2)
        bits_to_cut =  log(2^size(x, 2) / size(tt,2))/log(2);
        x = x(1:end-bits_to_cut);
    end
    index = bin2dec(num2str(x))+1;
    y = tt(index);
end