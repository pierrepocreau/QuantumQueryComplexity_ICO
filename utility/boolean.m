function handle = boolean(output)
    handle = @(x) eval(output, x);
end

function y = eval(output, x)
    % A function is defined by its 16 bits output vector.
    if size(output, 2) ~=  2^size(x, 2)
        bits_to_cut =  log(2^size(x, 2) / size(output,2))/log(2);
        x = x(1:end-bits_to_cut);
    end
    index = bin2dec(num2str(x))+1;
    y = output(index);
end