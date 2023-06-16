%% Test of the second function Ambainis propose in his paper.
%it evaluate to 0 if all variables are equals.
bits=2;
and_f = boolean([0 0 0 1]);

% Comparison with the polynomial expression

for x = dec2bin(0:2^bits-1)' - '0'
    assert(and_f(x') == pol(x'));
end
%%

function y = pol(x)
 x1 = x(2);
 x2 = x(1);
 y = mod(x1*x2, 2);
end

