%% Test of the second function Ambainis propose in his paper.
%it evaluate to 0 if all variables are equals.
bits=3;
amb2 = boolean([0 1 1 1 1 1 1 0]);

% Comparison with the polynomial expression

for x = dec2bin(0:2^bits-1)' - '0'
    assert(amb2(x') == pol(x'));
    assert(amb2(x') == comput_with_eq(x')); 
end
%%

function y = pol(x)
 x1 = x(3);
 x2 = x(2);
 x3 = x(1);
 y = mod(x1 + x2 + x3 -x1*x2 -x1*x3 - x2*x3, 2);
end

function y = comput_with_eq(x)
 x1 = x(3);
 x2 = x(2);
 x3 = x(1);
 a = x1 == x2;
 b = x2 == x3;
 if a && b
     y = 0;
 else
     y = 1;
 end
end