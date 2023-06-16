bits=4;
amb = boolean([0 1 1 1 0 1 0 0 0 0 1 0 1 1 1 0]);

                            %pos
assert(amb([1,1,0,0]) == 1) %13
assert(amb([0,0,1,0]) == 1) %3
assert(amb([1,0,1,0]) == 1) %11
assert(amb([1,1,1,0]) == 1) %15
assert(amb([0,0,0,1]) == 1) %2
assert(amb([0,1,0,1]) == 1) %6
assert(amb([1,1,0,1]) == 1) %14
assert(amb([0,0,1,1]) == 1) %4

% Comparison with the polynomial expression

for x = dec2bin(0:2^bits-1)' - '0'
    assert(amb(x') == pol(x'));
    assert(amb(x') == computing_with_eq(x'));
end
%%
computing_with_eq([1 1 0 0])
%%


function y = pol(x)
 x1 = x(4);
 x2 = x(3);
 x3 = x(2);
 x4 = x(1);
 y = mod(x1 + x2 + x3*x4 - x1*x4 - x2*x3 - x1*x2,2);
end

function y = computing_with_eq(x)
 x1 = x(4);
 x2 = x(3);
 x3 = x(2);
 x4 = x(1);

 a = x1 == x2;
 b = x3 == x4;
 if ~a && b
     y = 1;
 else
     y=1;
 end
end
