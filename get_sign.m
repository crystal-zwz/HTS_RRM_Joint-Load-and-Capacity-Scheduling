%{
Function: get_sign
Description: 求实数x的符号
Input: 实数x
Output: 实数x的符号y
Author: Marc Pony(marc_pony@163.com)
%}
function y = get_sign(x)
if abs(x) < 1.0e-8
    y = 0;
else
    if x < 0.0
        y = -1;
    else
        y = 1;
    end
end
end
