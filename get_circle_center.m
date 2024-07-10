%{
Function: get_circle_center
Description: 求三角形外接圆的圆心
Input: 平面三个点a，b，c
Output: 三角形外接圆的圆心center
Author: Marc Pony(marc_pony@163.com)
%}
function center = get_circle_center(a, b, c)
center = zeros(2, 1);
a1 = b(1) - a(1);
b1 = b(2) - a(2);
c1 = 0.5 * (a1 * a1 + b1 * b1);
a2 = c(1) - a(1);
b2 = c(2) - a(2);
c2 = 0.5 * (a2 * a2 + b2 * b2);
d = a1 * b2 - a2 * b1;
center(1) = a(1) + (c1 * b2 - c2 * b1) / d;
center(2) = a(2) + (a1 * c2 - a2 * c1) / d;
end
