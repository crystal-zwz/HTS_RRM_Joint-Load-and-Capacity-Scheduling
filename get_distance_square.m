%{
Function: get_distance_square
Description: 求平面两点之间距离的平方
Input: 平面两点a，b
Output: 平面两点之间距离的平方
Author: Marc Pony(marc_pony@163.com)
%}
function distanceSquare = get_distance_square(a, b)
distanceSquare = (a(1) - b(1))^2 + (a(2) - b(2))^2;
end
