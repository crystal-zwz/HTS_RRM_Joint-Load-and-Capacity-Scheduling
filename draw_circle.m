%{
Function: draw_circle
Description: 画圆周
Input: 圆心center，圆周半径radius，线型/颜色等设置参数options
Output: 无
Author: Marc Pony(marc_pony@163.com)
%}
function draw_circle(center, radius, options)
theta = 0.0 : 0.001 : 2.0 * pi;
x = center(1) + radius * cos(theta);
y = center(2) + radius * sin(theta);
plot(x, y, options)
end
