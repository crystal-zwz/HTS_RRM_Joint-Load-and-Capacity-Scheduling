%{
Function: min_cover_circle
Description: 求平面pointCount个点的最小覆盖圆
Input: 平面pointCount个点的坐标(x,y)，点个数pointCount
Output: 平面pointCount个点的最小覆盖圆圆心center，半径radius
Author: Marc Pony(marc_pony@163.com)
%}
function [center, radius] = min_cover_circle(x, y, pointCount)
p = [x(:)'; y(:)'];
p = p(:, randperm(pointCount)); %随机打乱数据
% disp(p)
center = p(:, 1);
radiusSquare = 0.0;
for i = 2 : pointCount
    if get_sign(get_distance_square(p(:, i), center) - radiusSquare) > 0
        center = p(:, i);
        radiusSquare = 0.0;
        for j = 1 : i - 1
            if get_sign(get_distance_square(p(:, j), center) - radiusSquare) > 0
                center = 0.5 * (p(:, i) + p(:, j));
                radiusSquare = get_distance_square(p(:, j), center);
                for k = 1 : j - 1
                    if get_sign(get_distance_square(p(:, k), center) - radiusSquare) > 0
                        center = get_circle_center(p(:, i), p(:, j), p(:, k));
                        radiusSquare = get_distance_square(p(:, i), center);
                    end
                end
            end
        end
    end
end
radius = sqrt(radiusSquare);
end

