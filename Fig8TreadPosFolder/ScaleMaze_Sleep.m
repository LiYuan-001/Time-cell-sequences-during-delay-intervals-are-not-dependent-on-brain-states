function [x_scaled,y_scaled, x_center,y_center,maze] = ScaleMaze_Sleep(x, y, target_maze) 

ind = isnan(x);

x_center = (min(x)+max(x))/2;
y_center = (min(y)+max(y))/2;
x = (x - x_center);
y = (y - y_center);

x_scaled = x * target_maze(1)/(max(x)-min(x));
y_scaled = y * target_maze(2)/(max(y)-min(y));

maze.center = [mean(x) mean(y)];

x_scaled(ind) = NaN;
y_scaled(ind) = NaN;




