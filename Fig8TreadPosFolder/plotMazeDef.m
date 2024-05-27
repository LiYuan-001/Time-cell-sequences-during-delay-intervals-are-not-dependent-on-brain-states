function plotMazeDef(maze)

plot(maze.regx.N1(:,1),maze.regx.N1(:,2))
hold on
plot(maze.regx.N2(:,1),maze.regx.N2(:,2))
plot(maze.regx.N3(:,1),maze.regx.N3(:,2))
plot(maze.regx.N4(:,1),maze.regx.N4(:,2))
plot(maze.regx.N5(:,1),maze.regx.N5(:,2))
plot(maze.regx.N6(:,1),maze.regx.N6(:,2))

plot(maze.locs.A23(:,1),maze.locs.A23(:,2))
hold on
plot(maze.locs.A34(:,1),maze.locs.A34(:,2))
plot(maze.locs.A45(:,1),maze.locs.A45(:,2))
plot(maze.locs.A56(:,1),maze.locs.A56(:,2))
plot(maze.locs.A16(:,1),maze.locs.A16(:,2))
plot(maze.locs.A12(:,1),maze.locs.A12(:,2))
plot(maze.locs.A25(:,1),maze.locs.A25(:,2))
hold off
end