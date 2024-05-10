function q=PlotWeb2D(Nod, I, J)
% Plotting points
h=plot(Nod(:,1),Nod(:,2),'.');
v=Nod(J,[1,2])-Nod(I,[1,2]);

% Plotting lines
hold on;
q = quiver(Nod(I,1),Nod(I,2),v(:,1),v(:,2),0);
q.ShowArrowHead = 'off';
hold off;
axis equal

end