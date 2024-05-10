function q=PlotWeb3D(Nod, I, J)
% Plotting points
h=plot3(Nod(:,1),Nod(:,2),Nod(:,3),'.');
v=Nod(J,:)-Nod(I,:);

% Plotting lines
hold on;
q = quiver3(Nod(I,1),Nod(I,2),Nod(I,3),v(:,1),v(:,2),v(:,3),0);
q.ShowArrowHead = 'off';
hold off;
axis equal

view(0,90)
end