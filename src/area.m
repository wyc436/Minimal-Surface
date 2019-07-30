function area = area(x, t)
% Calculating the area of mesh (x, t)

e1=x(t(:,1),:)-x(t(:,2),:);
e2=x(t(:,3),:)-x(t(:,2),:);
area=sum(sqrt(sum(cross(e1,e2).^2,2)));

