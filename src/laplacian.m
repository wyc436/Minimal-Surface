function L = laplacian(x, t)
% compute cotanget Laplacan of mesh (x, t)
nv = size(x, 1);

e1=x(t(:,1),:)-x(t(:,2),:);
e2=x(t(:,3),:)-x(t(:,1),:);
e3=x(t(:,2),:)-x(t(:,3),:);

alpha1=acos(-(dot(e1',e2'))'./sqrt(sum(e1.^2,2).*sum(e2.^2,2)));
alpha2=acos(-(dot(e2',e3'))'./sqrt(sum(e2.^2,2).*sum(e3.^2,2)));
alpha3=acos(-(dot(e1',e3'))'./sqrt(sum(e1.^2,2).*sum(e3.^2,2)));

L = sparse( t(:,[2 1 3]), t(:,[3 2 1]), cot([alpha1,alpha2,alpha3]), nv, nv );
L=L+L';
L=-spdiags(-sum(L,2),0,L);






