clc
clear

[x, t] = readObj('../meshes/Cat_head.obj');
[x, t]=subdivdision(x, t, 3);

%% draw 2 copies of the image
figure; set(gcf, 'Units', 'normalized', 'Position', [0.05,0.05,.8,.8]);
subplot(121); trimesh(t, x(:,1), x(:,2), x(:,3), 'edgecolor', 'k'); axis off; axis equal; title('input');
subplot(122); h = trimesh(t, x(:,1), x(:,2), x(:,3), 'edgecolor', 'k'); axis off; axis equal; title('output');

%% TODO: find boundary and interior vertices
nv=size(x,1);
nf=size(t,1);

%% Get boundary vertice id
Edge = sparse(t, t(:, [2 3 1]), true, nv, nv);
[edgepoints,e2] = find(xor( Edge, Edge'));

for i=1:5
    %% TODO: compute Laplacian
    L = laplacian(x, t);
    L(edgepoints(:),:)=0;
    ind=sub2ind([nv,nv],edgepoints(:),edgepoints(:));
    L(ind)=1;

    %% TODO: compute minimal surface using LOCAL approach
%     L=L./full(diag(L));
%     L = -spdiags(zeros(nv,1), 0, L);
%     L(ind)=1;
%     x2=0.7*x+0.3*L*x;
%     while norm(x-x2)>1e-3
%         x=x2;
%         x2=0.7*x+0.3*L*x;
%         set(h, 'Vertices', x); drawnow; pause(0.01);
%     end

    %% TODO: compute minimal surface using GLOBAL approach
    b=zeros(nv,3);
    b(edgepoints(:),:)=x(edgepoints(:),:);
    x=L\b;
end
%% draw mesh in minimal surface iterations
set(h, 'Vertices', x);