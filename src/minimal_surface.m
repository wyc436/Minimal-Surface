clc
clear

[x, t] = readObj('../meshes/Cat_head.obj');
% [x, t] = readObj('../meshes/Balls.obj');
% [x, t] = readObj('../meshes/Bunny_head.obj');

%% draw 2 copies of the image
figure; set(gcf, 'Units', 'normalized', 'Position', [0.05,0.05,.8,.8]);
subplot(121); trimesh(t, x(:,1), x(:,2), x(:,3), 'edgecolor', 'k'); axis off; axis equal; title('input');
subplot(122); h = trimesh(t, x(:,1), x(:,2), x(:,3), 'edgecolor', 'k'); axis off; axis equal; title('output');

%% TODO: find boundary and interior vertices
nv=size(x,1);
nf=size(t,1);

%Initializes the array of struct
[vertices(1:nv).numadjpoints]=deal(0);
[vertices(1:nv).neighbortriangles]=deal([]);
[vertices(1:nv).neighborpoints]=deal([]);

%Get the adjacency trangle facets of each point
for i=1:nf
    for k=1:3
        vertices(t(i,k)).neighbortriangles(end+1)=i;       
    end
end

%Get boundary vertices
for i=1:nv
    for j=1:size(vertices(i).neighbortriangles,2)
        temp=vertices(i).neighbortriangles(j);
        for k=1:3
            if t(temp,k)~=i && not(any(vertices(i).neighborpoints==t(temp,k)))
                vertices(i).neighborpoints=[vertices(i).neighborpoints,t(temp,k)];
                vertices(i).numadjpoints=vertices(i).numadjpoints+1;
            elseif t(temp,k)~=i 
                vertices(i).neighborpoints=setdiff(vertices(i).neighborpoints,t(temp,k));
                vertices(i).numadjpoints=vertices(i).numadjpoints-1;
            end
        end
    end
end
edgepoints=find([vertices(:).numadjpoints]);
    
%% TODO: compute Laplacian
L = laplacian(x, t,edgepoints);

%% TODO: compute minimal surface using LOCAL approach
L=L-diag(diag(L));
ind=sub2ind([nv,nv],edgepoints(:),edgepoints(:));
L(ind)=-1;
delta=x+L*x;
while any(delta(:)~=0)
    x=x-delta;
    delta=x+L*x;
end

%% TODO: compute minimal surface using GLOBAL approach
% b=zeros(nv,3);
% b(edgepoints(:),:)=x(edgepoints(:),:);
% x=L\b;

%% draw mesh in minimal surface iterations
set(h, 'Vertices', x);