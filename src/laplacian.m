function L = laplacian(x, t, edgepoints)
% compute cotanget Laplacan of mesh (x, t)

nv = size(x, 1);
nf = size(t, 1);
L = sparse(nv, nv);
ind=sub2ind([nv,nv],1:nv,1:nv);
L(ind)=1;

%Initializes the array of struct
[vertices(1:nv).neighbortriangles]=deal([]);
[vertices(1:nv).neighborpoints]=deal([]);

%Get the adjacency trangle facets of each point
for i=1:nf
    for j=1:3
        vertices(t(i,j)).neighbortriangles=[vertices(t(i,j)).neighbortriangles,i];
    end
end

%Get the adjacency point of each point
for i=1:nv
    temp=t(vertices(i).neighbortriangles(:),:);
    vertices(i).neighborpoints=union(vertices(i).neighborpoints,temp);
    vertices(i).neighborpoints=(setdiff(vertices(i).neighborpoints,i))';
end

% compute cotanget Laplacan of mesh (x, t)
for i=1:nv
    if all(edgepoints~=i) 
       numadjpoints=size(vertices(i).neighborpoints,2); 
       weights=zeros(1,numadjpoints);
       for j=1:numadjpoints 
           k=vertices(i).neighborpoints(j);
           tri=[];
           for l=1:size(vertices(i).neighbortriangles,2)
               m=vertices(i).neighbortriangles(l);
               if any(t(m,:)==i) && any(t(m,:)==k)
                   tri=[tri;t(m,:)];
               end
           end
           c=sum(tri(1,:),2)-i-k;
           d=sum(tri(2,:),2)-i-k;
           alpha=subspace((x(i,:)-x(c,:))',(x(k,:)-x(c,:))');
           beta=subspace((x(i,:)-x(d,:))',(x(k,:)-x(d,:))');
           weights(j)=cot(alpha)+cot(beta);
       end 
       L(i,vertices(i).neighborpoints(:))=-weights/sum(weights,2);
    end                    
end


