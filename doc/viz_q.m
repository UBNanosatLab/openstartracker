%Extract the X, Y and Z coordinates
%M=vertcat(sqrt(.5)*[1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1],sqrt(1/3)*[1 1 0;1 -1 0;-1 1 0;-1 -1 0;1 0 1;1 0 -1;-1 0 1;-1 0 -1;0 1 1;0 1 -1;0 -1 1;0 -1 -1],[.5 .5 .5;.5 .5 -.5;.5 -.5 .5;.5 -.5 -.5;-.5 .5 .5;-.5 .5 -.5;-.5 -.5 .5;-.5 -.5 -.5]);
sz=100;
[x,y,z]=sphere(sz);
M=[x(:),y(:),z(:)];
M=normr(horzcat(M,max(abs(M),[],2)));
%trisurf(delaunay(M(:,1:3)),M(:,1),M(:,2),M(:,3))
%surf(reshape(x(:),[sz+1,sz+1]),reshape(y(:),[sz+1,sz+1]),reshape(z(:),[sz+1,sz+1]))
surf(reshape(M(:,1),[sz+1,sz+1]),reshape(M(:,2),[sz+1,sz+1]),reshape(M(:,3),[sz+1,sz+1]))
pbaspect([1 1 1])
