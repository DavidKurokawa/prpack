function [x,ret] = pagerank(A)

assert(size(A,1)==size(A,2));
n=size(A,1);
[i,j,v] = find(A);
i = int32(i);
j = int32(j);
S = pagerank_solver(int32(n),i,j);
u = ones(n,1)./n;
v = ones(n,1)./n;
[x,ret] = S.solve(0.85,1e-10,u,v,'gs');
