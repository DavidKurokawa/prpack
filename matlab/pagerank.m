function [x,ret] = pagerank(A)

assert(size(A,1)==size(A,2));
n=size(A,1);
[i,j,v] = find(A);
i = int32(i-1);
j = int32(j-1);
S = pagerank_solver(int32(n),i,j);
u = [];
v = [];
[x,ret] = S.solve(0.85,1e-10,u,v,'gs');
