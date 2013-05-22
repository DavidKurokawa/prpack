n = 9914;
seed = 0;
rand('state',seed);
u = rand(n,1);
u = u/sum(u);

dlmwrite('csstan-u.vec',u,'precision','%.18e');

v1 = zeros(n,1);
v1(1) = 1;
dlmwrite('csstan-e1.vec',v1);

v2 = rand(n,1);
v2 = v2/sum(v2);
dlmwrite('csstan-v.vec',v2,'precision','%.18e');