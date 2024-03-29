clear
clc

w = 0.30000000;
K = [207.84609691 0.00000000 160.00000000; 0.00000000 207.84609691 120.00000000; 0.00000000 0.00000000 1.00000000];
q = [55.32846006 104.12408043 91.95464047 40.99451795; 127.18065960 130.50046628 191.07305342 182.49729124];
T = [0.99260833 0.01424885 0.12052252 0.01127236; -0.02872118 0.99245349 0.11921059 0.09066921; -0.11791438 -0.12179097 0.98552684 1.09961775; 0.00000000 0.00000000 0.00000000 1.00000000];
mu = 0.01000000;

M = T;

q = [q; 1 1 1 1];
p1 = [-0.5*w; -0.5*w];
p2 = [0.5*w; -0.5*w];
p3 = [0.5*w; 0.5*w];
p4 = [-0.5*w; 0.5*w];
p = [p1 p2 p3 p4; 0 0 0 0; 1 1 1 1];
PI_0 = [1 0 0 0; 0 1 0 0; 0 0 1 0];

[~,npoints]=size(q);

for x = 1:npoints
    [A{x},B{x},C{x}]=calc_test(M,K,PI_0,q(:,x),p(:,x));
end

sum_A = sum(cat(3,A{:}),3);
sum_B = sum(cat(3,B{:}),3);
sum_C = sum(cat(3,C{:}),3);
u = inv(sum_A + 0.01 * eye(6))*sum_B;
M = expm(ss61(u))*M;
disp(mat2str(M))
