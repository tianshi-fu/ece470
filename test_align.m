clear
clc

w = 0.13000000;
K = [207.84609691 0.00000000 160.00000000; 0.00000000 207.84609691 120.00000000; 0.00000000 0.00000000 1.00000000];
q = [117.67480026 142.53118100 126.82343763 103.92471824; 136.81132122 151.92401897 177.17978438 162.50916168];

q = [q; 1 1 1 1];
p1 = [-0.5*w; -0.5*w];
p2 = [0.5*w; -0.5*w];
p3 = [0.5*w; 0.5*w];
p4 = [-0.5*w; 0.5*w];
p = [p1 p2 p3 p4; 0 0 0 0; 1 1 1 1];
PI_0 = [1 0 0 0; 0 1 0 0; 0 0 1 0];

n = 2;
trials = 1;
[~,npoints]=size(q);


while trials < 51
rng shuffle
posi = randn(3,1);
[R,~]=qr(randn(3));
M = [R posi; 0 0 0 1];
tic
while n > 1e-10
for x = 1:npoints
    [A{x},B{x},C{x}]=calc_test(M,K,PI_0,q(:,x),p(:,x));
end

sum_A = sum(cat(3,A{:}),3);
sum_B = sum(cat(3,B{:}),3);
sum_C = sum(cat(3,C{:}),3);

u = inv(sum_A + 0.01 * eye(6))*sum_B;
M = expm(ss61(u))*M;
n = sum_C;

Ms{trials}=M;
nu(trials)=norm(u);

e = toc;
if toc > 0.5
    break
end


end
disp(trials)
errors(trials)=sum_C;
trials = trials+1;
end
[~,min_ind]=min(errors);
min_M = Ms{min_ind};
disp(mat2str(min_M))