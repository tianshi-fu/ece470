clear
clc


w = 0.12000000;
K = [207.84609691 0.00000000 160.00000000; 0.00000000 207.84609691 120.00000000; 0.00000000 0.00000000 1.00000000];
q = [135.74253750 150.92118362 155.00029328 136.48047414; 155.19908163 153.32591528 174.80770007 179.27370578];

q = [q; ones(1,4)];
p1 = [-0.5*w; -0.5*w];
p2 = [0.5*w; -0.5*w];
p3 = [0.5*w; 0.5*w];
p4 = [-0.5*w; 0.5*w];
p = [p1 p2 p3 p4; 0 0 0 0; 1 1 1 1];
PI_0 = [1 0 0 0; 0 1 0 0; 0 0 1 0];


[~,npoints]=size(q);
n = 2;
trials = 1; 



while trials < 51
rng shuffle
[R,~]=qr(randn(3));
posi = randn(3,1);
M = [R posi; 0 0 0 1];

tic

while n > 1e-10
  
   for c = 1:npoints
   [A{c},B{c},error{c}]=calc_b_eta(M,PI_0,K,q(:,c),p(:,c));
   end
   
   sum_A = sum(cat(3,A{:}),3);
   sum_B = sum(cat(3,B{:}),3);
   sum_error = sum(cat(3,error{:}),3);
   
   u = inv(sum_A+0.01*eye(6)) * sum_B;
    
   M = expm(ss61(u))*M;           
   n=sum_error;
   
   Ms{trials} = M;
    
   nu(trials) = norm(u);
  
   e = toc;
   
   if e > 0.5
       
       break
   end
end
disp(trials)
errors(trials) = sum_error;
trials = trials + 1;
end

[~,min_ind]=min(errors);
min_m = Ms{min_ind};
disp(mat2str(min_m))