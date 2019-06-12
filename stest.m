function A = stest(S,theta)
w = S(1:3);
v= S(4:6);
bw = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
ebwt = eye(3)+sin(theta)*bw+(1-cos(theta))*bw^2;
p = (eye(3)*theta+(1-cos(theta))*bw+(theta-sin(theta))*bw^2)*v;

A = [ebwt p; 0 0 0 1];


end