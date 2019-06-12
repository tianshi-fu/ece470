function Z = matexpo(S,theta)
w = S(1:3);
v = S(4:6);
matexp_w_theta = eye(3) + sin(theta)*ss31(w)+(1-cos(theta))*(ss31(w))^2;
p = (eye(3)*theta + (1-cos(theta))*ss31(w) + (theta - sin(theta))*(ss31(w))^2 )*v;
Z = [matexp_w_theta p; 0 0 0 1];
end