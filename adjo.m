function adj = adjo(M)
R = M(1:3,1:3);
p = M(1:3,4);
adj = [R zeros(3); ss31(p)*R R];
end
