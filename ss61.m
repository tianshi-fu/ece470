function B = ss61(V)
w = V(1:3);
v = V(4:6);
ssw = ss31(w);
B = [ssw v; 0 0 0 0];
end