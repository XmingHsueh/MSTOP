function objs = comp_multiplication(H,G)

[n,m] = size(H);
objs = zeros(n,m);

for i = 1:n
    objs(i,:) = H(i,:)*(1+G(i));
end