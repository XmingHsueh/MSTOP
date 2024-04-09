function objs = comp_addition(H,G)

[n,m] = size(H);
objs = zeros(n,m);

for i = 1:n
    objs(i,:) = H(i,:)+G(i);
end