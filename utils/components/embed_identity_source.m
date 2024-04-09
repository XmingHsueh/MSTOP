function Y = embed_identity_source(X,shift,embed_t,embed_t_inv,direction)

[n,d] = size(X);
Y = zeros(n,d);
for i = 1:n
    switch (direction)
        case 'x2u'
            Y(i,:) = embed_t_inv(X(i,:)-shift);
        case 'u2x'
            Y(i,:) = embed_t(X(i,:))+shift;
    end
end