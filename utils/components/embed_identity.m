function Y = embed_identity(X,shift,direction)

[n,d] = size(X);
Y = zeros(n,d);
for i = 1:n
    switch (direction)
        case 'x2u'
            Y(i,:) = X(i,:)-shift;
        case 'u2x'
            Y(i,:) = X(i,:)+shift;
    end
end