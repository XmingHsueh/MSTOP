function Y = embed_isometry(X,R,shift,u_opt,direction)

[n,d] = size(X);
Y = zeros(n,d);
for i = 1:n
    switch(direction)
        case 'x2u'
            Y(i,:) = (X(i,:)-shift-u_opt)*inv(R)+u_opt;
        case 'u2x'
            Y(i,:) = (X(i,:)-u_opt)*R+u_opt+shift;
    end
end