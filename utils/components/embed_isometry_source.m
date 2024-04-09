function Y = embed_isometry_source(X,R,shift,u_opt,embed_t,embed_t_inv,direction)

[n,d] = size(X);
Y = zeros(n,d);
for i = 1:n
    switch (direction)
        case 'x2u'
            Y(i,:) = embed_t_inv((X(i,:)-shift-u_opt)*inv(R)+u_opt);
        case 'u2x'
            Y(i,:) = (embed_t(X(i,:))-u_opt)*R+u_opt+shift;
    end
end