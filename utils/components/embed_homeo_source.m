function Y = embed_homeo_source(X,W,b,sigma,u_opt,shift,embed_t,embed_t_inv,direction)

[n,d] = size(X);
Y = zeros(n,d);
x_center = transpose(nn_feedforward(zeros(d,1),W,b,sigma));
for i = 1:n
    switch (direction)
        case 'x2u'
            Y(i,:) = embed_t_inv(transpose(nn_backward(X(i,:)'-shift'-u_opt'+x_center',W,b,sigma))+u_opt);
        case 'u2x'
            Y(i,:) = transpose(nn_feedforward(embed_t(X(i,:))'-u_opt',W,b,sigma))-x_center+u_opt+shift;
    end
end