function Y = embed_homeo(X,W,b,sigma,u_opt,shift,direction)

[n,d] = size(X);
Y = zeros(n,d);
x_center = transpose(nn_feedforward(zeros(d,1),W,b,sigma));
for i = 1:n
    switch(direction)
        case 'x2u'
            Y(i,:) = transpose(nn_backward(X(i,:)'-shift'-u_opt'+x_center',W,b,sigma))+u_opt;
        case 'u2x'
            Y(i,:) = transpose(nn_feedforward(X(i,:)'-u_opt',W,b,sigma))-x_center+u_opt+shift;
    end
end