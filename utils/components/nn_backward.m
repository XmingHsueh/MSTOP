function x_o = nn_backward(x,W,b,sigma)

no_layers = length(b);
eps = 1e-2;
x(x<0) = eps;
x(x>1) = 1-eps;
x_out = x;
for i =no_layers:-1:1
    if i == no_layers
        x_in = inv(W{i})*(-log(1./(x_out)-1)-b{i});
    else
        x_in = W{i}\(-log(1./(x_out/sigma+0.5)-1)-b{i});
    end
    x_out = x_in;
end
x_o = x_in;