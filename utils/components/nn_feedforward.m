function x_new = nn_feedforward(x,W,b,sigma)

no_layers = length(b);
x_in = x;
for i =1:no_layers
    if i<no_layers
        x_out = sigma*(1./(1+exp(-(W{i}*x_in+b{i})))-0.5);
    else
        x_out = 1./(1+exp(-(W{i}*x_in+b{i})));
    end
    x_in = x_out;
end
x_new = x_out;