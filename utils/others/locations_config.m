function [Xs,sims] = locations_config(xt,sim_distribution,k)

% unis = unifrnd(0,1,1,k);
unis = linspace(0,1,k); % k uniformly distributed real numbers within [0,1]
sims = zeros(1,k);

switch(sim_distribution)
    case 'H1'
        sims = ones(1,k);
    case 'H2'
        for i = 1:k
            sims(i) = (sqrt(unis(i))+1)/2;
        end
    case 'M1'
        for i = 1:k
            sims(i) = unis(i);
        end
    case 'M2'
        for i = 1:k
            sims(i) = sqrt(unis(i));
        end
    case 'M3'
        for i = 1:k
            sims(i) = 1-sqrt(1-unis(i));
        end
    case 'M4'
        for i = 1:k
            if unis(i)<=0.5
                sims(i) = sqrt(unis(i)/2);
            else
                sims(i) = (-sqrt(1-unis(i))+sqrt(2))/sqrt(2);
            end
        end
    case 'L1'
        sims = zeros(1,k);
    case 'L2'
        for i = 1:k
            sims(i) = (-sqrt(1-unis(i))+1)/2;
        end
end

dim = length(xt);
Xs = zeros(k,dim);
for i = 1:k
    r = rand(1,dim);
    x_gen = xt+(1-sims(i))/norm(r-xt,inf)*(r-xt);
    x_gen(x_gen>1) = 1;
    x_gen(x_gen<0) = 0;
    while sum(x_gen==0)+sum(x_gen==1)==dim
        r = rand(1,dim);
        x_gen = xt+(1-sims(i))/norm(r-xt,inf)*(r-xt);
        x_gen(x_gen>1) = 1;
        x_gen(x_gen<0) = 0;
    end
    Xs(i,:) = x_gen;
end