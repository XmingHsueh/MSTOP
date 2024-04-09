function objs = landscape_ackley(Xc,M,opt,xdc_lb,xdc_ub,lb,ub)
%Ackley function
%   - X: a set of design vectors
%   - M: rotation matrix
%   - opt: optimum
    [n,dim] = size(Xc);
    X = repmat(lb,n,1)+(repmat(ub,n,1)-repmat(lb,n,1)).*Xc;
    nd = length(xdc_lb);
    xd_lb = lb(1:nd)+(ub(1:nd)-lb(1:nd)).*xdc_lb;
    xd_ub = lb(1:nd)+(ub(1:nd)-lb(1:nd)).*xdc_ub;
    x_opt = opt.*(ub-lb)+lb;
    scale = ub(1)-lb(1);
    objs = zeros(n,1);
    for i = 1:n
        if sum(X(i,1:nd)>=xd_lb)+sum(X(i,1:nd)<=xd_ub) == 2*nd
            multiplier = 1;
        else
            multiplier = 1+sum(xd_lb(X(i,1:nd)<xd_lb)-X(i,X(i,1:nd)<xd_lb))/scale+...
                sum(X(i,X(i,1:nd)>xd_ub)-xd_ub(X(i,1:nd)>xd_ub))/scale;
        end
        var = (M(nd+1:end,nd+1:end)*(X(i,nd+1:end)-x_opt(nd+1:end))')';
        sum1 = 0; sum2 = 0;
        for j = 1: dim-nd
            sum1 = sum1 + var(j)*var(j);
            sum2 = sum2 + cos(2*pi*var(j));
        end
        avgsum1 = sum1/dim;
        avgsum2 = sum2/dim;
        objs(i) = multiplier*(-20*exp(-0.2*sqrt(avgsum1)) - exp(avgsum2) + 20 + exp(1));    
    end
end

