function objs = landscape_levy(Xc,opt,xdc_lb,xdc_ub,lb,ub)
%LEVY function
%   - X: a set of design vectors
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
    w = zeros(1,dim-nd);
    for j = nd+1:dim
        w(j-nd) = 1 + (X(i,j) - x_opt(j))/4;
    end
    term1 = (sin(pi*w(1)))^2;
    term3 = (w(dim-nd)-1)^2 * (1+(sin(2*pi*w(dim-nd)))^2);
    term2 = 0;
    for j = 1:dim-nd-1
        wi = w(j);
        term2 = term2 + (wi-1)^2 * (1+10*(sin(pi*wi+1))^2);
    end
    objs(i) = multiplier*(term1 + term2 + term3);
end
end
