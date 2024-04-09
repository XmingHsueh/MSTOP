function objs = landscape_ellipsoid(Xc,opt,xdc_lb,xdc_ub,lb,ub)
%Ellipsoid function
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
    sum_e = 0;
    for j = nd+1:dim
        sum_e = sum_e + (dim-j+1)*(X(i,j)-x_opt(j))^2;
    end
    objs(i) = multiplier*sum_e;
end

end

