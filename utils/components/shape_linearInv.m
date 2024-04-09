function objs = shape_linearInv(X,ud_lb,ud_ub)

[n,p_dim] = size(X);
m = p_dim+1;
objs = zeros(n,m);
Xd = (X-repmat(ud_lb,n,1))./(repmat(ud_ub,n,1)-repmat(ud_lb,n,1));
Xd(Xd<0) = 0;
Xd(Xd>1) = 1;

for i = 1:n
    if m == 2
        objs(i,1) = 1-Xd(i);
        objs(i,2) = Xd(i);
    elseif m == 3
        objs(i,1) = 1-Xd(i,1)*Xd(i,2);
        objs(i,2) = 1-Xd(i,1)*(1-Xd(i,2));
        objs(i,3) = Xd(i,1);
    else
        error('m is greater than three!');
    end
end