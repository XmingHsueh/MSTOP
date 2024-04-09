function objs = shape_concaveInv(X,ud_lb,ud_ub)

[n,p_dim] = size(X);
m = p_dim+1;
objs = zeros(n,m);
Xd = (X-repmat(ud_lb,n,1))./(repmat(ud_ub,n,1)-repmat(ud_lb,n,1));
Xd(Xd<0) = 0;
Xd(Xd>1) = 1;

for i = 1:n
    if m == 2
        objs(i,1) = 1-cos(0.5*pi*Xd(i));
        objs(i,2) = 1-sin(0.5*pi*Xd(i));
    elseif m == 3
        objs(i,1) = 1-cos(0.5*pi*Xd(i,1))*cos(0.5*pi*Xd(i,2));
        objs(i,2) = 1-cos(0.5*pi*Xd(i,1))*sin(0.5*pi*Xd(i,2));
        objs(i,3) = 1-sin(0.5*pi*Xd(i,1));
    else
        error('m is greater than three!');
    end
end