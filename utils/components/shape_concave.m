function objs = shape_concave(Xc,ud_lb,ud_ub)

[n,p_dim] = size(Xc);
m = p_dim+1;
objs = zeros(n,m);
Xd = (Xc-repmat(ud_lb,n,1))./(repmat(ud_ub,n,1)-repmat(ud_lb,n,1));
Xd(Xd<0) = 0;
Xd(Xd>1) = 1;

for i = 1:n
    if m == 2
        objs(i,1) = cos(0.5*pi*Xd(i));
        objs(i,2) = sin(0.5*pi*Xd(i));
    elseif m == 3
        objs(i,1) = cos(0.5*pi*Xd(i,1))*cos(0.5*pi*Xd(i,2));
        objs(i,2) = cos(0.5*pi*Xd(i,1))*sin(0.5*pi*Xd(i,2));
        objs(i,3) = sin(0.5*pi*Xd(i,1));
    else
        error('m is greater than three!');
    end
end