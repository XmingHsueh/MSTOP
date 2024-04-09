function objs = objs_evaluator(X,task,lb,ub)

[n,d] = size(X);
Xc = (X-repmat(lb,n,1))./(repmat(ub,n,1)-repmat(lb,n,1));
comp = task.func_comp;
embed_inv = task.func_embed_inv;
shape = task.func_shape;
landscape = task.func_landscape;

U = embed_inv(Xc);
Ud = U(:,1:task.m-1);
G = landscape(U);
H = shape(Ud);
objs = comp(H,G);