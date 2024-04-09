function problem = mop_config_source(comp,func_landscape,func_shape,location,sim_location,sim_shape,target,dim,m)

problem.dim = dim;
problem.m = m;
task = struct;
task.m = m;

switch(comp)
    case 'C1'
        task.func_comp = @(x,y)comp_multiplication(x,y);
    case 'C2'
        task.func_comp = @(x,y)comp_addition(x,y);
end

switch(func_landscape)
    case 'G1'
        problem.lb = -100*ones(1,dim);
        problem.ub = 100*ones(1,dim);
        u_opt = 0.5*ones(1,dim);
        task.ud_lb = target.ud_lb;
        task.ud_ub = target.ud_ub;
        task.func_landscape = @(u)landscape_sphere(u,u_opt,task.ud_lb,task.ud_ub,problem.lb,problem.ub);
    case 'G2'
        problem.lb = -50*ones(1,dim);
        problem.ub = 50*ones(1,dim);
        u_opt = 0.5*ones(1,dim);
        task.ud_lb = target.ud_lb;
        task.ud_ub = target.ud_ub;
        task.func_landscape = @(u)landscape_ellipsoid(u,u_opt,task.ud_lb,task.ud_ub,problem.lb,problem.ub);
    case 'G3'
        problem.lb = -30*ones(1,dim);
        problem.ub = 30*ones(1,dim);
        u_opt = 0.5*ones(1,dim);
        task.ud_lb = target.ud_lb;
        task.ud_ub = target.ud_ub;
        task.func_landscape = @(u)landscape_schwefel(u,u_opt,task.ud_lb,task.ud_ub,problem.lb,problem.ub);
    case 'G4'
        problem.lb = -5*ones(1,dim);
        problem.ub = 5*ones(1,dim);
        u_opt = 0.5*ones(1,dim);
        task.ud_lb = target.ud_lb;
        task.ud_ub = target.ud_ub;
        task.func_landscape = @(u)landscape_quartic(u,u_opt,task.ud_lb,task.ud_ub,problem.lb,problem.ub);
    case 'G5'
        problem.lb = -32*ones(1,dim);
        problem.ub = 32*ones(1,dim);
        u_opt = 0.5*ones(1,dim);
        task.ud_lb = target.ud_lb;
        task.ud_ub = target.ud_ub;
        rotation = eye(dim,dim);
        task.func_landscape = @(u)landscape_ackley(u,rotation,u_opt,task.ud_lb,task.ud_ub,problem.lb,problem.ub);
    case 'G6'
        problem.lb = -10*ones(1,dim);
        problem.ub = 10*ones(1,dim);
        u_opt = 0.5*ones(1,dim);
        task.ud_lb = target.ud_lb;
        task.ud_ub = target.ud_ub;
        rotation = eye(dim,dim);
        task.func_landscape = @(u)landscape_rastrigin(u,rotation,u_opt,task.ud_lb,task.ud_ub,problem.lb,problem.ub);
    case 'G7'
        problem.lb = -200*ones(1,dim);
        problem.ub = 200*ones(1,dim);
        u_opt = 0.5*ones(1,dim);
        task.ud_lb = target.ud_lb;
        task.ud_ub = target.ud_ub;
        rotation = eye(dim,dim);
        task.func_landscape = @(u)landscape_griewank(u,rotation,u_opt,task.ud_lb,task.ud_ub,problem.lb,problem.ub);
    case 'G8'
        problem.lb = -20*ones(1,dim);
        problem.ub = 20*ones(1,dim);
        u_opt = 0.5*ones(1,dim);
        task.ud_lb = target.ud_lb;
        task.ud_ub = target.ud_ub;
        task.func_landscape = @(u)landscape_levy(u,u_opt,task.ud_lb,task.ud_ub,problem.lb,problem.ub);
end

switch(func_shape)
    case 'H1'
        task.func_shape = @(x)shape_concave(x,task.ud_lb,task.ud_ub);
    case 'H2'
        task.func_shape = @(x)shape_concaveInv(x,task.ud_lb,task.ud_ub);
    case 'H3'
        task.func_shape = @(x)shape_convex(x,task.ud_lb,task.ud_ub);
    case 'H4'
        task.func_shape = @(x)shape_linear(x,task.ud_lb,task.ud_ub);
    case 'H5'
        task.func_shape = @(x)shape_linearInv(x,task.ud_lb,task.ud_ub);
end

no_points_ps = 100;
U = repmat(u_opt,no_points_ps,1);
U(2:no_points_ps,1:m-1) = lhsdesign_modified(no_points_ps-1,task.ud_lb,task.ud_ub);

switch(sim_shape)
    case 'E'
        shift = location-u_opt;
        task.func_embed_inv = @(x)embed_identity_source(x,shift,target.phi,target.phi_inv,'x2u');
        psx = embed_identity_source(U,shift,target.phi,target.phi_inv,'u2x');
        problem.location = embed_identity_source(u_opt,shift,target.phi,target.phi_inv,'u2x');
    case 'I'
        shift = location-u_opt;
        R = gramschmidt(rand(dim));
        task.func_embed_inv = @(x)embed_isometry_source(x,R,shift,u_opt,target.phi,target.phi_inv,'x2u');
        psx = embed_isometry_source(U,R,shift,u_opt,target.phi,target.phi_inv,'u2x');
        problem.location = embed_isometry_source(u_opt,R,shift,u_opt,target.phi,target.phi_inv,'u2x');
    case 'H'
        no_layers = 1;
        W = cell(no_layers,1);
        b = cell(no_layers,1);
        sigma = 5;
        for i = 1:no_layers
            W{i} = sigma*(rand(dim,dim)-0.5);
            b{i} = sigma*(rand(dim,1)-0.5);
        end
        shift = location-u_opt;
        task.func_embed_inv = @(x)embed_homeo_source(x,W,b,sigma,u_opt,shift,target.phi,target.phi_inv,'x2u');
        psx = embed_homeo_source(U,W,b,sigma,u_opt,shift,target.phi,target.phi_inv,'u2x');
        problem.location = embed_homeo_source(u_opt,W,b,sigma,u_opt,shift,target.phi,target.phi_inv,'u2x');
end

problem.func = @(x)objs_evaluator(x,task,problem.lb,problem.ub);
ps_valid = [];
for i =1:no_points_ps
    if sum(psx(i,:)<=1)==dim&&sum(psx(i,:)>=0)==dim
        ps_valid = [ps_valid;problem.lb+(problem.ub-problem.lb).*psx(i,:)];
    end
end
problem.ps = ps_valid;
problem.ud_lb = task.ud_lb;
problem.ud_ub = task.ud_ub;
problem.sim_location = sim_location;
problem.sim_shape = sim_shape;