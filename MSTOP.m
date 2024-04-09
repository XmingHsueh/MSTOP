% Author: Xiaoming Xue
% Email: xminghsueh@gmail.com
%
% ------------
% Description:
% ------------
% The class of black-box multiobjective sequential transfer optimization problems (MSTOPs).
%
% Configurable information of the target task:
% a) the composition method, C1-multiplication, C2-addition
% b) the shape function, H1-concave, H2-inverted concave, H3-convex, H4-linear, H5-inverted linear
% c) the landscape function, G1-Sphere, G2-Ellipsoid, G3-Schwefel 2.2,
% G4-Quartic, G5-Ackley, G6-Rastrigin, G7-Griewank, G8-Levy
% d) the embedding of Pareto manifold, E1-Identity, E2-Isometry, E3-Homeomorphism
% e) the problem dimension: a positive integer
% f) the number of objectives: 2 or 3
%
% Configurable information of the knowledge base:
% a) the number of source tasks, a positive integer
% b) the transfer scenario, Ta-intra-family transfer, Te-inter-family transfer
% c) the location similarity, HS-high similarity, MS-mixed similarity, LS-low similarity
% d) the shape similarity, E-Identity, I-Isometry, H-Homeomorphism
%
% ------------
% Reference:
% ------------
% X. Xue, L. Feng, C. Yang, S. Liu, L. Song, and K. C. Tan. "Multiobjective Sequential Transfer 
% Optimization: Benchmark Problems and Results", WCCI 2024.

classdef MSTOP

    % MSTOP parameters:
    % comp--->the composition, {C1, C2}
    % func_shape--->the shape function, {H1 to H5}
    % func_landscape--->the landscape function, {G1 to G8}
    % embed--->the embedding, {E1 to E3}
    % dim--->the problem dimension, a positive integer
    % m--->the number of objectives, 2 or 3
    % k--->the number of source tasks, a positive integer
    % trans_sce--->the transfer scenario, {Ta, Te}
    % sim_location--->the location similarity, {H1, H2, M1 to M4, L1, L2}
    % sim_shape--->the shape similarity, {E, I, H}
    %
    % optimizer--->the backbone optimizer, a valid search engine
    % FEsMax---the maximum function evaluations (FEs) for solving the source tasks, a positive integer
    % mode--->the mode of problem call, {gen: problem generation, opt: optimization}
    % target_problem--->the instantiated target task
    % source_problems---<read-only>the instantiated source tasks
    % knowledge_base--->the knowledge base with the solution data of the k source tasks
    % folder_stops--->the folder used for storing the generated MSTOPs
    % state_knowledgebase---<read-only>the availability of the specified STOP: 1->available; 0->unavailable

    properties
        comp = 'C1';
        func_shape = 'H1';
        func_landscape = 'G1';
        embed = 'E1';
        dim = 50;
        m = 2;
        k = 100;
        trans_sce = 'Ta'
        sim_location = 'H1';
        sim_shape = 'E';
        optimizer = 'NSGAII';
        FEsMax = 1e3;
        mode = 'gen';
        target_problem;
        knowledge_base = struct;
        folder_stops = 'benchmarks';
    end

    properties(SetAccess = protected)
        source_problems;
        state_knowledgebase;
    end

    methods

        function obj = MSTOP(varargin) % initialization
            isStr = find(cellfun(@ischar,varargin(1:end-1))&~cellfun(@isempty,varargin(2:end)));
            for i = isStr(ismember(varargin(isStr),{'comp','func_shape','func_landscape','embed',...
                    'dim','m','k','trans_sce','sim_location','sim_shape','optimizer','FEsMax','mode',...
                    'folder_stops'}))
                obj.(varargin{i}) = varargin{i+1};
            end

            % examine the availability of the specified MSTOP
            dir_sesto = [obj.folder_stops,'\',obj.comp,'-',obj.func_shape,'-',obj.func_landscape,...
                '-',obj.embed,'-d',num2str(obj.dim),'-m',num2str(obj.m),'-',obj.trans_sce,...
                '-',obj.sim_location,'-',obj.sim_shape,'-k',num2str(obj.k),'.mat'];
            obj.state_knowledgebase = sign(exist(dir_sesto,'file'));

            if obj.state_knowledgebase == 1 && strcmp(obj.mode,'opt') % load the MSTOP
                load([obj.folder_stops,'\',obj.comp,'-',obj.func_shape,'-',obj.func_landscape,...
                '-',obj.embed,'-d',num2str(obj.dim),'-m',num2str(obj.m),'-',obj.trans_sce,...
                '-',obj.sim_location,'-',obj.sim_shape,'-k',num2str(obj.k),'.mat']);
                obj.target_problem = target;
                obj.source_problems = sources;
                for i = 1:obj.k
                    obj.knowledge_base(i).solutions = knowledge(i).solutions;
                    obj.knowledge_base(i).objs = knowledge(i).objs;
                end
            elseif obj.state_knowledgebase == 0 % generate the MSTOP when it is not found
                obj = obj.Generation(); 
            end
        end

        function obj = Generation(obj) % problem constructor
            target = mop_config_target(obj.comp,obj.func_shape,obj.func_landscape,obj.embed,obj.dim,obj.m); % configure the target task
            obj.target_problem = target;
            [locations_source,sims_location] = locations_config(target.location,obj.sim_location,obj.k); % generate k location-related similarity values
            for i = 1:obj.k % configure the source tasks
                switch(obj.trans_sce)
                    case 'Ta' % intra-family transfer
                        comp_s = obj.comp;
                        func_landscape_s = obj.func_landscape;
                        func_shape_s = obj.func_shape;
                    case 'Te' % inter-family transfer
                        r = rand;
                        if r<0.5
                            comp_s = 'C1';
                        else
                            comp_s = 'C2';
                        end
                        r = randi(8);
                        while ceil(str2double(obj.func_landscape(2)))==r
                            r = randi(8);
                        end
                        func_landscape_s = ['G',num2str(r)];
                        r = randi(5);
                        while ceil(str2double(obj.func_shape(2)))==r
                            r = randi(5);
                        end
                        func_shape_s = ['H',num2str(r)];
                end
                sources(i) = mop_config_source(comp_s,func_landscape_s,func_shape_s,...
                    locations_source(i,:),sims_location(i),obj.sim_shape,target,obj.dim,obj.m);
            end
            h=waitbar(0,'Starting');
            for i = 1:obj.k % optimize the k source tasks
                mop.fun = sources(i).func;
                mop.dim = sources(i).dim;
                mop.lb = sources(i).lb;
                mop.ub = sources(i).ub;
                cmd_opt = ['[solutions,objs] = ',obj.optimizer,'(mop,obj.FEsMax,i);'];
                eval(cmd_opt);
                knowledge(i).solutions = solutions;
                knowledge(i).objs = objs;
                waitbar(i/obj.k,h,sprintf('Generation in progress: %.2f%%',i/obj.k*100));
            end
            close(h);
            save([obj.folder_stops,'\',obj.comp,'-',obj.func_shape,'-',obj.func_landscape,...
                '-',obj.embed,'-d',num2str(obj.dim),'-m',num2str(obj.m),'-',obj.trans_sce,...
                '-',obj.sim_location,'-',obj.sim_shape,'-k',num2str(obj.k),'.mat'],'target','sources','knowledge');
        end
    end
end