clc,clear

%% Specify a number of MSTOPs to be generated
comp_list = {'C1','C2'};
shape_list = {'H1','H2','H3','H4','H5'};
landscape_list = {'G1','G2','G3','G4','G5','G6','G7','G8'};
embed_list = {'E1','E2','E3'};
dim_list = [25 30 50];
m_list = [2 3];
trans_sce_list = {'Ta','Te'};
sim_location_list = {'H1','H2','M1','M2','M3','M4','L1','L2'};
sim_shape_list = {'E','I','H'};
k = 100;
specifications = [1 1 1 1 3 1 1 2 1 k; % MSTOP 1
    2 2 8 3 2 2 1 3 2 k; % MSTOP 2
    1 3 7 2 1 2 1 5 3 k; % MSTOP 3
    2 4 2 2 2 1 1 7 1 k; % MSTOP 4
    1 5 3 3 3 2 1 8 2 k; % MSTOP 5
    2 4 6 2 3 1 1 7 3 k; % MSTOP 6
    1 2 4 3 2 2 2 1 1 k; % MSTOP 7
    2 3 5 2 1 1 2 4 2 k; % MSTOP 8
    1 1 1 3 2 2 2 6 3 k; % MSTOP 9
    2 5 5 3 3 1 2 8 1 k; % MSTOP 10
    1 3 2 2 1 1 2 7 2 k; % MSTOP 11
    2 4 4 2 3 2 2 8 3 k]; % MSTOP 12
optimizer = 'NSGAII';
FEsMax = 1e4; % the maximum function evaluations
no_problems = size(specifications,1);

%% Generation of the STOPs
folder_problems = '.\benchmarks'; % the folder used for storing the generated STOPs
if ~isfolder(folder_problems)
    mkdir(folder_problems);
end
count = 0; % the number of available STOPs
for n = 1:size(specifications,1)
    MSTOP('comp',comp_list{specifications(n,1)},...
        'func_shape',shape_list{specifications(n,2)},...
        'func_landscape',landscape_list{specifications(n,3)},...
        'embed',embed_list{specifications(n,4)},...
        'dim',dim_list(specifications(n,5)),...
        'm',m_list(specifications(n,6)),...
        'trans_sce',trans_sce_list{specifications(n,7)},...
        'sim_location',sim_location_list{specifications(n,8)},...
        'sim_shape',sim_shape_list{specifications(n,9)},...
        'k',specifications(n,10),...
        'optimizer',optimizer,...
        'FEsMax',FEsMax,...
        'mode','geenration',...
        'folder_stops',folder_problems);
    count = count+1;
    fprintf('#%d of the 12 problems is ready!\n',count);
end
addpath(folder_problems);