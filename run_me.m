% run_me.m for TNF knockout proj
run_TNF_sens = 1;
draw_TNF_sens = 1;
run_sample_TNFo = 1;
draw_sample_TNFo = 1;
run_sample_p100o = 1;
draw_sample_p100o = 1;
 

%% user could change
vers = '20220606';
Num_sample = 10;
data_save_file_path = './raw_data/';
fig_save_path = './Figures/';
% data_save_file_path = '../raw_data/';
% fig_save_path = '../../TNF_Knock_out/Figures/';


%% initialize
if ~isfolder(data_save_file_path)
    mkdir(data_save_file_path)
end

if ~isfolder(fig_save_path)
    mkdir(fig_save_path)
end

addpath('./lib/')
addpath('./src/')


%% Scanning paramter in TNF receptor module
if run_TNF_sens
    Module_parameters_sim(data_save_file_path)
end

if draw_TNF_sens
    draw_Sens_TNF(data_save_file_path,fig_save_path)
end


%% Sampling different paramter (combinations) for TNF-/- and WT
if run_sample_TNFo
    
    TNFo_dual_para_sample(vers,data_save_file_path,Num_sample)
    TNFo_single_para_sample(vers,data_save_file_path,Num_sample)
end

if draw_sample_TNFo
    vers_fig = 'v1';
   
    single_or_dual = 'dual';
    draw_TNFo(vers,data_save_file_path, single_or_dual,vers_fig, fig_save_path )
    
    single_or_dual = 'single';
    draw_TNFo(vers,data_save_file_path, single_or_dual,vers_fig, fig_save_path )
end


%% Sampling for TNF-/- p100o and TNF-/-
if run_sample_p100o
    vers ='tnfo_20220606';% 
    NFkB_fold = 1;
    TNFo_dual_para_p100o(vers,data_save_file_path,Num_sample,NFkB_fold)

    vers = 'p100o_20220606';
    NFkB_fold = 1.25;
    TNFo_dual_para_p100o(vers,data_save_file_path,Num_sample,NFkB_fold)

end

if draw_sample_p100o
    vers_TNFo =  'tnfo_20220606';
    
    vers_fig = 'v1';
    vers_p100oTNFo = 'p100o_20220606';
    draw_p100o(vers_TNFo,vers_p100oTNFo,data_save_file_path,vers_fig, fig_save_path )
    
end