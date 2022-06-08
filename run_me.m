% run_me.m for TNF knockout proj

run_TNF_sens = 0;
draw_TNF_sens = 0;
run_sample_TNFo = 0;
draw_sample_TNFo = 1;
run_sample_p100o = 0;
draw_sample_p100o = 0;

addpath('./lib/')
addpath('./src/')


%% Scanning paramter in TNF receptor module
data_save_file_path = '../raw_data/';
fig_save_path = '../../TNF_Knock_out/Figures/';

if run_TNF_sens
    Module_parameters_sim(data_save_file_path)
end
if draw_TNF_sens
    draw_Sens_TNF(data_save_file_path,fig_save_path)
end

%% Sampling different paramter (combinations) for TNF-/- and WT
data_save_file_path = '../raw_data/';
vers = '20220606';
Num_sample = 2000;

if run_sample_TNFo
    
    TNFo_dual_para_sample_main(vers,data_save_file_path,Num_sample)
    TNFo_single_para_sample_main(vers,data_save_file_path,Num_sample)
end

if draw_sample_TNFo
    fig_save_path = '../../TNF_Knock_out/Figures/';
    vers_fig = 'v1';
   
    single_or_dual = 'dual';
    draw_TNFo(vers,data_save_file_path, single_or_dual,vers_fig, fig_save_path )
    
    single_or_dual = 'single';
    draw_TNFo(vers,data_save_file_path, single_or_dual,vers_fig, fig_save_path )
end

%% Sampling for TNF-/- p100o and TNF-/-
%
if run_sample_p100o
    vers = 'p100o_20220606';
    TNFo_dual_para_p100o_main(vers,data_save_file_path,Num_sample)
end
