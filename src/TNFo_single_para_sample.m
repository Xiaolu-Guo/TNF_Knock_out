% TNFo_single_para_sample_main
function [] = TNFo_single_para_sample(vers,data_save_file_path,Num_sample)
data_info.save_file_path = data_save_file_path;

% the paramters that will be sampled
gene_info.parameter_name_vec = {{'params54'},{'params61'},{'params62'},{'params63'},{'params65'}};
gene_info.para_val_vec = {8.224e-06, 0.125,1.875,320.3,1889};
gene_info.gene_type = {'wt','tnfo'};

% sampling ratio (lognormal distribution


rand_mat = randn(1,Num_sample);
fold_change_wt_log =log(10^0) + rand_mat* (log(2)/3);
fold_change_wt_vec_ele = exp(fold_change_wt_log);
fold_change_neg_tnfo_log = log(10^-0.1) + rand_mat* (log(2)/3);
fold_change_neg_tnfo_vec = exp(fold_change_neg_tnfo_log);
fold_change_pos_tnfo_log = log(10^0.1) + rand_mat* (log(2)/3);
fold_change_pos_tnfo_vec = exp(fold_change_pos_tnfo_log);

gene_info.gene_fold_change_wt_vec = cell(0);
gene_info.gene_fold_change_vec_genotype{1} = {fold_change_wt_vec_ele;
    fold_change_wt_vec_ele;
    fold_change_wt_vec_ele;
    fold_change_wt_vec_ele;
    fold_change_wt_vec_ele};
gene_info.gene_fold_change_vec_genotype{2} = {fold_change_pos_tnfo_vec;
    fold_change_neg_tnfo_vec;
    fold_change_pos_tnfo_vec;
    fold_change_neg_tnfo_vec;
    fold_change_pos_tnfo_vec};



% stimili info
sim_info.ligand = {'TNF'};
sim_info.dose_str = {'10ng/mL'};
sim_info.dose_val = {10};

% species that will be saved
% must be r x 1, for each cell i must be ri x 1
data_info.species_outputname = {'nucNFkB';'TNFR';'IKK'};
data_info.species_composition = {{'NFkBn';'IkBaNFkBn'};{'TNFR'};{'IKK'}};
data_info.save_file_name = strcat('tnfo_single_para_sample_',vers); % only beginning, no .mat

genotype_sim_save(sim_info,data_info,gene_info);

