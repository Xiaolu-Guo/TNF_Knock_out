% TNFo_dual_para_sample_main
function [] = TNFo_dual_para_p100o(vers,data_save_file_path,Num_sample)

data_info.save_file_path = data_save_file_path;

max_fold_sig = 2;
mean_fold = 1.25;

% the paramters that will be sampled
gene_info.parameter_name_vec = {{'params54','params61'}};
gene_info.para_val_vec = {[8.224e-06, 0.125]};
gene_info.gene_type = {'wt','tnfo'};

sigma_cor_pos = [1,0.8;0.8,1]* (log(max_fold_sig)/3);
sigma_cor_neg = [1,-0.8;-0.8,1]* (log(max_fold_sig)/3);
gene_info.para_val_vec_corr = {sigma_cor_neg,sigma_cor_pos};

tnf_fold_pos = [mean_fold,mean_fold];
tnf_fold_neg = [mean_fold,1/mean_fold];
gene_info.para_val_vec_tnf_fold = {tnf_fold_neg};

% sampling ratio (lognormal distribution


gene_info.gene_parameter_value_vec_genotype = cell(0);

for i_para_name = 1:length(gene_info.parameter_name_vec)
    mu_vec = gene_info.para_val_vec{i_para_name};
    mu_vec = log(mu_vec);
    
    sig_cor = gene_info.para_val_vec_corr{i_para_name};
    R_wt = mvnrnd(mu_vec,sig_cor,Num_sample);
    R_wt = exp(R_wt)';
    
    gene_info.gene_parameter_value_vec_genotype{1}{i_para_name} = R_wt;
    
    mu_vec = gene_info.para_val_vec{i_para_name}.*gene_info.para_val_vec_tnf_fold{i_para_name};
    mu_vec = log(mu_vec);
    sig_cor = gene_info.para_val_vec_corr{i_para_name};
    R_tnfo = mvnrnd(mu_vec,sig_cor,Num_sample);
    R_tnfo = exp(R_tnfo)';
    
    gene_info.gene_parameter_value_vec_genotype{2}{i_para_name} = R_tnfo;
end



% stimili info
sim_info.ligand = {'TNF'};
sim_info.dose_str = {'10ng/mL'};
sim_info.dose_val = {10};
sim_info.species_vec.name = {'NFkB'};
sim_info.species_vec.val = {0.08*10^0.1};

% species that will be saved
% must be r x 1, for each cell i must be ri x 1
data_info.species_outputname = {'nucNFkB';'TNFR';'IKK'};
data_info.species_composition = {{'NFkBn';'IkBaNFkBn'};{'TNFR'};{'IKK'}};
data_info.save_file_name = strcat('tnfo_dual_para_sample_',vers); % only beginning, no .mat

genotype_sim_save(sim_info,data_info,gene_info);




% current_folder = pwd;
% cd(current_folder)

