
function [] = draw_TNFo(vers,data_save_file_path, single_or_dual,vers_fig, fig_save_path )
data_info.save_file_path = data_save_file_path;

colors.WT = [93,147,191]/255;
colors.TNFKO = [233,72,73]/255;

switch single_or_dual
    case 'dual'
        data_info.save_file_name = strcat('tnfo_dual_para_sample_',vers); % only beginning, no .mat
        gene_info.parameter_name_vec = {{'params54','params61'},{'params54','params62'},{'params54','params63'},{'params54','params65'},...
            {'params61','params62'},{'params61','params63'},{'params61','params65'},...
            {'params62','params63'},{'params62','params65'},{'params63','params65'}};
        fig_name = 'tnfo_dual_para_sample_';
        fig_num_distri =[1,2,3];
    case 'single'
        data_info.save_file_name = strcat('tnfo_single_para_sample_',vers); % only beginning, no .mat
        gene_info.parameter_name_vec = {{'params54'},{'params61'},{'params62'},{'params63'},{'params65'}};
        fig_name = 'tnfo_single_para_sample_';
        fig_num_distri =[1];
end


genot_vec = {'wt','tnfo'};


for i_para_name = 1:length(gene_info.parameter_name_vec)
    
    %% load sim data
    sim_info.parameter_name = gene_info.parameter_name_vec{i_para_name};
    
    if isfield(data_info,'save_file_path') && isfield(data_info,'save_file_name')
        save_file_name = data_info.save_file_name;
        fig_save_name = strcat(fig_name,vers_fig);
        for i_para_name = 1:length(sim_info.parameter_name)
            save_file_name = strcat(save_file_name,'_',sim_info.parameter_name{i_para_name});
            fig_save_name = strcat(fig_save_name,'_',sim_info.parameter_name{i_para_name});
        end
        % save_file_name = strcat(save_file_name,'.mat');
        load(strcat(data_info.save_file_path,save_file_name,'.mat'),'sim_data_tbl','metrics','data');
        
    end
    
    %% get the data for drawing/plotting
    index_nuc_NFkB_wt = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='nucNFkB'...
        & sim_data_tbl.type == 'wt';
    
    index_nuc_NFkB_tnfo = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='nucNFkB'...
        & sim_data_tbl.type == 'tnfo';
    
    index_TNFR_wt = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='TNFR'...
        & sim_data_tbl.type == 'wt';
    
    index_TNFR_tnfo = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='TNFR'...
        & sim_data_tbl.type == 'tnfo';
    
    index_IKK_wt = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='IKK'...
        & sim_data_tbl.type == 'wt';
    
    index_IKK_tnfo = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='IKK'...
        & sim_data_tbl.type == 'tnfo';
    
    IKK_wt = mean(sim_data_tbl.trajectory(index_IKK_wt,:));
    IKK_tnfo = mean(sim_data_tbl.trajectory(index_IKK_tnfo,:));
    
    TNFR_wt = mean(sim_data_tbl.trajectory(index_TNFR_wt,1));
    TNFR_tnfo = mean(sim_data_tbl.trajectory(index_TNFR_tnfo,1));
    
    traj_wt = sim_data_tbl.trajectory(index_nuc_NFkB_wt,:);
    [~,traj_order] = sort(max(traj_wt,[],2),'descend');
    traj_wt = traj_wt(traj_order,:);
    traj_wt = traj_wt(:,1:5:end);
    traj_tnfo = sim_data_tbl.trajectory(index_nuc_NFkB_tnfo,:);
    [~,traj_order] = sort(max(traj_tnfo,[],2),'descend');
    traj_tnfo = traj_tnfo(traj_order,:);
    traj_tnfo = traj_tnfo(:,1:5:end);
    
    
    parameters_wt = sim_data_tbl.parameter_value(index_nuc_NFkB_wt,:);
    parameters_tnfo = sim_data_tbl.parameter_value(index_nuc_NFkB_tnfo,:);
    
    
    %% draw the picture and save
    
    if 0 %
        draw_osc_ratio_bar([1,2,3],colors,metrics,data)
        
        figure(1)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_heatmap_byospower'),'epsc')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_heatmap_byospower'),'svg')
        close
        
        figure(2)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_oscpercentage_diff_thresh'),'epsc')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_oscpercentage_diff_thresh'),'svg')
        close
        %
        figure(3)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_oscpercentage'),'epsc')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_oscpercentage'),'svg')
        close
        
    end
    
    if 1 %
        draw_traj_heatmap([1,2],traj_wt,traj_tnfo)
        
        figure(1)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_trajheatmap_wt'),'epsc')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_trajheatmap_wt'),'svg')
        close
        
        figure(2)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_trajheatmap_tnfo'),'epsc')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_trajheatmap_tnfo'),'svg')
        close
    end
    
    if 0 %parameter distribution
        
        draw_distribution(fig_num_distri ,colors, parameters_wt,parameters_tnfo,sim_info)
        figure(1)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_para1_distri'),'epsc')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_para1_distri'),'svg')
        close
        
        if length(fig_num_distri)>1
            figure(2)
            saveas(gcf,strcat(fig_save_path,fig_save_name,'_para2_distri'),'epsc')
            saveas(gcf,strcat(fig_save_path,fig_save_name,'_para2_distri'),'svg')
            close
            
            figure( 3)
            saveas(gcf,strcat(fig_save_path,fig_save_name,'_para_distri_2d'),'epsc')
            saveas(gcf,strcat(fig_save_path,fig_save_name,'_para_distri_2d'),'svg')
            close
        end
        
    end
    fig_num =1 ;
    if 0 %IKK
        
        
        
        draw_IKK(fig_num,colors,IKK_wt,IKK_tnfo)
        figure(fig_num)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_IKK'),'epsc')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_IKK'),'svg')
        close
    end
    
    if 0 %TNFR
        draw_TNFR(fig_num,colors,TNFR_wt,TNFR_tnfo)
        figure(fig_num)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_tnfr'),'epsc')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_tnfr'),'svg')
        close
        
    end
    
end
end

%% draw oscllatory cells

function [] = draw_osc_ratio_bar(fig_num,colors,metrics,data)

threshold = [0.5e-4, 1e-4, 1.5e-4, 2e-4, 2.5e-4];%[0.5e-4,1e-4,1.5e-4];(For Xiaolu's oscpower)
osc_ratio = NaN(2,length(threshold));
for i_genotype =1:2
    for i_threshold =1:length(threshold)
        osc_ratio(i_genotype,i_threshold) = sum(metrics{i_genotype}.oscpower > threshold(i_threshold) )/size(data.model_sim{i_genotype},1);
    end
    
end

[~,traj_order] = sort(metrics{1}.oscpower,'descend');
traj_wt = data.model_sim{1}(traj_order,:);

[~,traj_order] = sort(metrics{2}.oscpower,'descend');
traj_tnfo = data.model_sim{2}(traj_order,:);

if 0
    figure(fig_num(1))
    subplot(1,2,1)
    h1=heatmap(traj_wt,'ColorMap',parula,'GridVisible','off','XLabel',{},'YLabel',{},'ColorLimits',[-0.001,0.25]);
    XLabels = 0:1:size(traj_wt,2)-1;
    % Convert each number in the array into a string
    CustomXLabels = string(XLabels/12);
    % Replace all but the fifth elements by spaces
    CustomXLabels(:) = " ";
    % Set the 'XDisplayLabels' property of the heatmap
    % object 'h' to the custom x-axis tick labels
    h1.XDisplayLabels = CustomXLabels;
    
    YLabels = 1:size(traj_wt,1);
    % Convert each number in the array into a string
    YCustomXLabels = string(YLabels);
    % Replace all but the fifth elements by spaces
    YCustomXLabels(:) = " ";
    % Set the 'XDisplayLabels' property of the heatmap
    % object 'h' to the custom x-axis tick labels
    h1.YDisplayLabels = YCustomXLabels;
    colorbar('off')
    
    subplot(1,2,2)
    h2=heatmap(traj_tnfo,'ColorMap',parula,'GridVisible','off','XLabel',{},'YLabel',{},'ColorLimits',[-0.001,0.25]);
    
    XLabels = 0:1:size(traj_wt,2)-1;
    % Convert each number in the array into a string
    CustomXLabels = string(XLabels/12);
    % Replace all but the fifth elements by spaces
    CustomXLabels(:) = " ";
    % Set the 'XDisplayLabels' property of the heatmap
    % object 'h' to the custom x-axis tick labels
    h2.XDisplayLabels = CustomXLabels;
    
    YLabels = 1:size(traj_wt,1);
    % Convert each number in the array into a string
    YCustomXLabels = string(YLabels);
    % Replace all but the fifth elements by spaces
    YCustomXLabels(:) = " ";
    % Set the 'XDisplayLabels' property of the heatmap
    % object 'h' to the custom x-axis tick labels
    h2.YDisplayLabels = YCustomXLabels;
    colorbar('off')
    
    Set_figure_size_wide
    
end


figure(fig_num(2))

for i_osratio = 1:size(osc_ratio,2)
    
    subplot(1,length(threshold),i_osratio)
    
    b=bar(osc_ratio(:,i_osratio)*100);
    %xlabel()
    b.FaceColor = 'flat';
    b.CData(1,:) = colors.WT;
    b.CData(2,:) = colors.TNFKO;
    
    ylim([0,100])
    
    ytickformat(gca, '%g%%');
    yt = yticks;
    ytlable = cell(1,length(yt));
    for i_yt = 1:length(yt)
        ytlable{i_yt} = '';
    end
    ax = gca;
    ax.YGrid = 'on';
    %     ylabel('Osc ratio')
    set(gca,'YTickLabel',ytlable,'XTickLabel',{'',''})
    %     set(gca,'fontsize',7,'fontweight','b')
    
    
end
Set_figure_size_square


figure(fig_num(3))

i_osratio = ceil(size(osc_ratio,2)/2);
b=bar(osc_ratio(:,i_osratio)*100);
%xlabel()
b.FaceColor = 'flat';
b.CData(1,:) = colors.WT;
b.CData(2,:) = colors.TNFKO;

ylim([0,100])

ytickformat(gca, '%g%%');
yt = yticks;
ytlable = cell(1,length(yt));
for i_yt = 1:length(yt)
    ytlable{i_yt} = '';
end
ax = gca;
ax.YGrid = 'on';
set(gca,'YTickLabel',ytlable,'XTickLabel',{'',''})
Set_figure_size_square
end

%% draw heatmap
% to do: delete ylabel, xlabel only time 0,4,8
% proper figure size
% seperate the figures
% draw heatmap
function [] = draw_traj_heatmap(fig_num,traj_wt,traj_tnfo)

figure(fig_num(1))
h1=heatmap(traj_wt,'ColorMap',parula,'GridVisible','off','XLabel',{},'YLabel',{},'ColorLimits',[-0.001,0.25]);

XLabels = 0:1:size(traj_wt,2)-1;
% Convert each number in the array into a string
CustomXLabels = string(XLabels/12);
% Replace all but the fifth elements by spaces
CustomXLabels(:) = " ";
% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h1.XDisplayLabels = CustomXLabels;

YLabels = 1:size(traj_wt,1);
% Convert each number in the array into a string
YCustomXLabels = string(YLabels);
% Replace all but the fifth elements by spaces
YCustomXLabels(:) = " ";
% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h1.YDisplayLabels = YCustomXLabels;
colorbar('off')
Set_figure_size_square% clb=colorbar;
%
% clb.Label.String = 'NFkB(S.I.)';


figure(fig_num(2))
h2=heatmap(traj_tnfo,'ColorMap',parula,'GridVisible','off','XLabel',{},'YLabel',{},'ColorLimits',[-0.001,0.25]);

XLabels = 0:1:size(traj_wt,2)-1;
% Convert each number in the array into a string
CustomXLabels = string(XLabels/12);
% Replace all but the fifth elements by spaces
CustomXLabels(:) = " ";
% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h2.XDisplayLabels = CustomXLabels;

YLabels = 1:size(traj_wt,1);
% Convert each number in the array into a string
YCustomXLabels = string(YLabels);
% Replace all but the fifth elements by spaces
YCustomXLabels(:) = " ";
% Set the 'XDisplayLabels' property of the heatmap
% object 'h' to the custom x-axis tick labels
h2.YDisplayLabels = YCustomXLabels;
colorbar('off')
Set_figure_size_square% close

end

%% draw distribution

function [] = draw_distribution(fig_num,colors, parameters_wt,parameters_tnfo,sim_info)
%

figure(fig_num(1))
h1 = histogram(parameters_wt(:,1));hold on
h2 = histogram(parameters_tnfo(:,1));

h1.FaceColor = colors.WT;
h1.Normalization = 'countdensity';
%
h2.FaceColor = colors.TNFKO;
h2.Normalization = 'countdensity';

Set_figure_size_wide_short



if length(fig_num) >1
    figure(fig_num(2))
    h1 = histogram(parameters_wt(:,2));hold on
    h2 = histogram(parameters_tnfo(:,2));
    
    h1.FaceColor = colors.WT;
    h1.Normalization = 'countdensity';
    %
    h2.FaceColor = colors.TNFKO;
    h2.Normalization = 'countdensity';
    Set_figure_size_wide_short
    
end

if length(fig_num) == 3
    figure(fig_num(3))
    scatter(parameters_wt(:,1),parameters_wt(:,2),2,colors.WT,'filled');hold on
    scatter(parameters_tnfo(:,1),parameters_tnfo(:,2),2,colors.TNFKO,'filled')
    
    xlabel(replace(sim_info.parameter_name{1},'params','p'))
    % ylabel('countdensity')
    
    ylabel(replace(sim_info.parameter_name{2},'params','p'))
    
    set(gca,'fontsize',7,'fontweight','b')
    
    Set_figure_size_square
end

end

%%     draw IKK
function [] = draw_IKK(fig_num,colors,IKK_wt,IKK_tnfo)

figure(fig_num)
plot(((1:length(IKK_wt))-1),IKK_wt,'color',colors.WT,'LineWidth',1);hold on
plot(((1:length(IKK_tnfo))-1),IKK_tnfo,'color',colors.TNFKO,'LineWidth',1)

% legend('wt','tnfo')
% ylabel('IKK')
ylim([0,0.05])
yt = yticks;
ytlable = cell(1,length(yt));
for i_yt = 1:length(yt)
    ytlable{i_yt} = '';
end

xlim([0,96*5])
xt = 0:(24*5):(96*5);
xtlable = cell(1,length(xt));
for i_xt = 1:length(xt)
    xtlable{i_xt} = '';
end

%     ylabel('Osc ratio')
set(gca,'YTickLabel',ytlable,'XTick',xt,'XTickLabel',xtlable)

set(gca,'fontsize',7,'fontweight','b')
Set_figure_size_square

end

%% draw TNFR

function [] = draw_TNFR(fig_num,colors,TNFR_wt,TNFR_tnfo)
figure(fig_num)

b=bar([TNFR_wt/TNFR_wt,TNFR_tnfo/TNFR_wt]*100);
%xlabel()
b.FaceColor = 'flat';
b.CData(1,:) = colors.WT;
b.CData(2,:) = colors.TNFKO;

ytickformat(gca, '%g%%');
ax = gca;
ylim([0,150])
yt = yticks;
ytlable = cell(1,length(yt));
for i_yt = 1:length(yt)
    ytlable{i_yt} = '';
end
ax = gca;
ax.YGrid ='on';
%     ylabel('Osc ratio')
set(gca,'YTickLabel',ytlable,'XTickLabel',{'',''})
set(gca,'fontsize',7,'fontweight','b')

Set_figure_size_square

end
