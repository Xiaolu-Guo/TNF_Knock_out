
%% main function
function [] = draw_p100o(vers_TNFo,vers_p100oTNFo,data_save_file_path,vers_fig, fig_save_path )
data_info.save_file_path = data_save_file_path;

colors.WT = [93,147,191]/255;
colors.TNFKO = [233,72,73]/255;
colors.TNFKOP100KO = [119 172 48]/255;

data_info.save_file_name = strcat('tnfo_dual_para_sample_',vers_p100oTNFo); % only beginning, no .mat
gene_info.parameter_name_vec = {{'params54','params61'}};

fig_name = 'tnfo_dual_para_sample_';
fig_num_distri =[1,2,3];

genot_vec = {'wt','tnfo'};


for i_para_name = 1:length(gene_info.parameter_name_vec)
    
    %% load sim data
    sim_info.parameter_name = gene_info.parameter_name_vec{i_para_name};
    
    vers = vers_p100oTNFo;
    data_info.save_file_name = strcat('tnfo_dual_para_sample_',vers); % only beginning, no .mat
    
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
    metricsp100 = metrics;
    datap100 = data;
    index_nuc_NFkB_tnfo = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='nucNFkB'...
        & sim_data_tbl.type == 'tnfo';   
    traj_wt = sim_data_tbl.trajectory(index_nuc_NFkB_tnfo,:);
    [~,traj_order] = sort(max(traj_wt,[],2),'descend');
    traj_wt = traj_wt(traj_order,:);
    traj_tnfop100o = traj_wt(:,1:5:end);
    
    index_IkBat_p100otnfo = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='IkBat'...
        & sim_data_tbl.type == 'tnfo';
    IkBat_p100otnfo = mean(sim_data_tbl.trajectory(index_IkBat_p100otnfo,:));
    
    vers = vers_TNFo;
    data_info.save_file_name = strcat('tnfo_dual_para_sample_',vers); % only beginning, no .mat
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
    index_nuc_NFkB_tnfo = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='nucNFkB'...
        & sim_data_tbl.type == 'tnfo';   
    traj_tnfo = sim_data_tbl.trajectory(index_nuc_NFkB_tnfo,:);
    [~,traj_order] = sort(max(traj_tnfo,[],2),'descend');
    traj_tnfo = traj_tnfo(traj_order,:);
    traj_tnfo = traj_tnfo(:,1:5:end);
    
    index_IkBat_wt = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='IkBat'...
        & sim_data_tbl.type == 'wt';
    index_IkBat_tnfo = sim_data_tbl.dose_val == 10 & sim_data_tbl.species =='IkBat'...
        & sim_data_tbl.type == 'tnfo';
    IkBat_wt = mean(sim_data_tbl.trajectory(index_IkBat_wt,:));
    IkBat_tnfo = mean(sim_data_tbl.trajectory(index_IkBat_tnfo,:));
        
    %% draw the picture and save
    
    if 0 % osc ratio bar
        draw_osc_ratio_bar_p100o([1,2,3],colors,metrics,data,metricsp100,datap100)
        
        figure(1)
        %     saveas(gcf,strcat(fig_save_path,fig_save_name,'_heatmap_byospower'),'epsc')
        %     saveas(gcf,strcat(fig_save_path,fig_save_name,'_heatmap_byospower'),'svg')
        close
        
        figure(2)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_oscpercentage_p100o_diff_thresh'),'fig')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_oscpercentage_p100o_diff_thresh'),'svg')
        close
        
        figure(3)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_oscpercentage_p100o'),'fig')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_oscpercentage_p100o'),'svg')
        close
        
    end
    
    if 1 % draw heatmap
        draw_traj_heatmap([1,2],traj_tnfop100o,traj_tnfo)
        
        figure(1)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_fig4_trajheatmap_tnfop100o'),'fig')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_fig4_trajheatmap_tnfop100o'),'svg')
        close
        
        figure(2)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_fig4_trajheatmap_tnfo'),'fig')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_fig4_trajheatmap_tnfo'),'svg')
        close
    end
    
    if 0 % draw IkBat
        fig_num = 1 ;        
        draw_IkBat(fig_num,colors,IkBat_wt,IkBat_tnfo)
        figure(fig_num)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_IkBat_wt_tnfo'),'fig')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_IkBat_wt_tnfo'),'svg')
        close
        
        colors2.WT = colors.TNFKO;
        colors2.TNFKO = colors.TNFKOP100KO;
        draw_IkBat(fig_num,colors2,IkBat_tnfo,IkBat_p100otnfo)
        figure(fig_num)
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_IkBat_tnfo_p100otnfo'),'fig')
        saveas(gcf,strcat(fig_save_path,fig_save_name,'_IkBat_tnfo_p100otnfo'),'svg')
        close
    end
    
    if 0 %save oscpower for TNFKO and TNFKOP100KO
        colors_input = {colors.TNFKO, colors.TNFKOP100KO};
        i_genotype = 2;

        TNFo_osc = metrics{i_genotype}.oscpower;
        TNFop100o_osc = metricsp100{i_genotype}.oscpower;
        save('../raw_data/tnfop100o_oscpower.mat','TNFo_osc','TNFop100o_osc');
        
        vects= {TNFo_osc TNFop100o_osc};
        ax = gca;
        %ax=gca;
        id_begin =1;
        id_end = 2;
        ymax = max(cell2mat(vects(:)));
        ymin = min(cell2mat(vects(:)));
        ylimits = [0  1.1*ymax];
        data_type =1:2;
        
        % draw_violin_plot_p100o(colors_input,data_type,vects(id_begin:id_end), [0,0.25], 'Axes', ax, 'YLim', ylimits)
        
    end
end

end

%% draw oscllatory cells

function [] = draw_osc_ratio_bar_p100o(fig_num,colors,metrics,data,metricsp100,datap100)

threshold = [0.5e-4, 1e-4, 1.5e-4, 2e-4, 2.5e-4];%[0.5e-4,1e-4,1.5e-4];(For Xiaolu's oscpower)
osc_ratio = NaN(2,length(threshold));

i_genotype = 2;
for i_threshold =1:length(threshold)
    osc_ratio(1,i_threshold) = sum(metrics{i_genotype}.oscpower > threshold(i_threshold) )/size(data.model_sim{i_genotype},1);
end
for i_threshold =1:length(threshold)
    osc_ratio(2,i_threshold) = sum(metricsp100{i_genotype}.oscpower > threshold(i_threshold) )/size(datap100.model_sim{i_genotype},1);
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
    b.CData(1,:) = colors.TNFKO;
    b.CData(2,:) = colors.TNFKOP100KO;
    
    ylim([0,100])
    
    ytickformat(gca, '%g%%');
    yt = yticks;
    ytlable = cell(1,length(yt));
    for i_yt = 1:length(yt)
        ytlable{i_yt} = '';
    end
    ax = gca;
    % ax.YGrid = 'on';
    %     ylabel('Osc ratio')
    set(gca,'YTickLabel',ytlable,'XTickLabel',{'',''})
    %     set(gca,'fontsize',7,'fontweight','b')
    
    
end
Set_figure_size_square


figure(fig_num(3))

i_osratio = ceil(4);
b=bar(osc_ratio(:,i_osratio)*100);
%xlabel()
b.FaceColor = 'flat';
b.CData(1,:) = colors.TNFKO;
b.CData(2,:) = colors.TNFKOP100KO;

ylim([0,100])

ytickformat(gca, '%g%%');
yt = yticks;
ytlable = cell(1,length(yt));
for i_yt = 1:length(yt)
    ytlable{i_yt} = '';
end
ax = gca;
% ax.YGrid = 'on';
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
color_mp = colormap('parula');
%color_mp = color_mp(1:end-111,:);

h1=heatmap(traj_wt,'ColorMap',color_mp,'GridVisible','off','XLabel',{},'YLabel',{},'ColorLimits',[-0.001,0.4]);

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
h2=heatmap(traj_tnfo,'ColorMap',color_mp,'GridVisible','off','XLabel',{},'YLabel',{},'ColorLimits',[-0.001,0.4]);

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

%% draw IkBat

function [] = draw_IkBat(fig_num,colors,IKK_wt,IKK_tnfo)

figure(fig_num)
plot(((1:length(IKK_wt))-1),IKK_wt,'color',colors.WT,'LineWidth',1);hold on
plot(((1:length(IKK_tnfo))-1),IKK_tnfo,'color',colors.TNFKO,'LineWidth',1)

% legend('wt','tnfo')
% ylabel('IKK')
ylim([0,2e-4])
yt = yticks;
ytlable = cell(1,length(yt));
for i_yt = 1:length(yt)
    ytlable{i_yt} = '';
end

xlim([0,24*5])
xt = 0:(6*5):(24*5);
xtlable = cell(1,length(xt));
for i_xt = 1:length(xt)
    xtlable{i_xt} = '';
end

%     ylabel('Osc ratio')
set(gca,'YTickLabel',ytlable,'XTick',xt,'XTickLabel',xtlable)

set(gca,'fontsize',7,'fontweight','b')
Set_figure_size_square

end

%% draw IKK
function [] = draw_IKK(fig_num,colors,IKK_wt,IKK_tnfo)

figure(fig_num)
plot(((1:length(IKK_wt))-1),IKK_wt,'color',colors.WT,'LineWidth',1);hold on
plot(((1:length(IKK_tnfo))-1),IKK_tnfo,'color',colors.TNFKO,'LineWidth',1)

% legend('wt','tnfo')
% ylabel('IKK')
ylim([0,2e-4])
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

%% draw violin plot of oscpower

function [] = draw_violin_plot_p100o(colors_input,data_type,vects,places,varargin)

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [] = violin(vects, places, varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% VIOLIN creates a violin plot, spacing them according to a secondary vector (e.g. doses)
%
% INPUTS (required)
% vects          1xN cell array of vectors (1-D array of object measurements)
% places         1xN array directing placement of each Violin
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Color'        1x... cell vector specifying violin fill colors - cycles if length < N
% 'YLim'         2 element vector with graph [y_min, y_max]. Default is 5th and 95th percentile of all data
% 'ShowBins'     Show additional histogram  graph ('on' or 'off' - default is 'off')
% 'Area'         Total area of graph taken up by each shape (default = 0.01)
% 'BinScale'     Scaling factor (from default) for number of histogram bins (can specify scalar or a vector of length N)
% 'Bins'         Vector of bin centers - if provided, 'BinScale' will be ignored
% 'XSpace'       Axis whitespace before the first violin plot, and after the last (default = 0.1)
% 'Axes'         Axes handle of axes where new violin figure will be plotted (default: create new figure)
% 'Smoothing'    Smoothing of violin shapes (such that histogram bars are evident, or smoothed out). Default = 'on'
% 'Connect'      Add connecting line between shapes (default = 'on')
% 'LineWidth'    Line width around shape - default = 0.5 (can be 0)
% 'MarkerSize'   Size of median marker (if shown). Default = 7.
%
% OUTPUTS
% violin        Axes handle of violin figure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% 1) Vector data (must be cell matric)
addRequired(p,'vects',@iscell);
% 2) X axis Placement
valid_places = @(x) assert(numel(x)==numel(vects), '2nd argument must be same size as 1st');
addRequired(p,'places',valid_places);

% Optional parameters
colors = setcolors;
default_color = colors.peacock(end:-1:1);

valid_color = @(x) assert(iscell(x)&&length(x{1})==3, 'Specify colors with a cell matrix of RGB triplets');
addParameter(p,'Color', default_color,valid_color);

all = cell2mat(vects(:));
valid_ylim = @(x) assert(length(x)==2,'YLim must be a 2 element vector');
addParameter(p,'YLim', prctile(all(:),[5 95]),valid_ylim);
expectedFlags = {'on','off'};
addParameter(p,'ShowBins','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'Area',0.01,@isnumeric);
addParameter(p,'XSpace',0.1,@isnumeric);
addParameter(p,'BinScale',1,@isnumeric);
valid_width = @(x) assert(isnumeric(x)&&(x>=0),'Line width/marker size must be >= 0');
addParameter(p,'LineWidth',0.5,valid_width);
valid_bins = @(x) assert((length(x)>1) && isnumeric(x) && issorted(x), 'Bins must be monotonically increasing vector');
addParameter(p,'Bins',nan,valid_bins);
addParameter(p,'Axes',nan,@ishandle);
addParameter(p,'Smoothing','on', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'Connect','on', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'MarkerSize',7,valid_width);


% Parse inputs, save some to variables
parse(p,vects,places, varargin{:})
xspace = p.Results.XSpace;
bin_scale = p.Results.BinScale;
bins = p.Results.Bins;
ylim = p.Results.YLim;
% colors = p.Results.Color;
colors = colors_input;

% Create figure (if axes wasn't provided)
if ~ishandle(p.Results.Axes)
    viofig = figure('Position', [500, 1031, 800, 300], 'PaperPositionMode','auto');
    violin = axes('Parent',viofig);
else
    violin = p.Results.Axes;
end

% Get medians of all sets
medians = cellfun(@nanmedian,vects);

% Create figures; set XLim
if strcmp(p.Results.ShowBins,'on') % Diagnostic output: histogram overlaid with spline fit
    figure('Position',[500,357, 350 100*length(vects)]);
    ha = tight_subplot(length(vects),1);
end
if length(places)>1
    x_lim = [min(places)-(1.5*xspace*range(places)),max(places)+1.5*xspace*range(places)];
else
    x_lim = sort([places*0.9,places*1.1],'ascend');
end
tot_area = abs(diff(ylim(1:2)))*abs(diff(x_lim));

% Make bins
bin_scale = repmat(bin_scale,1,length(vects));
bin_scale = bin_scale(1:length(vects));

% Loop through sets, generate shapes
for i = 1:length(vects)
    if isnan(bins)
        bin_width = 2*iqr(all)*((numel(all)/length(vects))^(-1/3))/bin_scale(i);
        x = prctile(all,0):bin_width:prctile(all,100); %changed from 99.5 percentile to 100 for upper bound; changed from 1 to 0 for lower bound
    else
        bin_width = bins(2)-bins(1);
        x = bins(:)';
    end
    
    
    % Generate histogram data
    y = hist(vects{i},x);
    
    % Cap histogram with zero values (keep spline from spiking @ end) and interpolate to get shape
    if y(1)==0
        pos = find(y>0,1,'first');
        y(1:(pos-1)) = [];
        x(1:(pos-1)) = [];
    end
    if y(end)==0
        pos = find(y>0,1,'last');
        y((pos+1):end) = [];
        x((pos+1):end) = [];
    end
    
    x = [min(x)-bin_width, x, max(x)+bin_width];
    y = [0 y/sum(y) 0];
    
    if strcmpi(p.Results.Smoothing,'on')
        xx = min(x):bin_width/10:max(x);
        yy = spline(x, y, xx);
        yy(yy<0) = 0;
    else
        xx = sort([x-bin_width/2.001, x+bin_width/2.001]);
        yy = repmat(y,2,1);
        yy = yy(:)';
    end
    
    % (optionally) show subplot of bins+spline fit
    if strcmp(p.Results.ShowBins,'on')
        hold(ha(i),'on')
        bar(ha(i),x,y,'FaceColor',[45 191 104]/255,'EdgeColor','none')
        set(ha(i),'XLim',ylim,'YLim',[0 .5])
        plot(ha(i),xx,yy,'LineWidth',2,'Color',[0 0 0])
        hold(ha(i),'off')
    end
    
    % Scale shape width so total area is consistent
    obj_width = p.Results.Area*tot_area/(sum(y)*diff(x(1:2))*2);
    
    % Make main violin plot
    if p.Results.LineWidth==0
        lnstyl = 'none';
        lnwid = 0.5;
    else
        lnstyl = '-';
        lnwid = p.Results.LineWidth;
    end
    
    hold(violin,'on')
    
    fill([places(i)+obj_width*yy,places(i)-obj_width*yy(end:-1:1)],[xx,xx(end:-1:1)],...
        colors{data_type(mod(i-1,length(data_type))+1)},'LineWidth',1,'Parent',violin,'LineStyle',lnstyl,'LineWidth',lnwid,...
        'EdgeColor',[0.4431 0.4510 0.4627])
    
    hold(violin,'off')
end
% Plot medians and set graph properties
hold(violin,'on')
if strcmpi(p.Results.Connect,'on')
    lnstyl = '-o';
else
    lnstyl = 'o';
end

% plot(violin, places,medians,lnstyl,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],...
%     'Color', [0 0 0],'LineWidth',1.2,'MarkerSize',p.Results.MarkerSize)
scatter(violin, places,medians,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0,0,0],'SizeData',p.Results.MarkerSize*3);
hold(violin,'off')
set(violin,'YLim',ylim,'XLim',x_lim);

end
