function cats = osc_cats(peakfreq,off_times ,varargin )
% Calculates the percentage of trajectories that fall into different
% oscillation categories
%----------------------------------------------------------------------
% frac = osc_cats(peakfreq,off_times,'cutoff_fq', 0.42 )
%--------------------------------------------------------------------------
% INPUT:
%   REQUIRED:
%       peakfreq:           peak frequencies of nfkb trajectories
%       off_times:          time of signal termination
%   OPTIONAL:
%       'cutoff_fq':          freq cutoff for oscillatory & non-oscillatory
% OUTPUT:
%   tbl:     fraction of cells that are 'on', 'oscillatory', 'transient',
%           'persistent'
p=inputParser; 
addRequired(p, 'peakfreq', @isnumeric); 
addRequired(p, 'off_times', @isnumeric);
addParameter(p, 'cutoff_fq', 0.42, @isnumeric); 
parse(p,peakfreq,off_times ,varargin{:}); 
peakfreq = p.Results.peakfreq; 
cutoff_fq = p.Results.cutoff_fq; 
off_times = p.Results.off_times; 
% cats = cell(size(peakfreq)); 
% cats(off_times==0) = {'off'};
% cats((off_times>0)&(peakfreq<cutoff_fq)) = {'non_osc'};
% cats((peakfreq>=cutoff_fq))={'osc'};

cats = zeros(size(peakfreq)); 
cats(off_times==0) = 1;
cats((off_times>0)&(peakfreq<cutoff_fq)) = 3;
cats((peakfreq>=cutoff_fq))=2;
varNames = {'off', 'osc', 'non_osc'};

cats= categorical(cats,1:3, varNames); 

% valueSet = 1:numel(varNames);
% cats = categorical(grp2idx(cats), valueSet, varNames);


end