
function [proj_path,proj_name] = subfunc_get_proj_path(proj_num,proj_type)
% relative path to Postdoc projects

if nargin<2
    proj_type = 'XGES';
end

proj_name = strcat(proj_type,sprintf( '%04d', proj_num));

proj_path = './NFkB_para_estm_project/SAEM_proj_2022/';

proj_path = strcat(proj_path,proj_name,'/');

end

