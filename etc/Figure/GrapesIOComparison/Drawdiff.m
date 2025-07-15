function [] = Drawdiff(var,lvl)
%DRAWDIFF Summary of this function goes here
%   Detailed explanation goes here
ps_t  = ncread(['../../../output/GroundTruth/grapesinput_',var,'.nc'],var);
ps_rc  = ncread(['../../../output/BeforeRev/grapesinput_',var,'.nc'],var);
ps_rv  = ncread(['../../../output/grapesinput_',var,'.nc'],var);


ps_diff = ps_t-ps_rc;
ps_diff2 = ps_t-ps_rv;
 
t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
sgtitle(var);

nexttile
pcolor(reshape(ps_t(lvl, :,:),67,901).')
shading flat
% axis equal
colorbar
nexttile
pcolor(reshape(ps_rc(lvl, :,:),67,901).')
shading flat
% axis equal
colorbar
nexttile
pcolor(reshape(ps_diff(lvl, :,:),67,901).')
shading flat
% axis equal
colorbar
nexttile
pcolor(reshape(ps_diff2(lvl, :,:),67,901).')
shading flat
% axis equal
colorbar

set(gcf, 'Position',[100,100, 1500, 900])

cc = 1;
end

