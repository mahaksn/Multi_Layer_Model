clear all; close all; clc;
save_file=1;

folder=['D:\20081\MATLAB Output\Paper-2_Outputs\' ...
    'Paper-2_Multi_layer_domain\FEM_Meshes\'];
time = datestr(datetime('now'),'yyyymmdd_HHMMSS');
prefix = [folder,time];

%%
LS=[0,100,0.1,0.2];
[pS,tS,pS_1DB,tS_1DB]=Rectangle2D(LS,0.03,50,'edge',[1,3],[prefix,'_S']);
% whos tS_1DB tS
bottom_pS = pS_1DB(1,:);
bottom_tS = tS_1DB(1:2,:);
top_pS    = pS_1DB(2,:);
top_tS    = tS_1DB(3:4,:);
% none
if save_file
save('pS_2D.mat','-mat','pS')
save('tS_2D.mat','-mat','tS')
save('bottom_pS_1D.mat','-mat','bottom_pS')
save('bottom_tS_1D.mat','-mat','bottom_tS')
save('top_pS_1D.mat','-mat','top_pS')
save('top_tS_1D.mat','-mat','top_tS')
end

%%
LC=[0,100,0,0.1];
[pC,tC,pC_1DB,tC_1DB]=Rectangle2D(LC,0.03,50,'edge',[1,3],[prefix,'_C']);
bottom_pC = pC_1DB(1,:);
bottom_tC = tC_1DB(1:2,:);
top_pC    = pC_1DB(2,:);
top_tC    = tC_1DB(3:4,:);
if save_file
save('pC_2D.mat','-mat','pC')
save('tC_2D.mat','-mat','tC')
save('bottom_pC_1D.mat','-mat','bottom_pC')
save('bottom_tC_1D.mat','-mat','bottom_tC')
save('top_pC_1D.mat','-mat','top_pC')
save('top_tC_1D.mat','-mat','top_tC')
end

%%
function [p,t,bp,bt]=Rectangle2D(L,dx,d,RegionType,RegionID,prefix)
Lxa=L(1); Lxb=L(2);
Lya=L(3); Lyb=L(4);
model=createpde();
R = [3,4,Lxa,Lxb,Lxb,Lxa,Lya,Lya,Lyb,Lyb]';
g = decsg(R);
geom = geometryFromEdges(model, g);
mesh = generateMesh(model,'Hmax',dx,"GeometricOrder","linear");
p=mesh.Nodes;
t=mesh.Elements; 
for i=1:length(RegionID)
bp(i,:) = findNodes(mesh,"region",RegionType,RegionID(i));
size(bp);
if RegionType=="Edge" || RegionType=="edge"
    if RegionID(i)==3
        bp(i,:)=flip(bp(i,:));
    elseif RegionID(i)==1
        bp(i,:)=bp(i,:);
    else
        error("Boundary edge is not top or bottom")
    end
else
    error("RegionType is not a string or not 'Edge' or 'edge'")
end
bt(2*i-1:2*i,:)=[bp(i,1:end-1);bp(i,2:end)];
end
fig=figure;
subplot(1,3,1)
pdegplot(model,EdgeLabels="on")
axis equal
xlim([Lxa-d,Lxb+d])
ylim([Lya-d,Lyb+d])
subplot(1,3,2)
pdemesh(model,NodeLabels="on",ElementLabels="on")
axis equal
xlim([Lxa-d,Lxb+d])
ylim([Lya-d,Lyb+d])
subplot(1,3,3)
hold on;
for k = 1:size(bt,1)
    idx = bt(k,:);
    plot(p(1,idx), p(2,idx), 'g*-', 'LineWidth', 1);
end
axis equal
xlim([Lxa-d,Lxb+d])
ylim([Lya-d,Lyb+d])
saveas(fig,[prefix,'_mesh.png']);
saveas(fig,[prefix,'_mesh.fig']);

end
