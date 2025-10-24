clear all; close all; clc;
tic

%% pre-setting
showanimation=1;
makegif=1;
drawperframe=10;
pattern_detected = false;
tol=1e-7;
% dx=0.115;
% dy=0.115;
T=2;

%% Load Mesh Data: p-nodes; t-elements;
LS=[0,0.1];
LC=[0,0.1,0,0.1];
% [pC,tC,top_pC,top_tC]=Rectangle2D(LC,0.03,0,'edge',3);
load pC_2D.mat
load tC_2D.mat
load top_pC_1D.mat
load top_tC_1D.mat
nC=length(pC);

Ne =  length(top_pC);
pS = linspace(LS(1), LS(2), Ne);
tS = [1:Ne-1; 2:Ne]'; 
nS=length(pS);

%% Parameters
a = 2; b=3; etaS = 0.1; etaC=0.1; ep=0.01;
D_S=1; D_C=9;
%% Steady states
uS0=a; uC0=b/a;
%% Define the Source Term and Coupling Term
f_S=@(u,v,eta) a - u - ep*u.^3 + eta*(u.^2.*v - b*u);
f_C=@(v) ep*(1-v);

G_C = @(v, u) b*u - u.^2.*v;

%% Assemble Matrices
ord=3;tic
[S_S,M_S]=AssembleGlobalMatrices1D(pS,tS,ord);
[S_C,M_C]=AssembleGlobalMatrices(pC,tC,ord);
M_big = blkdiag(M_S,M_C);
S_big = blkdiag(D_S*S_S, D_C*S_C);
toc

%% saving folder and name
folder=['D:\20081\MATLAB Output\' ...
    'Paper-2_Multi_layer_domain\FEM_2_1D2D\'];
time = datestr(datetime('now'),'yyyymmdd_HHMMSS');
prefix = [folder,time];

%% saving code
options = struct('format','pdf', ...
    'outputDir',folder, ...
    'codeToEvaluate','0;');
% publish(pdfname,options);
diary([prefix,'.txt']);
fprintf('saving to %s\n',folder);

%%
FileNameAndLocation = mfilename('fullpath');   % current script full path without extension
[filepath,name,~] = fileparts(FileNameAndLocation);
ext='.m';
% Original script full filename
origFile = fullfile(filepath,[name ext]);

% Backup filename (adds "backup" + version + .txt)
newbackup = fullfile(folder, sprintf('%s_%s.txt',time,name));

% Check if backup already exists
A = exist(newbackup,'file');
if (A ~= 0)
    warning('Backup already exists for the current version')
else
    % Create backup by copying the current .m file
    copyfile(origFile, newbackup);
    fprintf('Backup created: %s\n', newbackup);
end

%% Steady states
fprintf(['\nInitial steady state for\n' ...
    'surface layer is: uS0=%.4f\n'...
    'core layer is: uC0=%.4f\n'], ...
    uS0,uC0);

%% time discretization
% T=1000;
dt=0.01; %1e-4;
nt=T/dt+1;
stopti=nt; % stop patterning at nt
nFrame=ceil((T/dt)/drawperframe);

%%
xxC=linspace(LC(1),LC(2),length(top_pC));
yyC=linspace(LC(3),LC(4),length(top_pC));
[xqC,yqC]=meshgrid(xxC,yyC);

%% Preallocate storage
max_steps = nt;
uS_store = zeros(max_steps,nS);
tt = zeros(max_steps,1);

%%
for e=1:length(etaS)
evalS=round(etaS(e),4);
evalC=round(etaC(e),4);
fprintf('\nFor etaS=%.4f, etaC=%.4f,\n',evalS, evalC);
eval=['_',num2str(evalS),'_',num2str(evalC)];

pattern_end=10;
stopti=nt; % stop patterning at nt

%% Initial condition
pit=0; rng(pit); %change random seed
PerturbationsS = 1+0.1*randn(nS,1);
PerturbationsC = 1+0.1*randn(nC,1);
uS=uS0 * PerturbationsS;
uC=uC0 * PerturbationsC;

%% set up figure
giffile = [prefix,eval,'_pattern','.gif'];
fig = figure('Color','w');%,'WindowState', 'maximized');

subplot(2,1,1)
hold on
uS_fig=plot(pS, uS, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
xlim([LS(1),LS(2)]);
% ylim([]);
xlabel('$x$','Interpreter','latex');
ylabel('$u_s$','Interpreter','latex');
figtitle=title('t=0');
ax = gca;
ax.FontSize = 14;
hold off

vqC=griddata(pC(1,:),pC(2,:),uC,xqC,yqC);
subplot(2,1,2)
hold on
uC_fig=surf(xqC,yqC,vqC,'LineStyle','none');
xlim([LC(1),LC(2)])
ylim([LC(3),LC(4)])
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$u_c$','Interpreter','latex');
ax = gca; 
ax.FontSize = 14;
view(2)
colorbar
hold off  

%% Simulation iteration
for ti = 1:nt
    t = dt * (ti - 1);
    tt(ti) = t;

%% --- Right-hand side assembly ---
ord=3;
F_S = ReactKineInt1D(pS,tS,uS,uC,etaS(e),f_S,ord);
F_C = ReactKineInt(pC,tC,uC,f_C,ord);

ord=2;
bC_top = NLBoundFluxInt(top_tC, tS, pC, pS, uC, uS, G_C, etaC(e), ord);

%% Matrix system
AuS = M_S + dt/2*D_S*S_S;
BuS = (M_S - dt/2*D_S*S_S)*uS + dt*(F_S);
uS_new = AuS\BuS;
AuC = M_C + dt/2*D_C*S_C;
BuC = (M_C - dt/2*D_C*S_C)*uC + dt*(F_C + bC_top);
uC_new = AuC\BuC;

%%
U_big = [uS; uC];
F_big = [F_S; F_C];
B_big = [F_S*0; bC_top];
A = M_big + 0.5*dt*S_big;
B = (M_big - 0.5*dt*S_big)*U_big + dt*(F_big + B_big);

%% Stability and NaN check
if any(isnan(A(:))) || any(isnan(B(:)))
 [any(isnan(A(:)))   %1
any(isnan(B(:)))    %2
any(isnan(uS(:)))   %3
any(isnan(uC(:)))   %4
any(isnan(B_big(:)))    %5
any(isnan(F_big(:)))    %6
any(isnan(bS_bottom(:)))    %7
any(isnan(bC_top(:)))]    %8
   error('Matrix A or RHS B contains NaN at step %d (time = %.5f)', ti, t);
end
if condest(A) > 1e12
    warning('Matrix A is ill-conditioned at step %d (time = %.5f), condest = %.2e', ti, t, condest(A));
end

    %% Pattern detection using norms
    errS = norm(uS_new - uS)/norm(uS);
    errC = norm(uC_new - uC)/norm(uC);
    maxErr = max([errS, errC]);
   
%% Update and save
uS = uS_new;
uC = uC_new;

uS_store(ti,:) = uS;

%% Plot and gif
    if mod(ti, drawperframe) == 1
        % vqS = griddata(pS(1,:),pS(2,:),uS,xqS,yqS);
        vqC = griddata(pC(1,:),pC(2,:),uC,xqC,yqC);
        if showanimation
            uS_fig.YData = uS;
            uC_fig.ZData = vqC;
            figtitle.String = ['t = ', num2str(t, '%.4f')];
            drawnow;
        end
        if makegif
            frame = getframe(fig);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if ti == 1
                imwrite(imind, cm, giffile, 'gif', 'Loopcount', inf);
            else
                imwrite(imind, cm, giffile, 'gif', 'WriteMode', 'append', 'DelayTime', 0);
            end
        end
    end
    if ~pattern_detected && maxErr < tol
        fprintf('\nPattern fully developed at t = %.5f\n', t);
        pattern_detected = true;
        stopti = ti + round(pattern_end / dt);
    end

    if ti == stopti
        break;
    end
end
%% saving final pattern
saveas(fig,[prefix,eval,'_final.png']);
saveas(fig,[prefix,eval,'_final.fig']);

end

%% --- Space-Time Visualization ---
figure;
imagesc(1:nS, tt(1:ti), uS_store(1:ti,:));
ylabel('Time'); xlabel('Node Index'); title('u_S Space-Time'); colorbar;

%%
fprintf('\nDone!\n');
toc
diary OFF