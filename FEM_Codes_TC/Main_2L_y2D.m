clear all; close all; clc;
tic

%% pre-setting
showanimation=1;
makegif=1;
drawperframe=500;
pattern_detected = false;
tol=1e-9;
% dx=0.115;
% dy=0.115;
T=5000;
dt=0.001;

%% Load Mesh Data: p-nodes; t-elements;
load pC_2D.mat
load tC_2D.mat
load pS_2D.mat
load tS_2D.mat
load top_pC_1D.mat
load top_tC_1D.mat
load bottom_pC_1D.mat
load bottom_tC_1D.mat
load top_pS_1D.mat
load top_tS_1D.mat
load bottom_pS_1D.mat
load bottom_tS_1D.mat
nC=length(pC);
nS=length(pS);
% nC = size(pC,2);
% nS = size(pS,2);

LS=[0,100,0.25,0.5];
LC=[0,100,0,0.25];

%% Parameters
% Hs=abs(LS(3)-LS(4)); Hc=abs(LC(3)-LC(4)); k=2;
% D_S=1; D_C=1;
% etaS=10;%[6 10 20]; 
% etaC=etaS;

HS=abs(LS(3)-LS(4)); HC=abs(LC(3)-LC(4));
% a = 2.1; b=1.075;
% a=1.82; b=1;
% a=2.1; b=1.1;
a=0.05; b=0.5;
etaS = 0.25;%[0.1 1 10 0.01 0.001]; 
etaC=etaS;
D_S=1; D_C=20;

%% Steady states
% uS0=a; uC0=b/a;
uS0=a+b; uC0=b/(a+b)^2;

%% Define the Source Term and Coupling Term
% f_S=@(u) 0*u;
% f_C=@(v) 0*v;

f_S=@(u) a - u;
f_C=@(v) b - 0*v;
% 
G_S=@(us,uc) us.^2.*uc;
G_C=@(uc,us) us.^2.*uc;
% G_C=@(uc,us) us.*uc.^2;

% G_S=@(us,uc) us.*uc;
% G_C=@(uc,us) us.*uc;

%% saving folder and name
folder=['D:\20081\MATLAB Output\Paper-2_Outputs\' ...
    'Paper-2_Multi_layer_domain\FEM_2_y2D\'];
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
% fprintf(['\nInitial steady state for\n' ...
%     'surface layer is: uS0=%.4f\n'...
%     'core layer is: uC0=%.4f\n'], ...
%     uS0,uC0);

%% time discretization
% T=1000;
% dt=0.01; %1e-4;
nt=T/dt+1;
stopti=nt; % stop patterning at nt
nFrame=ceil((T/dt)/drawperframe);

%% Assemble Matrices
ord=3;tic
[S_S,M_S]=AssembleGlobalMatrices_corrected2(pS,tS,ord);
[S_C,M_C]=AssembleGlobalMatrices_corrected2(pC,tC,ord);
M_big = blkdiag(M_S,M_C);
S_big = blkdiag(D_S*S_S, D_C*S_C);
toc

%%
xxS=linspace(LS(1),LS(2),length(bottom_pS));
yyS=linspace(LS(3),LS(4),length(bottom_pS));
[xqS,yqS]=meshgrid(xxS,yyS);
xxC=linspace(LC(1),LC(2),length(top_pC));
yyC=linspace(LC(3),LC(4),length(top_pC));
[xqC,yqC]=meshgrid(xxC,yyC);

%% Preallocate storage
max_steps = nt;
tt = zeros(max_steps,1);  % store time values

% csh = cosh(HS/sqrt(D_S));
% ssh = sinh(HS/sqrt(D_S));
% A_sol_prev=1e-6;

%%
for e=1:length(etaS)
evalS=round(etaS(e),4);
evalC=round(etaC(e),4);
fprintf('\nFor etaS=%.4f, etaC=%.4f,\n',evalS, evalC);
eval=['_',num2str(evalS),'_',num2str(evalC)];
pattern_detected = false;

pattern_end=5;
stopti=nt; 

%% Initial condition
rng(0); % reproducible random seed
Perturbations2 = 1+0.1*randn(nS,1);
Perturbations1 = 1+0.1*randn(nC,1);
uS=uS0 * Perturbations2;
uC=uC0 * Perturbations1;

% Aa = b*HC/(sqrt(D_S)*sinh(HS/sqrt(D_S)));
% uS0 = Aa*cosh((pS(2,:)' - (HC+HS))/sqrt(D_S)) + a;
% % C2 = -b*HC/(etaS(e)*( (-b*HC/sqrt(D_S))*coth(HS/sqrt(D_S)) + a )^2) + b/(2*D_C)*HC^2;
% % uC0 = - (b/(2*D_C))*pC(2,:)'.^2 + C2;
% Sval = a + (b*HC/sqrt(D_S))*coth(HS/sqrt(D_S));
% C2 = b*HC/( etaS(e) * Sval^2 ) + b/(2*D_C)*HC^2;
% uC0 = - (b/(2*D_C))*pC(2,:)'.^2 + C2;
% 
% eps = 1e-2;
% pert_interior_C = eps * randn(nC,1);
% pert_interior_S = eps * randn(nS,1);
% 
% uC = uC0 + pert_interior_C;
% uS = uS0 + pert_interior_S;

%% set up figure
giffile = [prefix,eval,'_pattern','.gif'];
fig = figure('Color','w');%,'WindowState', 'maximized');

vqS=griddata(pS(1,:),pS(2,:),uS,xqS,yqS);
subplot(2,1,1)
hold on
uS_fig=surf(xqS,yqS,vqS,'linestyle','none');
figtitle=title('t=0');
xlim([LS(1),LS(2)])
ylim([LS(3),LS(4)])
xlabel('$x$','Interpreter','latex');
ylabel('$u_s,~y$','Interpreter','latex');
zlabel('$u_s$','Interpreter','latex');
ax = gca; 
ax.FontSize = 14;
view(2)
colorbar
hold off

vqC=griddata(pC(1,:),pC(2,:),uC,xqC,yqC);
subplot(2,1,2)
hold on
uC_fig=surf(xqC,yqC,vqC,'LineStyle','none');
xlim([LC(1),LC(2)])
ylim([LC(3),LC(4)])
xlabel('$x$','Interpreter','latex');
ylabel('$u_c,~y$','Interpreter','latex');
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
F_S = ReactKineInt(pS,tS,uS,f_S,ord);
F_C = ReactKineInt(pC,tC,uC,f_C,ord);
% F_S=0*uS;
% F_C=0*uC;

ord=2;
bS_bottom = NLBoundFluxInt(bottom_tS, top_tC, pS, pC, uS, uC, G_S, etaS(e), ord);
bC_top    = NLBoundFluxInt(top_tC, bottom_tS, pC, pS, uC, uS, G_C, -etaC(e), ord);

if any(isinf(bS_bottom(:))) || any(isinf(bC_top(:)))
    error('Inf at the boundary conditions  at step %d (time = %.5f)', ti, t);
end

%% Matrix system
AuS = M_S + dt/2*D_S*S_S;
BuS = (M_S - dt/2*D_S*S_S)*uS + dt*(F_S + bS_bottom);
% % --- Surface top: uS = k ---
% AuS(top_pS,:) = 0;                               % zero the row
% AuS(sub2ind(size(AuS), top_pS, top_pS)) = 1;     % put 1 on the diagonal
% BuS(top_pS) = k;                                 % enforce value

uS_new = AuS\BuS;

AuC = M_C + dt/2*D_C*S_C;
BuC = (M_C - dt/2*D_C*S_C)*uC + dt*(F_C + bC_top);
% % --- Core bottom: uC = 0 ---
% AuC(bottom_pC,:) = 0;
% AuC(sub2ind(size(AuC), bottom_pC, bottom_pC)) = 1;
% BuC(bottom_pC) = 0;                              % enforce value

uC_new = AuC\BuC;

%%
U_big = [uS; uC];
F_big = [F_S; F_C];
B_big = [bS_bottom; bC_top];
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

    %% Plot and gif
    if mod(ti, drawperframe) == 1
        vqS = griddata(pS(1,:),pS(2,:),uS,xqS,yqS);
if any(isnan(vqS(:)))
    warning('griddata returned NaN values for surface interpolation.');
end
        vqC = griddata(pC(1,:),pC(2,:),uC,xqC,yqC);
        if showanimation
            uS_fig.ZData = vqS;
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

%%
fprintf('\nDone!\n');
toc
diary OFF