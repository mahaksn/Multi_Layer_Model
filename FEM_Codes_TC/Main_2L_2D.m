clear all; close all; clc;
tic

%% pre-setting
showanimation=1;
makegif=1;
drawperframe=500;
pattern_detected = false;
tol=1e-8;
% dx=0.115;
% dy=0.115;
T=5000;
dt=0.001;
pdfname="Main_2L_2D";

%% Load Mesh Data: p-nodes; t-elements;
load pC_2D.mat
load tC_2D.mat
load pS_2D.mat
load tS_2D.mat
load top_pC_1D.mat
load top_tC_1D.mat
load bottom_pS_1D.mat
load bottom_tS_1D.mat
nC=length(pC);
nS=length(pS);

LS=[0,100,0.1,0.2];
LC=[0,100,0,0.1];

%% Parameters
% a = 0.01; b=3; etaS = 0.1; etaC=1; ep=0.001;
% D_S=1; D_C=100;
% a = 2; b=3; etaS = 0.5; etaC=0.5; ep=0.01;
% D_S=1; D_C=9;
% a = 0.05; b=0.5; etaS = 0.05; etaC=etaS;
% D_S=1; D_C=20;

% a=-1; b=-3/2; H=2; D=0.08; D_S=D; D_C=1;
% C=1; etaS=1; etaC=etaS;
% Del=5*C^2-8*H+12;

HS=abs(LS(3)-LS(4)); HC=abs(LC(3)-LC(4));
% a = 2.1; b=1.075; 
% a=0.1; b=0.115;
a=2.1; b=1.1;
etaS = [0.1 1];
etaC=etaS;
D_S=1; D_C=20;

%% Steady states
% uS0=1+a+b; uC0=3+b/(a+b)^2;
% uS0=a; uC0=b;
% uS0=a; uC0=b/a;
% uS0=a+b; uC0=b/(a+b)^2;
% uS0=(-5*C+sqrt(5*Del))/(4*H+4);
% uC0=(-5*C+sqrt(5*Del))/(10);

%% Define the Source Term and Coupling Term
% f_S=@(u) a - u - ep*u.^3;
% f_C=@(v) ep*(1-v);%b - v;

% G_S = @(u, v) u.^2.*v - b*u;  % Nonlinear function G_S(uS, uB)
% G_C = @(v, u) b*u - u.^2.*v;  % Nonlinear function G_C(uC, uB)

% G_S = @(u,v) u-v;
% G_C = @(u,v) u-v;

% f_S=@(u) a - u;
% f_C=@(v) b + 0*v;
% 
% G_S = @(u, v) u.^2.*v;
% G_C = @(v, u) u.^2.*v;

% f_S=@(u,v,eta) a - u;
% f_C=@(u,v,eta) b + 0*v;
% 
% G_S = @(u, v) u.^2.*v;
% G_C = @(u, v) u.^2.*v;

% f_S=@(u) u;
% f_C=@(v) b*v;
% 
% G_S = @(u, v) a*v-C*u.*v-u.*v.^2;
% G_C = @(v, u) H*u+C*u.*v+u.*v.^2;

f_S=@(u) a - u;
f_C=@(v) b + 0*v;

G_S=@(us,uc) us.^2.*uc;
G_C=@(uc,us) us.*uc.^2;

fc_S=@(us,uc,eta) a - us + eta*(us.^2.*uc);
fc_C=@(us,uc,eta) b - eta*(us.*uc.^2);

%% saving folder and name
folder=['D:\20081\MATLAB Output\Paper-2_Outputs\' ...
    'Paper-2_Multi_layer_domain\FEM_2_2D\'];
time = datestr(datetime('now'),'yyyymmdd_HHMMSS');
prefix = [folder,time];

%% saving code
% options = struct('format','pdf','outputDir',folder,'evalCode',false);
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
[S_S,M_S]=AssembleGlobalMatrices(pS,tS,ord);
[S_C,M_C]=AssembleGlobalMatrices(pC,tC,ord);
% M_big = blkdiag(M_S,M_C);
% S_big = blkdiag(D_S*S_S, D_C*S_C);
toc

%%
% xxS=LS(1):dx:LS(2);
% yyS=LS(3):dy:LS(4);
xxS=linspace(LS(1),LS(2),length(bottom_pS));
yyS=linspace(LS(3),LS(4),length(bottom_pS));
[xqS,yqS]=meshgrid(xxS,yyS);
% xxC=LC(1):dx:LC(2);
% yyC=LC(3):dy:LC(4);
xxC=linspace(LC(1),LC(2),length(top_pC));
yyC=linspace(LC(3),LC(4),length(top_pC));
[xqC,yqC]=meshgrid(xxC,yyC);

%% Preallocate storage
max_steps = nt;
% uS_store = zeros(N3, max_steps);
% uB_store = zeros(N2, max_steps);
% uC_store = zeros(N1, max_steps);
tt = zeros(max_steps,1);  % store time values

csh = cosh(HS/sqrt(D_S));
ssh = sinh(HS/sqrt(D_S));
A_sol_prev=1e-6;

%%
for e=1:length(etaS)
evalS=round(etaS(e),4);
evalC=round(etaC(e),4);
fprintf('\nFor etaS=%.4f, etaC=%.4f,\n',evalS, evalC);
eval=['_',num2str(evalS),'_',num2str(evalC)];
pattern_detected = false;

pattern_end=10;
stopti=nt; % stop patterning at nt

%% Initial condition
pit=0; rng(pit); %change random seed
% Perturbations2 = 0.01*(2*rand(N2,1)-1);
% Perturbations1 = 0.01*(2*rand(N1,1)-1);
% uS=uS0 + Perturbations2;
% uC=uC0 + Perturbations1;
% Perturbations2 = 1+0.1*randn(nS,1);
% Perturbations1 = 1+0.1*randn(nC,1);
% uS=uS0 * Perturbations2;
% uC=uC0 * Perturbations1;

% Define cubic equation in A
eqn = @(A) b*HC * etaS(e) * (A*csh + a)^3 - (A.^2) * D_S * (ssh^2);

% Solve numerically
A_sol = fzero(eqn, A_sol_prev);
A_sol_prev = A_sol;   % continuation method

fprintf('Solution for A = %.6f\n', A_sol);

C2 = sqrt(D_S)*(A_sol*ssh) / (-etaS(e)*(A_sol*csh + a)^2) + b*HC^2/(2*D_C);
fprintf('C2 = %.6f\n', C2);
uS0 = A_sol*cosh((pS(2,:)' - (HC+HS))/sqrt(D_S)) + a;
uC0 = - (b/(2*D_C))*pC(2,:)'.^2 + C2;

eps = 1e-2;
pert_interior_C = eps * randn(nC,1);   % excludes endpoints y(1) and y(end)
pert_interior_S = eps * randn(nS,1);

uC = uC0 + pert_interior_C;
uS = uS0 + pert_interior_S;

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
ylabel('$y$','Interpreter','latex');
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

   %% % if any(isnan(uS(:)))
% disp(['min(uS) = ', num2str(min(uS(:)))])
% disp(['max(uS) = ', num2str(max(uS(:)))])
% disp(['any NaN in uS: ', num2str(any(isnan(uS(:))))])
%     % end

%% --- Right-hand side assembly ---
% ord=3;
% F_S = ReactKineInt(pS,tS,uS,f_S,ord);
% F_C = ReactKineInt(pC,tC,uC,f_C,ord);
% 
% ord=2;
% bS_bottom = NLBoundFluxInt(bottom_tS, top_tC, pS, pC, uS, uC, G_S, etaS, ord);
% bC_top    = NLBoundFluxInt(top_tC, bottom_tS, pC, pS, uC, uS, G_C, etaC, ord);

ord=3;
F_S = ReactKineInt2D(pS,tS,pC,tC,uS,uC,etaS(e),fc_S,ord);
F_C = ReactKineInt2D(pC,tC,pS,tS,uC,uS,etaC(e),fc_C,ord);

%% Matrix system
AuS = M_S + dt/2*D_S*S_S;
BuS = (M_S - dt/2*D_S*S_S)*uS + dt*F_S;% + dt*bS_bottom;
uS_new = AuS\BuS;
AuC = M_C + dt/2*D_C*S_C;
BuC = (M_C - dt/2*D_C*S_C)*uC + dt*F_C;% + dt*bC_top;
uC_new = AuC\BuC;

%%
% U_big = [uS; uC];
% F_big = [F_S; F_C];
% B_big = [bS_bottom; bC_top];
% A = M_big + 0.5*dt*S_big;
% B = (M_big - 0.5*dt*S_big)*U_big + dt*(F_big + B_big);

%% Stability and NaN check
% if any(isnan(A(:))) || any(isnan(B(:)))
%  [any(isnan(A(:)))   %1
% any(isnan(B(:)))    %2
% any(isnan(uS(:)))   %3
% any(isnan(uC(:)))   %4
% any(isnan(B_big(:)))    %5
% any(isnan(F_big(:)))    %6
% any(isnan(bS_bottom(:)))    %7
% any(isnan(bC_top(:)))]    %8
%    error('Matrix A or RHS B contains NaN at step %d (time = %.5f)', ti, t);
% end
% if condest(A) > 1e12
%     warning('Matrix A is ill-conditioned at step %d (time = %.5f), condest = %.2e', ti, t, condest(A));
% end

%% Solve system
% U = A \ B;
% if any(isnan(U)) || any(isinf(U))
%  [any(isnan(A(:)))   %1
% any(isnan(B(:)))    %2
% any(isnan(uS(:)))   %3
% any(isnan(uB(:)))   %4
% any(isnan(uC(:)))   %5
% any(isnan(B_big(:)))    %6
% any(isnan(F_big(:)))    %7
% any(isnan(bS_bottom(:)))    %8
% any(isnan(bB_bottom(:)))    %9
% any(isnan(bB_top(:)))    %10
% any(isnan(bC_top(:)))]    %11
%     error('Solution exploded at step %d (time = %.5f)', ti, t);
% end
% if length(U) ~= N1+N2+N3
%     error('Inconsistent solution vector size!');
% end
% uS_new = U(1:N3);
% uB_new = U(N3+1:N3+N2);
% uC_new = U(N3+N2+1:end);

    %% Pattern detection using norms
    errS = norm(uS_new - uS)/norm(uS);
    errC = norm(uC_new - uC)/norm(uC);
    maxErr = max([errS, errC]);

%% Update and save
uS = uS_new;
uC = uC_new;

% uS_store(:,ti) = uS;
% uB_store(:,ti) = uB;
% uC_store(:,ti) = uC;
% toc
    %% Plot and gif
    if mod(ti, drawperframe) == 1
        vqS = griddata(pS(1,:),pS(2,:),uS,xqS,yqS);
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

%% --- Space-Time Visualization ---
% figure;
% subplot(3,1,1);
% imagesc(tt(1:ti), 1:N3, uS_store(:,1:ti));
% xlabel('Time'); ylabel('Node Index'); title('u_S Space-Time'); colorbar;
% 
% subplot(3,1,2);
% imagesc(tt(1:ti), 1:N2, uB_store(:,1:ti));
% xlabel('Time'); ylabel('Node Index'); title('u_B Space-Time'); colorbar;
% 
% subplot(3,1,3);
% imagesc(tt(1:ti), 1:N1, uC_store(:,1:ti));
% xlabel('Time'); ylabel('Node Index'); title('u_C Space-Time'); colorbar;

%%
fprintf('\nDone!\n');
toc
diary OFF