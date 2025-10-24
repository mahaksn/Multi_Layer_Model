clear all; close all; clc;
tic

%% pre-setting
showanimation=1;
makegif=1;
drawperframe=500;
pattern_detected = false;
tol=1e-7;
T=1000;

%% Load Mesh Data: p-nodes; t-elements;
LS=[0.2;0.1]; 
LC=[0.1;0];
Ne = 200;

% Node coordinates
pS = linspace(LS(1), LS(2), Ne + 1)';
pC = linspace(LC(1), LC(2), Ne + 1)';
% Element connectivity
tS = [1:Ne; 2:Ne+1]'; 
tC = [1:Ne; 2:Ne+1]'; 

nC=length(pC);
nS=length(pS);

%% Parameters
% a = 0.1; b=1; etaS = [0.01 0.05 0.1 0.5 1 5 10];
% etaC=etaS;
% D_S=0.1; D_C=10;

% Hs=abs(LS(1)-LS(2)); Hc=abs(LC(1)-LC(2)); k=2;
% D_S=1; D_C=1;
% etaS=[6 10 20]; etaC=[6 10 20];

% a=-1; b=-3/2; H=2; D=0.08; D_S=D; D_C=1;
% C=1; etaS=1; etaC=etaS;
% Del=5*C^2-8*H+12;

HS=abs(LS(1)-LS(2)); HC=abs(LC(1)-LC(2));
a = 2.1; b=1.075; etaS = 0.1;%[0.01 0.05 0.1 0.5 1 5 10];
etaC=etaS;
D_S=1; D_C=20;

%% Steady states
% uS0=a; uC0=b/a;
% uS0=(-5*C+sqrt(5*Del))/(4*H+4);
% uC0=(-5*C+sqrt(5*Del))/(10);

% uS0=a+b; uC0=b/(a+b)^2;

%% Define the Source Term and Coupling Term
% f_S=@(u,v,eta) 0*u + 0*v;
% f_C=@(u,v,eta) 0*u + 0*v;
% 
% G_S=@(us,uc) us*uc;
% G_C=@(uc,us) us*uc;

% f_S=@(u,v,eta) u;
% f_C=@(u,v,eta) b*v;
% 
% G_S = @(u, v) a*v-C*u.*v-u.*v.^2;
% G_C = @(u, v) H*u+C*u.*v+u.*v.^2;

f_S=@(u,v,eta) a - u;
f_C=@(u,v,eta) b + 0*v;

G_S = @(u, v) u.^2.*v;
G_C = @(u, v) u.^2.*v;
% G_C = @(u, v) u.*v.^2;

%% Assemble Matrices
ord=3;
[S_S,M_S]=AssembleGlobalMatrices1D(pS,tS,ord);
[S_C,M_C]=AssembleGlobalMatrices1D(pC,tC,ord);
M_big = blkdiag(M_S,M_C);
S_big = blkdiag(D_S*S_S, D_C*S_C);

%% saving folder and name
folder=['D:\20081\MATLAB Output\Paper-2_Outputs\' ...
    'Paper-2_Multi_layer_domain\FEM_2_y1D\'];
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
dt=0.001;%1e-4;
nt=T/dt+1;
stopti=nt; % stop patterning at nt
nFrame=ceil((T/dt)/drawperframe);

%% Preallocate storage
max_steps = nt;
uS_store = zeros(max_steps,nS);
uC_store = zeros(max_steps,nC);
tt = zeros(max_steps,1); % store time values

%%
for e=1:length(etaS)
evalS=round(etaS(e),4);
evalC=round(etaC(e),4);
fprintf('\nFor etaS=%.4f, etaC=%.4f,\n',evalS, evalC);
eval=['_',num2str(evalS),'_',num2str(evalC)];

pattern_end=10;
stopti=nt; % stop patterning at nt

%% Initial condition
rng(0); % reproducible random seed
% denom = etaS(e) * Hs * Hc;
% fprintf('\nThe value of eta*k*Hc is %.4f\n',etaS(e)*k*Hc);
% B = (etaS(e) * k * Hc - D_C) / denom;
% A = (D_S / D_C) * B;
% uC_analytic = A * pC;                    % on [0,Hc]
% uS_analytic = B * (pS - (Hc + Hs)) + k;  % on [Hc,Hc+Hs]
% eps_pert = 1e-6;
% 
% uC = uC_analytic .* (1 + eps_pert * randn(size(uC_analytic)));
% uS = uS_analytic .* (1 + eps_pert * randn(size(uS_analytic)));

% pit=0; rng(pit); %change random seed
% % Perturbations2 = 0.01*(2*rand(N2,1)-1);
% % Perturbations1 = 0.01*(2*rand(N1,1)-1);
% % uS=uS0 + Perturbations2;
% % uC=uC0 + Perturbations1;
% PerturbationsS = 1+0.1*randn(nS,1);
% PerturbationsC = 1+0.1*randn(nC,1);
% uS=uS0 * PerturbationsS;
% uC=uC0 * PerturbationsC;

A = -b*HC/(sqrt(D_S)*sinh(HS/sqrt(D_S)));
uS0 = A*cosh((pS - (HC+HS))/sqrt(D_S)) + a;
C2 = b*HC/(etaS(e)*( (-b*HC/sqrt(D_S))*coth(HS/sqrt(D_S)) + a )^2) + b/(2*D_C)*HC^2;
uC0 = - (b/(2*D_C))*pC.^2 + C2;

eps = 1e-4;
pert_interior_C = eps * randn(nC,1);
pert_interior_S = eps * randn(nS,1);
uC = uC0 + pert_interior_C;
uS = uS0 + pert_interior_S;
% 
% % % enforce no-flux at boundaries by reflecting the adjacent interior value
% % % (ensures discrete derivative ~ 0 at boundary)
% % uC(1)       = uC(2);         % duC/dy ~ (uC(2)-uC(1))/dy = 0
% % uC(end)     = uC(end-1);
% % uS(1)       = uS(2);
% % uS(end)     = uS(end-1);

%% set up figure
giffile = [prefix,eval,'_pattern','.gif'];
fig = figure('Color','w');%,'WindowState', 'maximized');

subplot(2,1,1)
hold on
uS_fig=plot(pS, uS, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
xlim([min(LS),max(LS)]);
% ylim([]);
xlabel('$y$','Interpreter','latex');
ylabel('$u_s$','Interpreter','latex');
figtitle=title('t=0');
ax = gca; 
ax.FontSize = 14;
hold off  

subplot(2,1,2)
hold on
uC_fig=plot(pC, uC, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
xlim([min(LC),max(LC)]);
% ylim([]);
xlabel('$y$','Interpreter','latex');
ylabel('$u_c$','Interpreter','latex');
ax = gca; 
ax.FontSize = 14;
hold off  

%% Simulation iteration
for ti = 1:nt
    t = dt * (ti - 1);
    tt(ti) = t;

%% --- Right-hand side assembly ---
ord=2;
% Compute RHS from reactions
F_S = ReactKineInt1D(pS,tS,uS,uC,etaS(e),f_S,ord);
F_C = ReactKineInt1D(pC,tC,uS,uC,etaC(e),f_C,ord);
% F_S(end)
% F_C(1)
% Add Neumann boundary contributions (flux)
g_valS = G_S(uS(end), uC(1));
g_valC = G_C(uS(end), uC(1));
F_S(end) = F_S(end) - etaS(e) * g_valS;
F_C(1)   = F_C(1)   + etaC(e) * g_valC;
% F_S(end)
% F_C(1)

%%
% AuS = M_S + dt/2*D_S*S_S;
% BuS = (M_S - dt/2*D_S*S_S)*uS + dt*(F_S);
% uS_new = AuS\BuS;
% AuB = M_B + dt/2*D_B*S_B;
% BuB = (M_B - dt/2*D_B*S_B)*uB + dt*(F_B);
% uB_new = AuB\BuB;
% AuC = M_C + dt/2*D_C*S_C;
% BuC = (M_C - dt/2*D_S*S_S)*uC + dt*(F_C);
% uC_new = AuC\BuC;

%% Matrix system
U_big = [uS; uC];
% F_big = [F_S; F_C];
A = M_big + 0.5*dt*S_big;
B = (M_big - 0.5*dt*S_big)*U_big; % + dt*(F_big);

%% Stability and NaN check
if any(isnan(A(:))) || any(isnan(B(:)))
 [any(isnan(A(:)))   %1
any(isnan(B(:)))    %2
any(isnan(uS(:)))   %3
any(isnan(uC(:)))   %5
any(isnan(B_big(:)))    %6
any(isnan(F_big(:)))    %7
any(isnan(bS_bottom(:)))    %8
any(isnan(bC_top(:)))]    %11
    error('Matrix A or RHS B contains NaN at step %d (time = %.5f)', ti, t);
end

if condest(A) > 1e12
    warning('Matrix A is ill-conditioned at step %d (time = %.5f), condest = %.2e', ti, t, condest(A));
end

%% Solve system
U = A \ B;
if any(isnan(U)) || any(isinf(U))
    error('Solution exploded at step %d (time = %.5f)', ti, t);
end
if length(U) ~= nC+nS
    error('Inconsistent solution vector size!');
end
uS_new = U(1:nS);
uC_new = U(nS+1:end);

    %% Pattern detection using norms
    errS = norm(uS_new - uS)/norm(uS);
    errC = norm(uC_new - uC)/norm(uC);
    maxErr = max([errS, errC]);

    %% Update and save
uS = uS_new;
uC = uC_new;

uS_store(ti,:) = uS;
uC_store(ti,:) = uC;
% toc
    %% Plot and gif
    if mod(ti, drawperframe) == 1
        if showanimation
            uS_fig.YData = uS;
            uC_fig.YData = uC;
            figtitle.String = ['t = ', num2str(t, '%.2f')];
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
%%
    if ~pattern_detected && maxErr < tol
        fprintf('\nPattern fully developed at t = %.2f\n', t);
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

%%
% % pointwise differences
% err_uC = uC - uC_analytic;
% err_uS = uS - uS_analytic;
% 
% % norms
% L2_uC = norm(err_uC) / sqrt(length(err_uC));
% L2_uS = norm(err_uS) / sqrt(length(err_uS));
% Linf_uC = max(abs(err_uC));
% Linf_uS = max(abs(err_uS));
% 
% % relative norms (scale by analytic magnitude)
% relL2_uC = L2_uC / (norm(uC_analytic)/sqrt(length(uC_analytic)) + eps);
% relL2_uS = L2_uS / (norm(uS_analytic)/sqrt(length(uS_analytic)) + eps);
% 
% fprintf('L2 rel error uC = %.3e, uS = %.3e\n', relL2_uC, relL2_uS);
% fprintf('Linf error uC = %.3e, uS = %.3e\n', Linf_uC, Linf_uS);

%% --- Space-Time Visualization ---
fig_xt=figure;
subplot(1,2,1);
imagesc(1:nS, tt(1:ti), uS_store(1:ti,:));
ylabel('Time'); xlabel('Node Index'); title('u_S Space-Time'); colorbar;

subplot(1,2,2);
imagesc(1:nC, tt(1:ti), uC_store(1:ti,:));
ylabel('Time'); xlabel('Node Index'); title('u_C Space-Time'); colorbar;

%% saving final pattern
saveas(fig_xt,[prefix,eval,'_timeplot.png']);
saveas(fig_xt,[prefix,eval,'_timeplot.fig']);

end
%%
toc
diary OFF

% [uC_analytic uC uC_analytic .* (1 + eps_pert * randn(size(uC_analytic)))]
% 
% [uS_analytic uS uS_analytic .* (1 + eps_pert * randn(size(uS_analytic)))]

fprintf('\nDone!\n');
