clear all; close all; clc;

%%
syms yS yC

LS=[0.2;0.1]; 
LC=[0.1;0];

HS=abs(LS(1)-LS(2)); HC=abs(LC(1)-LC(2));
a = 2.1; b=1.1; etaS = [0.01 0.05 0.1 0.5 1 5 9 10];%0.1;
etaC=etaS;
D_S=1; D_C=20;

% for e=1:length(etaS)
% A = -b*HC/(sqrt(D_S)*sinh(HS/sqrt(D_S)));
% uS0 = A*cosh((yS - (HC+HS))/sqrt(D_S)) + a;
% C2 = b*HC/(etaS(e)*( (-b*HC/sqrt(D_S))*coth(HS/sqrt(D_S)) + a )^2) + (b/(2*D_C))*HC^2;
% uC0 = - (b/(2*D_C))*yC.^2 + C2;
% 
% % [A uS0 C2 uC0]
% [A C2]
% [subs(uS0,yS,0:0.1:1)' subs(uC0,yC,0:0.1:1)']
% end

% Constants
csh = cosh(HS/sqrt(D_S));
ssh = sinh(HS/sqrt(D_S));
A_sol_prev = 1e-6;

for e=1:length(etaS)
    evalS=round(etaS(e),4);
evalC=round(etaC(e),4);
fprintf('\nFor etaS=%.4f, etaC=%.4f,\n',evalS, evalC);
eval=['_',num2str(evalS),'_',num2str(evalC)];

% Define cubic equation in A
eqn = @(A) b*HC * etaS(e) * (A*csh + a)^3 - (A.^2) * D_S * (ssh^2);

% Solve numerically
A_sol = fzero(eqn, A_sol_prev);
A_sol_prev = A_sol;   % continuation method

fprintf('Solution for A = %.6f\n', A_sol);

C2 = sqrt(D_S)*(A_sol*ssh) / (-etaS(e)*(A_sol*csh + a)^2) + b*HC^2/(2*D_C);
fprintf('C2 = %.6f\n', C2);
uS0 = A_sol*cosh((yS - (HC+HS))/sqrt(D_S)) + a;
uC0 = - (b/(2*D_C))*yC.^2 + C2;

[double(subs(uS0,yS,0:0.01:0.1))' double(subs(uC0,yC,0:0.01:0.1))']
end
