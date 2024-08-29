clear; clc;
yalmip('clear');

%% set up the initial data for the problem
K = [-1 1 0; 2 -3 1; 3 0 -3];
alpha = [1; 0.5; -1];
A = diag(sign(alpha)) * K * diag(alpha);

%%
for tau= 1e-6:0.01:1 % create a loop for the time-delayed of the problem

% variable declaration
P = sdpvar(3,3,'symmetric');
Z = sdpvar(3,3,'symmetric');
Q = sdpvar(3,3, 'symmetric');

% set up LMIs
M=[P+Z, -Z; -Z, tau*Q+Z];
N= [ Q- 1/tau *Z, P*A+1/tau *Z, zeros(3); A'*P+1/tau * Z, -Q-1/tau * Z, tau*A'*Z; zeros(3), tau*Z*A, -tau*Z];

Constraints = [P>= eye(3)*1e-9, Z>= eye(3)*1e-9, M>=eye(6)*1e-9, N<= -1e-9 *eye(9)];

sol = optimize(Constraints);

if(sol.problem ~= 0)
    break
end
end

%% 
disp(value(tau-0.01))
