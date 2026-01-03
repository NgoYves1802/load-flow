function ga_ieee14_demo()
clc; clear;
%% -------------------- Network Data (p.u.) --------------------
% Bus: 1 slack, 2,3,6,8 PV, others PQ
% V1 fixed:
V1 = 1.06;
% Loads (Pd, Qd) in p.u.
Pd = zeros(14,1); Qd = zeros(14,1);
Pd(2)=0.217; Qd(2)=0.127;
Pd(3)=0.942; Qd(3)=0.19;
Pd(4)=0.478; Qd(4)=-0.039;
Pd(5)=0.076; Qd(5)=0.016;
Pd(6)=0.112; Qd(6)=0.075;
Pd(9)=0.295; Qd(9)=0.166;
Pd(10)=0.09; Qd(10)=0.058;
Pd(11)=0.035; Qd(11)=0.018;
Pd(12)=0.061; Qd(12)=0.016;
Pd(13)=0.135; Qd(13)=0.058;
Pd(14)=0.149; Qd(14)=0.05;
% Active generation at PV buses (p.u.)
Pg = zeros(14,1);
Pg(2) = 0.40;
Pg(3) = 0.00;
Pg(6) = 0.00;
Pg(8) = 0.00;
% Qg limits for PV buses
Qg_min = zeros(14,1); Qg_max = zeros(14,1);
Qg_min(2) = -0.40; Qg_max(2) = 0.50;
Qg_min(3) = 0.00; Qg_max(3) = 0.40;
Qg_min(6) = -0.06; Qg_max(6) = 0.24;
Qg_min(8) = -0.06; Qg_max(8) = 0.24;
% Lines and transformers: [from to R X B] (B total charging)
lines = [
    1 2 0.01938 0.05917 0.05280
    1 5 0.05403 0.22304 0.04920
    2 3 0.04699 0.19797 0.04380
    2 4 0.05811 0.17632 0.03400
    2 5 0.05695 0.17388 0.03460
    3 4 0.06701 0.17103 0.01280
    4 5 0.01335 0.04211 0.00000
    4 7 0.00000 0.20912 0.00000
    4 9 0.00000 0.55618 0.00000
    5 6 0.00000 0.25202 0.00000
    6 11 0.09498 0.19890 0.00000
    6 12 0.12291 0.25581 0.00000
    6 13 0.06615 0.13027 0.00000
    7 8 0.00000 0.17615 0.00000
    7 9 0.00000 0.11001 0.00000
    9 10 0.03181 0.08450 0.00000
    9 14 0.12711 0.27038 0.00000
    10 11 0.08205 0.19207 0.00000
    12 13 0.22092 0.19988 0.00000
    13 14 0.17093 0.34802 0.00000
];
Ybus = buildYbus(14, lines);
% Add shunt at bus 9
Ybus(9,9) = Ybus(9,9) + 1j*0.19;
% Net injection specifications (P,Q) = Pg - Pd
Pspec = zeros(14,1); Qspec = zeros(14,1);
pv_buses = [2 3 6 8];
pq_buses = setdiff(2:14, pv_buses);
% Set Pspec for slack and PV buses
Pspec(1) = Pg(1) - Pd(1); % Will be determined by power balance
for i = pv_buses
    Pspec(i) = Pg(i) - Pd(i);
end
% Set P and Q spec for PQ buses
for i = pq_buses
    Pspec(i) = -Pd(i);
    Qspec(i) = -Qd(i);
end
%% -------------------- GA Parameters --------------------
Npop = 40;
Ng = 8000;
Pc = 0.85; % crossover probability
Pm = 0.20; % mutation probability (per gene)
sigma_mut = 0.01; % mutation std dev
elite = 1;
Vmin = 0.95; Vmax = 1.05;
% Penalties
alpha = 1e4; % for delta P
beta = 1e4; % for delta Q
gamma = 1e4; % for Qg violations
% Chromosome: [V2 V3 ... V14]
D = 13;
%% -------------------- Initialization --------------------
pop = Vmin + (Vmax - Vmin)*rand(Npop, D);
fit = zeros(Npop,1);
bestFit = inf; bestX = []; bestInfo = struct();
%% -------------------- GA Loop --------------------
for gen = 1:Ng
    % Evaluate fitness
    for i = 1:Npop
        [fit(i), info] = fitnessIEEE14(pop(i,:), V1, Ybus, Pspec, Qspec, Qg_min, Qg_max, alpha, beta, gamma, lines, pv_buses, Qd);
        if fit(i) < bestFit
            bestFit = fit(i);
            bestX = pop(i,:);
            bestInfo = info;
        end
    end
    % Display progress
    fprintf('Gen %3d | bestFit = %.6f | Ploss = %.6f\n', gen, bestFit, bestInfo.Ploss);
    % New population
    newPop = zeros(size(pop));
    % Elitism
    if elite
        newPop(1,:) = bestX;
        startIdx = 2;
    else
        startIdx = 1;
    end
    % Fill population
    k = startIdx;
    while k <= Npop
        % Tournament selection
        p1 = pop(tournament(fit),:);
        p2 = pop(tournament(fit),:);
        % Crossover
        if rand < Pc
            [c1,c2] = arithmeticCrossover(p1,p2);
        else
            c1 = p1; c2 = p2;
        end
        % Mutation
        c1 = gaussianMutation(c1, Pm, sigma_mut);
        c2 = gaussianMutation(c2, Pm, sigma_mut);
        % Enforce bounds
        c1 = min(max(c1, Vmin), Vmax);
        c2 = min(max(c2, Vmin), Vmax);
        newPop(k,:) = c1;
        if k+1 <= Npop
            newPop(k+1,:) = c2;
        end
        k = k + 2;
    end
    pop = newPop;
end
%% -------------------- Results --------------------
fprintf('\n=== Best solution found ===\n');
fprintf('V2 to V14: ');
fprintf('%.4f ', bestX);
fprintf('\n');
fprintf('Ploss (p.u.) = %.6f\n', bestInfo.Ploss);
fprintf('Angles (rad):\n');
disp(bestInfo.theta');
fprintf('P mismatches (excluding slack):\n');
disp(bestInfo.dP(2:end)');
fprintf('Q mismatches (PQ buses only):\n');
disp(bestInfo.dQ');
end
%% ==================== FITNESS ====================
function [F, info] = fitnessIEEE14(X, V1, Ybus, Pspec, Qspec, Qg_min, Qg_max, alpha, beta, gamma, lines, pv_buses, Qd)
% X = [V2 V3 ... V14]
V = [V1; X(:)];
nb = 14;
% Estimate angles theta2..theta14 by solving P equations (Newton)
theta0 = zeros(nb,1);
theta = solveAnglesNewtonP(V, theta0, Ybus, Pspec);
% Calculate P,Q from (V,theta)
[Pcalc, Qcalc] = calcPQ(V, theta, Ybus);
% Active losses
Ploss = calcLosses(V, theta, lines);
% Mismatches
dP = Pspec - Pcalc;
% slack (bus 1) not penalized
dP(1) = 0;
% Q mismatch only on PQ buses
dQ = Qspec - Qcalc;
dQ(1) = 0; 
dQ(pv_buses) = 0; % PV buses controlled
% Qg for PV buses = Qcalc(pv) + Qd(pv)
penQg = 0;
for pv = pv_buses
    Qg = Qcalc(pv) + Qd(pv);
    if Qg < Qg_min(pv)
        penQg = penQg + (Qg_min(pv) - Qg)^2;
    elseif Qg > Qg_max(pv)
        penQg = penQg + (Qg - Qg_max(pv))^2;
    end
end
F = Ploss + alpha*sum(dP.^2) + beta*sum(dQ.^2) + gamma*penQg;
info.Ploss = Ploss;
info.theta = theta;
info.Pcalc = Pcalc; 
info.Qcalc = Qcalc;
info.dP = dP; 
info.dQ = dQ; 
end
%% ==================== NEWTON on P equations ====================
function theta = solveAnglesNewtonP(V, theta0, Ybus, Pspec)
% Solves P(V,theta) = Pspec for theta2..theta14 (theta1=0)
nb = length(V);
theta = theta0;
theta(1) = 0;
maxIt = 30;
tol = 1e-8;
pqAngles = (2:nb); % variables theta2..theta14
for it = 1:maxIt
    [Pcalc, ~] = calcPQ(V, theta, Ybus);
    mismatch = Pspec(pqAngles) - Pcalc(pqAngles);
    if norm(mismatch, inf) < tol
        return;
    end
    % Jacobian dP/dtheta (dimension 13x13)
    J = jacobianP_theta(V, theta, Ybus, pqAngles);
    % Check for singularity
    if rcond(J) < 1e-14
        warning('Jacobian is nearly singular');
        return;
    end
    % Update
    dth = J \ mismatch;
    theta(pqAngles) = theta(pqAngles) + dth;
end
% If not converged, return last value (GA will penalize)
end
function J = jacobianP_theta(V, theta, Ybus, idx)
% J(m,n) = dP_i/dtheta_k for i=idx(m), k=idx(n)
nb = length(V);
G = real(Ybus); 
B = imag(Ybus);
m = length(idx);
J = zeros(m,m);
% Pre-calculate Q for diagonal terms
[~, Q] = calcPQ(V, theta, Ybus);
for a = 1:m
    i = idx(a);
    for b = 1:m
        k = idx(b);
        if i == k
            % dP_i/dtheta_i = -Q_i - B_ii*V_i^2
            J(a,b) = -Q(i) - B(i,i)*V(i)^2;
        else
            % dP_i/dtheta_k = V_i V_k (G_ik sin(theta_i-theta_k) - B_ik cos(theta_i-theta_k))
            thik = theta(i) - theta(k);
            J(a,b) = V(i)*V(k)*( G(i,k)*sin(thik) - B(i,k)*cos(thik) );
        end
    end
end
end
%% ==================== P,Q Calculation ====================
function [P, Q] = calcPQ(V, theta, Ybus)
nb = length(V);
G = real(Ybus); 
B = imag(Ybus);
P = zeros(nb,1); 
Q = zeros(nb,1);
for i = 1:nb
    for k = 1:nb
        th = theta(i) - theta(k);
        P(i) = P(i) + V(i)*V(k)*( G(i,k)*cos(th) + B(i,k)*sin(th) );
        Q(i) = Q(i) + V(i)*V(k)*( G(i,k)*sin(th) - B(i,k)*cos(th) );
    end
end
end
%% ==================== Losses ====================
function Ploss = calcLosses(V, theta, lines)
Ploss = 0;
for e = 1:size(lines,1)
    i = lines(e,1); 
    j = lines(e,2);
    R = lines(e,3); 
    X = lines(e,4);
    % B not used for losses
    y = 1/(R + 1i*X);
    g = real(y);
    th = theta(i) - theta(j);
    Ploss = Ploss + g*( V(i)^2 + V(j)^2 - 2*V(i)*V(j)*cos(th) );
end
end
%% ==================== Ybus ====================
function Ybus = buildYbus(nb, lines)
Ybus = zeros(nb,nb);
for e = 1:size(lines,1)
    i = lines(e,1); 
    j = lines(e,2);
    R = lines(e,3); 
    X = lines(e,4); 
    B = lines(e,5);
    z = R + 1i*X;
    if abs(z) > 0
        y = 1 / z;
    else
        y = 0;
    end
    y_sh = 1i * B / 2;
    Ybus(i,i) = Ybus(i,i) + y + y_sh;
    Ybus(j,j) = Ybus(j,j) + y + y_sh;
    Ybus(i,j) = Ybus(i,j) - y;
    Ybus(j,i) = Ybus(j,i) - y;
end
end
%% ==================== GA Operators ====================
function idx = tournament(fit)
% Tournament selection of size 2
n = length(fit);
a = randi(n); 
b = randi(n);
if fit(a) < fit(b)
    idx = a; 
else
    idx = b; 
end
end
function [c1,c2] = arithmeticCrossover(p1,p2)
% Arithmetic crossover
lam = rand;
c1 = lam*p1 + (1-lam)*p2;
c2 = (1-lam)*p1 + lam*p2;
end
function c = gaussianMutation(c, Pm, sigma)
for g = 1:length(c)
    if rand < Pm
        c(g) = c(g) + sigma*randn;
    end
end
end