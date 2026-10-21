%% ============================================================
% Compute:
% 1) Fleiss' kappa
% 2) proportion of full disagreement (3 raters all different)
%% ============================================================

clear; clc;

file = 'D:\online_learning\production\all_rater.xlsx';

%% -------- read data --------
T = readtable(file,'ReadVariableNames',false);

r1 = lower(string(T{:,2}));
r2 = lower(string(T{:,3}));
r3 = lower(string(T{:,4}));

ratings = [r1 r2 r3];

N = size(ratings,1);
n = 3; % raters

%% -------- category mapping --------
cats = unique(ratings(:));
cats(cats=="") = [];

K = length(cats);

R = zeros(N,n);

for i = 1:N
    for j = 1:n
        R(i,j) = find(cats==ratings(i,j));
    end
end

%% -------- Fleiss kappa --------
M = zeros(N,K);

for i = 1:N
    for j = 1:n
        M(i,R(i,j)) = M(i,R(i,j)) + 1;
    end
end

p = sum(M)/(N*n);
P = (sum(M.^2,2)-n)/(n*(n-1));

Pbar = mean(P);
Pe = sum(p.^2);

fleiss_kappa = (Pbar-Pe)/(1-Pe);

%% -------- full disagreement --------
all_diff = false(N,1);

for i = 1:N
    all_diff(i) = (length(unique(R(i,:)))==3);
end

prop_all_diff = mean(all_diff);

%% -------- print results --------
fprintf('\n================ RESULTS ================\n');
fprintf('Total tokens                  : %d\n',N);
fprintf('Fleiss'' kappa                 : %.4f\n',fleiss_kappa);
fprintf('Fully disagreed tokens        : %d\n',sum(all_diff));
fprintf('Proportion fully disagreed    : %.4f (%.2f%%)\n',...
        prop_all_diff,100*prop_all_diff);