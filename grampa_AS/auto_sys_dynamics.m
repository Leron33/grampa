clear all;
load('.\mat_files\auto_sys_data.mat');

n = 10000;
deg = sum(auto_sys_mat{1});
[~, I] = maxk(deg, n);
A = auto_sys_mat{1}(I, I);

corr_sp = zeros(9, 1);
val_sp = zeros(9, 1);
corr_dp = zeros(9, 1);
val_dp = zeros(9, 1);
% corr_qp = zeros(9, 1);
% val_qp = zeros(9, 1);
val_truth = zeros(9, 1);

for i = 1:9
    disp(i);
    B = auto_sys_mat{i}(I, I);
    val_truth(i) = sum(dot(A, B));
    
    P_rnd = sparse(eye(n));
    P_rnd = P_rnd(:, randperm(n));
    B = P_rnd * B * P_rnd';

    %% Robust spectral
    tic;
    P_sp = sparse(matching_robust_spectral(full(A), full(B), 1));
%     P_sp = sparse(matching_robust_spectral_sparse(A, B, 1));
%     A1 = full(A); A1 = A1 - mean(mean(A1)); 
%     B1 = full(B); B1 = B1 - mean(mean(B1)); 
%     P_sp = sparse(matching_robust_spectral(A1, B1, 1));
    corr_sp(i) = sum(dot(P_rnd, P_sp))/n;
    val_sp(i) = sum(dot(P_sp * A * P_sp', B));
    toc;
    
    %% Degree profile
    tic;
    P_dp = sparse(matching_deg_pro(full(A), full(B)));
    corr_dp(i) = sum(dot(P_rnd, P_dp))/n;
    val_dp(i) = sum(dot(P_dp * A * P_dp', B));
    toc;
    
    %% Full constrained quadratic programming
%     tic;
%     P_qp = matching_full_qp(A, B);
%     corr_qp(i) = sum(dot(P_rnd, P_qp))/n;
%     val_qp(i) = sum(dot(P_qp * A * P_qp', B));
%     toc;
    
    save(strcat('.\mat_files\auto_sys_dynamics_', int2str(i), '.mat'),  'corr_sp', 'corr_dp', 'val_sp', 'val_dp', 'val_truth', 'B', 'P_rnd', 'P_sp', 'P_dp');
    clear B P_rnd P_sp P_dp;
end

save('.\mat_files\auto_sys_dynamics.mat');
% save('.\mat_files\auto_sys_dynamics.mat', 'corr_sp', 'corr_dp', 'val_sp', 'val_dp', 'val_truth');