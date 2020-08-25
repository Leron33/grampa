% Compare spectral methods for graph matching 
clear all;

%% initialization
p = 0.5;
num_run = 10;
vec_dim = 100:10:100;
len_dim = length(vec_dim);
vec_noise = 0:0.05:0.5;
len_noise = length(vec_noise);

topeig_corr = zeros(len_dim, len_noise, num_run);
isorank_corr = zeros(len_dim, len_noise, num_run);
eigalign_corr = zeros(len_dim, len_noise, num_run);
lowrank_corr = zeros(len_dim, len_noise, num_run);
umeyama_corr = zeros(len_dim, len_noise, num_run);
robust_corr = zeros(len_dim, len_noise, num_run);

%% Iteration over independent samples 
for ind_run = 1:num_run
    fprintf('Iteration %i \n', ind_run);

    %% Iteration over dimensions 
    for ind_dim = 1:len_dim
        n = vec_dim(ind_dim);
        fprintf('Matrix dimension %i \n', n);
        
        %% Iteration over noise levels 
        for ind_noise = 1:len_noise
            tic;
            
            sigma = vec_noise(ind_noise); disp(sigma);
            [A, B, A0, B0, P_rnd] = generate_er(n, p, sigma);

            %% Top eigenvector alignment 
            P = matching_top_eigvec(A, B);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            topeig_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            
            %% IsoRank 
            P = matching_isorank(A0, B0, 0.85);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            isorank_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            
            %% EigenAlign
            P = matching_eigenalign(A, B, 0.1/sqrt(n*p*(1-p)));
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            eigalign_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            
            %% LowRankAlign
            P = matching_lowrankalign(A, B, 2);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            lowrank_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;

            %% Umeyama's method 
            P = matching_umeyama(A, B);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            umeyama_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;

            %% Robust spectral method 
            P = matching_robust_spectral(A, B, 0.2);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            robust_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            
            toc;
        end
    end
end

topeig_corr_mean = mean(topeig_corr, 3);
isorank_corr_mean = mean(isorank_corr, 3);
eigalign_corr_mean = mean(eigalign_corr, 3);
lowrank_corr_mean = mean(lowrank_corr, 3);
umeyama_corr_mean = mean(umeyama_corr, 3);
robust_corr_mean = mean(robust_corr, 3);

clear -regexp _corr$;
save('.\mat_files\comparison_spectral_methods.mat');