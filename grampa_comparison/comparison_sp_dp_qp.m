% Compare gradient descent methods for graph matching 
clear all;

%% initialization
p = 0.5;
num_run = 10;
vec_dim = 500:10:500;
len_dim = length(vec_dim);
vec_noise = 0:0.05:0.8;
len_noise = length(vec_noise);

full_qp_corr = ones(len_dim, len_noise, num_run);
full_qp_run = zeros(len_dim, len_noise, num_run);
deg_pro_corr = ones(len_dim, len_noise, num_run);
deg_pro_run = zeros(len_dim, len_noise, num_run);
robust_corr = ones(len_dim, len_noise, num_run);
robust_run = zeros(len_dim, len_noise, num_run);

%% Iteration over independent samples 
for ind_run = 1:num_run
    fprintf('Iteration %i \n', ind_run);

    %% Iteration over dimensions 
    for ind_dim = 1:len_dim
        n = vec_dim(ind_dim);
        fprintf('Matrix dimension %i \n', n);
        
        %% Iteration over noise levels 
        for ind_noise = 3:len_noise
            sigma = vec_noise(ind_noise); disp(sigma);
            [A, B, A0, B0, P_rnd] = generate_er(n, p, sigma);
            
            %% Full QP 
            tic;
            P = matching_full_qp(A0, B0);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            full_qp_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            full_qp_run(ind_dim, ind_noise, ind_run) = toc;
            
            %% Degree profile
            tic;
            P = matching_deg_pro(A0, B0);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            deg_pro_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            deg_pro_run(ind_dim, ind_noise, ind_run) = toc;
                  
            %% Robust spectral method 
            tic;
            P = matching_robust_spectral(A, B, 0.2);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            robust_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            robust_run(ind_dim, ind_noise, ind_run) = toc;
        end
    end
end

full_qp_corr_mean = mean(full_qp_corr, 3);
deg_pro_corr_mean = mean(deg_pro_corr, 3);
robust_corr_mean = mean(robust_corr, 3);

full_qp_run_mean = mean(mean(mean(full_qp_run(1, 3:end, :))));
deg_pro_run_mean = mean(mean(mean(deg_pro_run(1, 3:end, :))));
robust_run_mean = mean(mean(mean(robust_run(1, 3:end, :))));

disp([robust_run_mean, deg_pro_run_mean, full_qp_run_mean]);

clear -regexp _corr$ _run$;
save('.\mat_files\comparison_sp_dp_qp.mat');