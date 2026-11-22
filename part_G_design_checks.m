%% Part G: Design checks - R choice, gate confidence level, bias random-walk level
% Extra runs used to justify the design choices in the report. The baseline
% filter (Parts C-E) is not changed here. Outlier set for scoring = the
% true outliers identified in Part D (TRUTH USED FOR EVALUATION ONLY).

known_outliers = false(1, M);
known_outliers(true_outlier_idx) = true;

%% (a) R from point estimates vs. upper 95% bounds
R_point = diag([sigma_n^2, sigma_e^2]);
P0_point = P0;
P0_point(1,1) = sigma_n^2 / n_gnss_static;
P0_point(2,2) = sigma_e^2 / n_gnss_static;

out_Rpoint = ekf6_run(t_imu, u_imu, t_gnss, z_gnss, x0, P0_point, Q, R_point, true, gamma);

fprintf('\n--- Part G(a): R from point estimates vs. upper 95%% bounds (gated) ---\n');
fprintf('%14s %10s %14s %10s %10s %10s %8s %8s\n', 'R', 'RMSE [m]', 'MaxOutErr[m]', 'MeanNIS', 'MeanNEES', 'NIS in95', '#Rej', '#FA');
print_check_row('point est.', out_Rpoint, x_true_all, eval_mask, outage_mask, t_gnss, t_imu, known_outliers);
print_check_row('upper 95%', out_Dgated, x_true_all, eval_mask, outage_mask, t_gnss, t_imu, known_outliers);

%% (b) Gate confidence level
p_list = [0.95, 0.99, 0.999];
fprintf('\n--- Part G(b): gate confidence level ---\n');
fprintf('%14s %10s %14s %10s %10s %10s %8s %8s\n', 'p (gamma)', 'RMSE [m]', 'MaxOutErr[m]', 'MeanNIS', 'MeanNEES', 'NIS in95', '#Rej', '#FA');
for pi_ = 1:numel(p_list)
    gamma_p = -2*log(1 - p_list(pi_));
    out_p = ekf6_run(t_imu, u_imu, t_gnss, z_gnss, x0, P0, Q, R, true, gamma_p);
    print_check_row(sprintf('%.3f (%.2f)', p_list(pi_), gamma_p), out_p, x_true_all, eval_mask, outage_mask, t_gnss, t_imu, known_outliers);
end

% How loose can the gate be before the weakest outlier gets through?
d2_min_outlier = min(out_Dgated.d2(true_outlier_idx));
p_pass = 1 - exp(-d2_min_outlier/2);
fprintf('Smallest outlier d^2 = %.2f -> it would pass for gamma > %.2f, i.e. p > %.5f\n', ...
    d2_min_outlier, d2_min_outlier, p_pass);

%% (c) Bias random-walk level (DRIFT_FRACTION)
drift_list = [0.1, 0.5, 2];
fprintf('\n--- Part G(c): bias random-walk level (fraction of sigma over the mission) ---\n');
fprintf('%14s %10s %14s %10s %10s %10s %8s %8s\n', 'fraction', 'RMSE [m]', 'MaxOutErr[m]', 'MeanNIS', 'MeanNEES', 'NIS in95', '#Rej', '#FA');
for di = 1:numel(drift_list)
    Qd = Q;
    Qd(5,5) = (drift_list(di)*sigma_a)^2 / N;
    Qd(6,6) = (drift_list(di)*sigma_w)^2 / N;
    out_d = ekf6_run(t_imu, u_imu, t_gnss, z_gnss, x0, P0, Qd, R, true, gamma);
    print_check_row(sprintf('%.1f', drift_list(di)), out_d, x_true_all, eval_mask, outage_mask, t_gnss, t_imu, known_outliers);
end

%% --- local function ---
function print_check_row(label, out, x_true_all, eval_mask, outage_mask, t_gnss, t_imu, known_outliers)
    err_h = hypot(out.X(1,:) - x_true_all(1,:), out.X(2,:) - x_true_all(2,:));
    rmse = sqrt(mean(err_h(eval_mask).^2));
    max_out = max(err_h(outage_mask));

    err = out.X - x_true_all;
    err(3,:) = wrap_angle(err(3,:));
    nees = nan(1, size(err,2));
    for k = find(eval_mask)
        nees(k) = err(:,k)' * (out.Pfull(:,:,k) \ err(:,k));
    end
    nis_vals = out.d2(out.accepted & t_gnss >= 20);
    nis_in95 = mean(nis_vals >= chi2inv(0.025,2) & nis_vals <= chi2inv(0.975,2));

    rejected = ~out.accepted;
    n_fa = sum(rejected & ~known_outliers);
    fprintf('%14s %10.3f %14.3f %10.3f %10.3f %9.1f%% %8d %8d\n', label, rmse, max_out, ...
        mean(nis_vals), mean(nees(eval_mask)), 100*nis_in95, sum(rejected), n_fa);
end
