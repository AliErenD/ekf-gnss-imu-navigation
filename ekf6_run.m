function out = ekf6_run(t_imu, u_imu, t_gnss, z_gnss, x0, P0, Q, R, gate_on, gamma)
%EKF6_RUN Multi-rate EKF for the 6-state Dubins navigation model.
% Predicts at every IMU sample; updates with GNSS once its time stamp is
% reached (nearest-IMU-sample matching, see k_gnss below). No ground
% truth is used anywhere in this function.
%
% out.X        : 6xN estimated state history (a posteriori where updated)
% out.Pdiag    : 6xN diagonal of P at each IMU time (for 3-sigma plots)
% out.Pfull    : 6x6xN full covariance history (needed for NEES in Part E)
% out.d2       : 1xM Mahalanobis distance^2 for every GNSS epoch
% out.accepted : 1xM logical, true if the epoch passed the gate (or gate_on=false)
% out.K_last   : last Kalman gain applied (for inspecting outage behaviour)
% out.K_hist   : 6x2xM Kalman gain applied at each epoch (NaN if rejected)
% out.Ppos_prior : 1xM P11+P22 of the PRIOR covariance P_k|k-1 at each epoch

N = numel(t_imu);
M = numel(t_gnss);
dt = t_imu(2) - t_imu(1);
H = [1 0 0 0 0 0; 0 1 0 0 0 0];

% Map each GNSS sample to the nearest IMU index (max error = dt/2).
k_gnss = round((t_gnss - t_imu(1)) / dt) + 1;
k_gnss = min(max(k_gnss, 1), N);
gnss_at_k = zeros(1, N);
gnss_at_k(k_gnss) = 1:M;   % if two GNSS samples map to the same k (shouldn't happen here), the later one wins

X = zeros(6, N);
Pfull = zeros(6, 6, N);
Pdiag = zeros(6, N);
d2 = nan(1, M);
accepted = false(1, M);
K_last = zeros(6,2);
K_hist = nan(6, 2, M);
Ppos_prior = nan(1, M);

X(:,1) = x0;
Pfull(:,:,1) = P0;
Pdiag(:,1) = diag(P0);

% Handle a GNSS sample that maps to the very first IMU index (k=1),
% which the k=2:N predict/update loop below would otherwise never see.
j1 = gnss_at_k(1);
if j1 > 0
    z = z_gnss(:,j1);
    innov = z - H*x0;
    S = H*P0*H' + R;
    d2(j1) = innov' * (S\innov);
    Ppos_prior(j1) = P0(1,1) + P0(2,2);
    do_update = true;
    if gate_on
        do_update = d2(j1) <= gamma;
    end
    accepted(j1) = do_update;
    if do_update
        K = P0*H' / S;
        K_last = K;
        K_hist(:,:,j1) = K;
        X(:,1) = x0 + K*innov;
        X(3,1) = wrap_angle(X(3,1));
        IKH = eye(6) - K*H;
        Pfull(:,:,1) = IKH*P0*IKH' + K*R*K';
        Pdiag(:,1) = diag(Pfull(:,:,1));
    end
end

for k = 2:N
    x_prev = X(:,k-1);
    P_prev = Pfull(:,:,k-1);
    am = u_imu(1,k-1); wm = u_imu(2,k-1);

    %% Predict
    psi_prev = x_prev(3); V_prev = x_prev(4); ba = x_prev(5); bw = x_prev(6);
    x_pred = x_prev + dt * [V_prev*cos(psi_prev); V_prev*sin(psi_prev); wm-bw; am-ba; 0; 0];
    x_pred(3) = wrap_angle(x_pred(3));

    F = eye(6);
    F(1,3) = -V_prev*sin(psi_prev)*dt;  F(1,4) = cos(psi_prev)*dt;
    F(2,3) =  V_prev*cos(psi_prev)*dt;  F(2,4) = sin(psi_prev)*dt;
    F(3,6) = -dt;
    F(4,5) = -dt;

    P_pred = F*P_prev*F' + Q;

    %% Update (if a GNSS epoch maps to this IMU step)
    j = gnss_at_k(k);
    if j > 0
        z = z_gnss(:,j);
        innov = z - H*x_pred;
        S = H*P_pred*H' + R;
        d2(j) = innov' * (S\innov);
        Ppos_prior(j) = P_pred(1,1) + P_pred(2,2);

        do_update = true;
        if gate_on
            do_update = d2(j) <= gamma;
        end
        accepted(j) = do_update;

        if do_update
            K = P_pred*H' / S;
            K_last = K;
            K_hist(:,:,j) = K;
            x_upd = x_pred + K*innov;
            x_upd(3) = wrap_angle(x_upd(3));
            IKH = eye(6) - K*H;
            P_upd = IKH*P_pred*IKH' + K*R*K';   % Joseph form
        else
            x_upd = x_pred;
            P_upd = P_pred;
        end
    else
        x_upd = x_pred;
        P_upd = P_pred;
    end

    X(:,k) = x_upd;
    Pfull(:,:,k) = P_upd;
    Pdiag(:,k) = diag(P_upd);
end

out.X = X;
out.Pdiag = Pdiag;
out.Pfull = Pfull;
out.d2 = d2;
out.accepted = accepted;
out.k_gnss = k_gnss;
out.K_last = K_last;
out.K_hist = K_hist;
out.Ppos_prior = Ppos_prior;
end
