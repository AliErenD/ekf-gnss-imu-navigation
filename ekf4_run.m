function out = ekf4_run(t_imu, u_imu, t_gnss, z_gnss, x0, P0, Q, R, ba_fix, bw_fix, gate_on, gamma)
%EKF4_RUN Multi-rate EKF for the 4-state model x=[pn,pe,psi,V].
% The IMU bias is NOT estimated: am,wm are corrected with the fixed
% static estimates (ba_fix, bw_fix) from Part A before being used.

N = numel(t_imu);
M = numel(t_gnss);
dt = t_imu(2) - t_imu(1);
H = [1 0 0 0; 0 1 0 0];

k_gnss = round((t_gnss - t_imu(1)) / dt) + 1;
k_gnss = min(max(k_gnss, 1), N);
gnss_at_k = zeros(1, N);
gnss_at_k(k_gnss) = 1:M;

X = zeros(4, N);
Pfull = zeros(4, 4, N);
Pdiag = zeros(4, N);
d2 = nan(1, M);
accepted = false(1, M);
consec_rejected = zeros(1, M); % running count of consecutive rejections (for the lockout check)

X(:,1) = x0;
Pfull(:,:,1) = P0;
Pdiag(:,1) = diag(P0);
run_reject = 0;

j1 = gnss_at_k(1);
if j1 > 0
    z = z_gnss(:,j1);
    innov = z - H*x0;
    S = H*P0*H' + R;
    d2(j1) = innov' * (S\innov);
    do_update = true;
    if gate_on
        do_update = d2(j1) <= gamma;
    end
    accepted(j1) = do_update;
    if ~do_update, run_reject = run_reject + 1; consec_rejected(j1) = run_reject; end
    if do_update
        K = P0*H' / S;
        X(:,1) = x0 + K*innov;
        X(3,1) = wrap_angle(X(3,1));
        IKH = eye(4) - K*H;
        Pfull(:,:,1) = IKH*P0*IKH' + K*R*K';
        Pdiag(:,1) = diag(Pfull(:,:,1));
    end
end

for k = 2:N
    x_prev = X(:,k-1);
    P_prev = Pfull(:,:,k-1);
    am = u_imu(1,k-1) - ba_fix;
    wm = u_imu(2,k-1) - bw_fix;

    psi_prev = x_prev(3); V_prev = x_prev(4);
    x_pred = x_prev + dt * [V_prev*cos(psi_prev); V_prev*sin(psi_prev); wm; am];
    x_pred(3) = wrap_angle(x_pred(3));

    F = eye(4);
    F(1,3) = -V_prev*sin(psi_prev)*dt;  F(1,4) = cos(psi_prev)*dt;
    F(2,3) =  V_prev*cos(psi_prev)*dt;  F(2,4) = sin(psi_prev)*dt;

    P_pred = F*P_prev*F' + Q;

    j = gnss_at_k(k);
    if j > 0
        z = z_gnss(:,j);
        innov = z - H*x_pred;
        S = H*P_pred*H' + R;
        d2(j) = innov' * (S\innov);

        do_update = true;
        if gate_on
            do_update = d2(j) <= gamma;
        end
        accepted(j) = do_update;
        if do_update, run_reject = 0; else, run_reject = run_reject + 1; end
        consec_rejected(j) = run_reject;

        if do_update
            K = P_pred*H' / S;
            x_upd = x_pred + K*innov;
            x_upd(3) = wrap_angle(x_upd(3));
            IKH = eye(4) - K*H;
            P_upd = IKH*P_pred*IKH' + K*R*K';
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
out.consec_rejected = consec_rejected;
out.k_gnss = k_gnss;
end
