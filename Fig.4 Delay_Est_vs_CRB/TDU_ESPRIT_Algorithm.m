function [miu_est, niu_est] = TDU_ESPRIT_Algorithm(y_com_angle, K_miu_Re, K_miu_Im, K_niu_Re, K_niu_Im)
    
    % Transform complex-valued data matrix into the real-valued matrix
    [M, N] = size(y_com_angle);
    y_com_real_tmp = [y_com_angle, flipud(eye(M))*conj(y_com_angle)*flipud(eye(N))];   % Applying the forward backward averaging technique.
    y_com_real = real(Q(M)'*y_com_real_tmp*Q(2*N));
    
    % Select the first column of Left Singular Matrix for SVD of y_com_real
    [e_s, ~, ~] = svds(y_com_real, 1);	% Partial SVD
    
    % TLS
    U_miu_1 = K_miu_Re*e_s; U_miu_2 = K_miu_Im*e_s; e_miu = [U_miu_1, U_miu_2];
    U_niu_1 = K_niu_Re*e_s; U_niu_2 = K_niu_Im*e_s; e_niu = [U_niu_1, U_niu_2];
    % Diagonal elements of Sn(S1,S2) in ascending order when E_1'*E_1 is a conjugate symmetry matrix for EVD.
    [U_miu, ~] = eig(e_miu'*e_miu); e_miu_12 = U_miu(1,1); e_miu_22 = U_miu(2,1);
    [U_niu, ~] = eig(e_niu'*e_niu); e_niu_12 = U_niu(1,1); e_niu_22 = U_niu(2,1);
    phi_miu = - e_miu_12/e_miu_22;
    phi_niu = - e_niu_12/e_niu_22;
    
    miu_est = 2*atan(phi_miu);
    niu_est = 2*atan(phi_niu);
end

function Q_N = Q(N)
% Generating a unitary matrix Q_N, which defined as follows:
    if mod(N,2) == 0
        k = N/2;
        I_k = eye(k);
        J_k = flipud(eye(k));
        Q_N = [I_k,1j*I_k;...
               J_k,-1j*J_k]/sqrt(2);
    else
        k = (N-1)/2;
        I_k = eye(k);
        J_k = flipud(eye(k));
        zero_col = zeros(k,1);
        Q_N = [I_k,zero_col,1j*I_k;...
               zero_col.',sqrt(2),zero_col.';...
               J_k,zero_col,-1j*J_k]/sqrt(2);
    end
end