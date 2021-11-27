function [K_miu_Re,K_miu_Im,K_niu_Re,K_niu_Im] = A0_Gene_2D_ESPRIT_Para(M_h_sub,M_v_sub)
% Generate the ESPRIT Parameters
    
    % % Set 2D selection matrices
    J_M_miu_2 = Select_M2(M_h_sub);
    J_M_niu_2 = Select_M2(M_v_sub);

    J_miu_2 = kron(eye(M_v_sub),J_M_miu_2);
    J_niu_2 = kron(J_M_niu_2,eye(M_h_sub));

    % % N_¦Ì corresponding pairs of transformed selection matrices
    K_miu = Q(size(J_miu_2,1))'*J_miu_2*Q(size(J_miu_2,2));
    K_niu = Q(size(J_niu_2,1))'*J_niu_2*Q(size(J_niu_2,2));

    K_miu_Re = 2*real(K_miu); K_miu_Im = 2*imag(K_miu);
    K_niu_Re = 2*real(K_niu); K_niu_Im = 2*imag(K_niu);
    
end

function J_Mr2 = Select_M2(Mr)
% Generating selection matrices J_(n)2
J_Mr2 = [zeros(Mr-1,1),eye(Mr-1)];
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