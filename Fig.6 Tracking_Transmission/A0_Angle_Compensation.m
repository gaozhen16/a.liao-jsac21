function a_vec_angle_comp_kl = A0_Angle_Compensation(I_h,I_v,Azi_error,Ele_error,Azi_est,Ele_est,k,K,fs,fc)
% Generate array response vector of UPA for angle compensation

    i_h = (0:1:I_h-1)'; i_v = (0:1:I_v-1)';
    a_vec_h_error_k = exp(-1i*pi*((k-1)/K-0.5)*fs/fc*i_h*sin(Azi_error)*cos(Ele_error));
    a_vec_v_error_k = exp(-1i*pi*((k-1)/K-0.5)*fs/fc*i_v*sin(Ele_error));
    a_vec_angle_error_k1 = kron(a_vec_v_error_k,a_vec_h_error_k);
    
    a_vec_h_est_k = exp(1i*pi*((k-1)/K-0.5)*fs/fc*i_h*sin(Azi_est)*cos(Ele_est));
    a_vec_v_est_k = exp(1i*pi*((k-1)/K-0.5)*fs/fc*i_v*sin(Ele_est));
    a_vec_angle_est_k1 = kron(a_vec_v_est_k,a_vec_h_est_k);
    
    a_vec_angle_comp_kl = a_vec_angle_error_k1.*a_vec_angle_est_k1;
    
end