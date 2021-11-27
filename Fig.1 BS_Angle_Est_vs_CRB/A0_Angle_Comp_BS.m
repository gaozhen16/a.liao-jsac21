function a_vec_angle_comp_kl = A0_Angle_Comp_BS(I_BS_h,I_BS_v,Azi_BS_error,Ele_BS_error,Azi_BS_est,Ele_BS_est,k,K,fs,fc)
% Generate array response vector of UPA for angle compensation at the BS

    i_BS_h = (0:1:I_BS_h-1)'; i_BS_v = (0:1:I_BS_v-1)';
    a_vec_BS_h_error_k = exp(-1i*pi*((k-1)/K-0.5)*fs/fc*i_BS_h*sin(Azi_BS_error)*cos(Ele_BS_error));
    a_vec_BS_v_error_k = exp(-1i*pi*((k-1)/K-0.5)*fs/fc*i_BS_v*sin(Ele_BS_error));
    a_vec_angle_error_k1 = kron(a_vec_BS_v_error_k,a_vec_BS_h_error_k);
    
    a_vec_BS_h_est_k = exp(1i*pi*((k-1)/K-0.5)*fs/fc*i_BS_h*sin(Azi_BS_est)*cos(Ele_BS_est));
    a_vec_BS_v_est_k = exp(1i*pi*((k-1)/K-0.5)*fs/fc*i_BS_v*sin(Ele_BS_est));
    a_vec_angle_est_k1 = kron(a_vec_BS_v_est_k,a_vec_BS_h_est_k);
    
    a_vec_angle_comp_kl = a_vec_angle_error_k1.*a_vec_angle_est_k1;
    
end