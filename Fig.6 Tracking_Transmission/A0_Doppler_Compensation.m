function a_vec_Dop_Comp_kl = A0_Doppler_Compensation(N_Dop,Doppler_z_init_error,lambda_z,T_sym,k,K,fs)
% Generate Doppler compensation matrix
    
    radial_vt_BS_set_init_eeror = Doppler_z_init_error*lambda_z;
    n_Dop = (0:1:N_Dop-1).';
    a_vec_Dop_Comp_kl = exp(-1i*2*pi*((k-1)/K-0.5)*fs/3e8*n_Dop*T_sym*radial_vt_BS_set_init_eeror(ll));
    
end