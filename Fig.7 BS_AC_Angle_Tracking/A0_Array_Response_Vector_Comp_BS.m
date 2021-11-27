function a_vec_BS_Comp_k = A0_Array_Response_Vector_Comp_BS(Index_BS_Ant_Comp,N_BS_h,N_BS_v,N_BS_Comp,N_BS_Comp_h,N_BS_Comp_v,M_BS_Comp,...
    M_BS_Comp_h,Azi_BS_est,Ele_BS_est,k,K,fs,fc)
% Generate array response vector of UPA for compensation
    
    n_BS_h = (0:1:N_BS_h-1)'; n_BS_v = (0:1:N_BS_v-1)';
    if N_BS_Comp_h == 1 && N_BS_Comp_v == 1
        a_vec_BS_h_k = exp(-1i*pi*((k-1)/K-0.5)*fs/fc*n_BS_h*sin(Azi_BS_est)*cos(Ele_BS_est)); % /sqrt(N_BS_h)
        a_vec_BS_v_k = exp(-1i*pi*((k-1)/K-0.5)*fs/fc*n_BS_v*sin(Ele_BS_est)); % /sqrt(N_BS_v)
        a_vec_BS_Comp_k = kron(a_vec_BS_v_k,a_vec_BS_h_k);
    else
        a_vec_BS_h_k = exp(-1i*pi*((k-1)/K-0.5)*fs/fc*n_BS_h*sin(Azi_BS_est)*cos(Ele_BS_est)); % /sqrt(N_BS_h)
        a_vec_BS_v_k = exp(-1i*pi*((k-1)/K-0.5)*fs/fc*n_BS_v*sin(Ele_BS_est)); % /sqrt(N_BS_v)
        a_vec_Comp_k_temp = kron(a_vec_BS_v_k,a_vec_BS_h_k);
        A_Mtx_BS_k_temp = zeros(N_BS_h,N_BS_v);
        for m_BS_comp = 1:M_BS_Comp
            m_BS_Comp_v = ceil(m_BS_comp/M_BS_Comp_h);
            m_BS_Comp_h = m_BS_comp - (m_BS_Comp_v-1)*M_BS_Comp_h;
            A_Mtx_BS_k_temp((m_BS_Comp_h-1)*N_BS_Comp_h+1:m_BS_Comp_h*N_BS_Comp_h,(m_BS_Comp_v-1)*N_BS_Comp_v+1:m_BS_Comp_v*N_BS_Comp_v) = ...
                a_vec_Comp_k_temp(Index_BS_Ant_Comp(round(N_BS_Comp/2),m_BS_comp))*ones(N_BS_Comp_h,N_BS_Comp_v);
        end
        a_vec_BS_Comp_k = reshape(A_Mtx_BS_k_temp,N_BS_h*N_BS_v,1);
    end
    
end