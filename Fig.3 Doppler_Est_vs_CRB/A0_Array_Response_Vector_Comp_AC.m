function a_vec_AC_Comp_k = A0_Array_Response_Vector_Comp_AC(ii_AC,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp_ll,N_AC_h,N_AC_v,...
    M_AC_h,M_AC_v,N_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp,M_AC_Comp_h,Azi_AC_est,Ele_AC_est,k,K,fs,fc)
% Generate array response vector of UPA for compensation
    
    n_AC_h = (0:1:N_AC_h-1)'; n_AC_v = (0:1:N_AC_v-1)';
    if N_AC_Comp_h == 1 && N_AC_Comp_v == 1
        a_vec_AC_h_k = exp(-1i*pi*((k-1)/K-0.5)*fs/fc*n_AC_h*sin(Azi_AC_est)*cos(Ele_AC_est)); % /sqrt(N_AC_h)
        a_vec_v_AC_k = exp(-1i*pi*((k-1)/K-0.5)*fs/fc*n_AC_v*sin(Ele_AC_est)); % /sqrt(N_AC_v)
        a_vec_AC_Comp_temp = kron(a_vec_v_AC_k,a_vec_AC_h_k);
        a_vec_AC_Comp_k = zeros(size(a_vec_AC_Comp_temp));
        a_vec_AC_Comp_k(Index_W_RF_Sub_Array(:,ii_AC)) = a_vec_AC_Comp_temp(Index_W_RF_Sub_Array(:,ii_AC));
    else
        a_vec_AC_h_k = exp(-1i*pi*((k-1)/K-0.5)*fs/fc*n_AC_h*sin(Azi_AC_est)*cos(Ele_AC_est)); % /sqrt(N_AC_h)
        a_vec_v_AC_k = exp(-1i*pi*((k-1)/K-0.5)*fs/fc*n_AC_v*sin(Ele_AC_est)); % /sqrt(N_AC_v)
        a_vec_Comp_k_temp = kron(a_vec_v_AC_k,a_vec_AC_h_k);
        i_AC_v = ceil(ii_AC/I_AC_h);
        i_AC_h = ii_AC-(i_AC_v-1)*I_AC_h;
        A_Mtx_k_sub = zeros(M_AC_h,M_AC_v);
        for m_AC_comp = 1:M_AC_Comp
            m_AC_Comp_v = ceil(m_AC_comp/M_AC_Comp_h);
            m_AC_Comp_h = m_AC_comp - (m_AC_Comp_v-1)*M_AC_Comp_h;
            A_Mtx_k_sub((m_AC_Comp_h-1)*N_AC_Comp_h+1:m_AC_Comp_h*N_AC_Comp_h,(m_AC_Comp_v-1)*N_AC_Comp_v+1:m_AC_Comp_v*N_AC_Comp_v) = ...
                a_vec_Comp_k_temp(Index_AC_Ant_Comp_ll(round(N_AC_Comp/2),m_AC_comp))*ones(N_AC_Comp_h,N_AC_Comp_v);
        end
        A_Mtx_k = zeros(N_AC_h,N_AC_v);
        A_Mtx_k((i_AC_h-1)*M_AC_h+1:i_AC_h*M_AC_h,(i_AC_v-1)*M_AC_v+1:i_AC_v*M_AC_v) = A_Mtx_k_sub;
        a_vec_AC_Comp_k = reshape(A_Mtx_k,N_AC_h*N_AC_v,1);
    end
    
end