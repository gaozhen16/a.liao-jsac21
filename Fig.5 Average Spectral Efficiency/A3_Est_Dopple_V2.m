function [Est_niu_Doppler,Est_Doppler_set,Y_dop_bar_set] = A3_Est_Dopple_V2(H_Beam_Alignment_Doppler,G_ls_init,alpha_init,Doppler_init_k_set,...
    tau_init,fs,k_index_com,N_OFDM_Doppler,Power_Tx_ii,sigma_no,awgn_en,Ptx_gain_ii,T_sym,s_dop_set,Adjust_num,Doppler_z_init_error,lambda_z)
    
    [K_sub,L] = size(k_index_com); K = K_sub*L;
    Y_dop_bar_set = zeros(N_OFDM_Doppler,K_sub,L);
    for ll_dop = 1:L
        for kk_dop_ll = 1:K_sub
            k_index_dop_ll = k_index_com(kk_dop_ll,ll_dop);
            for nn_dop_ll = 1:N_OFDM_Doppler
                if size(H_Beam_Alignment_Doppler,1) == 1
                    Y_dop_bar_set(nn_dop_ll,kk_dop_ll,ll_dop) = alpha_init(ll_dop)*exp(1i*2*pi*(nn_dop_ll-1)*T_sym*Doppler_init_k_set(kk_dop_ll,ll_dop))*...
                        exp(-1i*2*pi*tau_init(ll_dop)*(-1/2+(k_index_dop_ll-1)/K)*fs)*H_Beam_Alignment_Doppler(ll_dop)*s_dop_set(kk_dop_ll,ll_dop)...
                        *Ptx_gain_ii + awgn_en*sigma_no*(normrnd(0,1) + 1i*normrnd(0,1))/sqrt(2);
                else
                    Y_dop_bar_set(nn_dop_ll,kk_dop_ll,ll_dop) = alpha_init(ll_dop)*exp(1i*2*pi*(nn_dop_ll-1)*T_sym*Doppler_init_k_set(kk_dop_ll,ll_dop))*...
                        exp(-1i*2*pi*tau_init(ll_dop)*(-1/2+(k_index_dop_ll-1)/K)*fs)*H_Beam_Alignment_Doppler(kk_dop_ll,ll_dop)*...
                        s_dop_set(kk_dop_ll,ll_dop)*Ptx_gain_ii + awgn_en*sigma_no*(normrnd(0,1) + 1i*normrnd(0,1))/sqrt(2);
                end
            end
        end
    end
    
    n_Dop = (1:N_OFDM_Doppler).'; radial_vt_BS_set_init_error = Doppler_z_init_error*lambda_z;
    Est_Doppler_set = zeros(L,Adjust_num);
    for ad_num = 1:Adjust_num
        if  ad_num == 1
            for ll_dop_est1 = 1:L
                Est_Doppler_set(ll_dop_est1,ad_num) = TLS_ESPRIT_Algorithm(Y_dop_bar_set(:,:,ll_dop_est1),1, T_sym);
            end
        elseif ad_num == 2
            for ll_dop_est2 = 1:L
                Mtx_a_vec_Dop_Adjust_ll_2 = zeros(N_OFDM_Doppler,K_sub);
                for kk_dop2 = 1:K_sub
                    Mtx_a_vec_Dop_Adjust_ll_2(:,kk_dop2) = exp(1i*2*pi*(n_Dop-1)*T_sym*radial_vt_BS_set_init_error(ll_dop_est2)*...
                        ((k_index_com(kk_dop2,ll_dop_est2)-1)/K-0.5)*fs/3e8);
                end
                Y_dop_bar_ll_2 = conj(Mtx_a_vec_Dop_Adjust_ll_2).*Y_dop_bar_set(:,:,ll_dop_est2);
                Est_Doppler_set(ll_dop_est2,ad_num) = TLS_ESPRIT_Algorithm(Y_dop_bar_ll_2,1, T_sym);
            end
            radial_vt_BS_set_init_est = Est_Doppler_set(:,ad_num)*lambda_z;
        else
            for ll_dop_est3 = 1:L
                Mtx_a_vec_Dop_Adjust_ll_3 = zeros(N_OFDM_Doppler,K_sub);
                for kk_dop3 = 1:K_sub
                    Mtx_a_vec_Dop_Adjust_ll_3(:,kk_dop3) = exp(1i*2*pi*(n_Dop-1)*T_sym*radial_vt_BS_set_init_est(ll_dop_est3)*...
                        ((k_index_com(kk_dop3,ll_dop_est3)-1)/K-0.5)*fs/3e8);
                end
                Y_dop_bar_ll_3 = conj(Mtx_a_vec_Dop_Adjust_ll_3).*Y_dop_bar_set(:,:,ll_dop_est3);
                Est_Doppler_set(ll_dop_est3,ad_num) = TLS_ESPRIT_Algorithm(Y_dop_bar_ll_3,1, T_sym);
            end
        end
    end
    Est_niu_Doppler = 2*pi*Est_Doppler_set*T_sym;
    
end