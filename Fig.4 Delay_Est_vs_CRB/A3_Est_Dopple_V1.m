function [Est_niu_Doppler,Est_Doppler_set] = A3_Est_Dopple_V1(H_Beam_Alignment_Doppler,G_ls_init,alpha_init,Phase_Doppler_Dop_set,...
    tau_init,fs,k_index_com,N_OFDM_Doppler,Power_Tx_ii,sigma_no,awgn_en,Ptx_gain_ii,T_sym,s_dop_set)
    
    [K_sub,L] = size(k_index_com); K = K_sub*L;
    Y_dop_bar_set = zeros(N_OFDM_Doppler,K_sub,L);
    for ll_dop = 1:L
        for kk_dop_ll = 1:K_sub
            k_index_dop_ll = k_index_com(kk_dop_ll,ll_dop);
            for nn_dop_ll = 1:N_OFDM_Doppler
                if size(H_Beam_Alignment_Doppler,1) == 1
                    Y_dop_bar_set(nn_dop_ll,kk_dop_ll,ll_dop) = alpha_init(ll_dop)*exp(1i*Phase_Doppler_Dop_set(ll_dop,nn_dop_ll))*...
                        exp(-1i*2*pi*tau_init(ll_dop)*(-1/2+(k_index_dop_ll-1)/K)*fs)*H_Beam_Alignment_Doppler(ll_dop)*s_dop_set(kk_dop_ll,ll_dop)...
                        *Ptx_gain_ii + awgn_en*sigma_no*(normrnd(0,1) + 1i*normrnd(0,1))/sqrt(2);
                else
                    Y_dop_bar_set(nn_dop_ll,kk_dop_ll,ll_dop) = alpha_init(ll_dop)*exp(1i*Phase_Doppler_Dop_set(ll_dop,nn_dop_ll))*...
                        exp(-1i*2*pi*tau_init(ll_dop)*(-1/2+(k_index_dop_ll-1)/K)*fs)*H_Beam_Alignment_Doppler(kk_dop_ll,ll_dop)*...
                        s_dop_set(kk_dop_ll,ll_dop)*Ptx_gain_ii + awgn_en*sigma_no*(normrnd(0,1) + 1i*normrnd(0,1))/sqrt(2);
                end
            end
        end
    end
    Est_Doppler_set = zeros(L,1);
    for ll_dop_est = 1:L
        Est_Doppler_set(ll_dop_est) = TLS_ESPRIT_Algorithm(Y_dop_bar_set(:,:,ll_dop_est),1, T_sym);
    end
    Est_niu_Doppler = 2*pi*Est_Doppler_set*T_sym;
    
end