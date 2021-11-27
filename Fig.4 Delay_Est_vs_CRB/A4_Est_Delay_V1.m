function [Est_miu_tau,Est_tau_set,Y_tau_bar_set] = A4_Est_Delay_V1(H_Beam_Alignment_Delay_Doppler,G_ls_init,alpha_init,Phase_Delay_Dop_Comp_set,...
    tau_init,fs,k_index_com,N_OFDM_Delay,Power_Tx_ii,sigma_no,awgn_en,D_subc,Ptx_gain_ii,s_tau_set)
    
    [K_sub,L] = size(k_index_com); K = K_sub*L;
    Y_tau_bar_set = zeros(K_sub,N_OFDM_Delay,L);
    for ll_tau = 1:L
        for nn_tau_ll = 1:N_OFDM_Delay
            for kk_tau_ll = 1:K_sub
                k_index_tau_ll = k_index_com(kk_tau_ll,ll_tau);
                if size(H_Beam_Alignment_Delay_Doppler,1) == 1
                    Y_tau_bar_set(kk_tau_ll,nn_tau_ll,ll_tau) = alpha_init(ll_tau)*exp(1i*Phase_Delay_Dop_Comp_set(ll_tau,nn_tau_ll))*...
                        exp(-1i*2*pi*tau_init(ll_tau)*(-1/2+(k_index_tau_ll-1)/K)*fs)*H_Beam_Alignment_Delay_Doppler(ll_tau)*...
                        s_tau_set(nn_tau_ll,ll_tau)*Ptx_gain_ii + awgn_en*sigma_no*(normrnd(0,1) + 1i*normrnd(0,1))/sqrt(2);
                else
                    Y_tau_bar_set(kk_tau_ll,nn_tau_ll,ll_tau) = alpha_init(ll_tau)*exp(1i*Phase_Delay_Dop_Comp_set(ll_tau,nn_tau_ll,kk_tau_ll))*...
                        exp(-1i*2*pi*tau_init(ll_tau)*(-1/2+(k_index_tau_ll-1)/K)*fs)*H_Beam_Alignment_Delay_Doppler(kk_tau_ll,ll_tau)*...
                        s_tau_set(nn_tau_ll,ll_tau)*Ptx_gain_ii + awgn_en*sigma_no*(normrnd(0,1) + 1i*normrnd(0,1))/sqrt(2);
                end
            end
        end
    end
    
    Est_tau_set = zeros(L,1);
    for ll_tau_est = 1:L
        Est_tau_set(ll_tau_est) = TLS_ESPRIT_Algorithm(Y_tau_bar_set(:,:,ll_tau_est),1, D_subc, fs, K);
    end
    Est_miu_tau = -2*pi*fs*Est_tau_set/K;
    
end