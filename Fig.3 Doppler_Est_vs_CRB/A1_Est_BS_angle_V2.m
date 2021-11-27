function [Est_miu_BS_set,Est_niu_BS_set,Est_Azi_BS_set,Est_Ele_BS_set] = A1_Est_BS_angle_V2(H_Beam_Alignment_BS_angle,G_ls_init,alpha_init,...
    k_index_com,fs,fc,Power_Tx_ii,sigma_no,awgn_en,K_miu_BS_Re,K_miu_BS_Im,K_niu_BS_Re,K_niu_BS_Im,Ptx_gain_ii,...
    Delta_BS,miu_BS_set_init_error,niu_BS_set_init_error,I_BS_h_bar,I_BS_v_bar,Azi_BS_set_init_error,Ele_BS_set_init_error,Adjust_num,s_BS_angle_set)

    [K_sub,N_OFDM_BS_angle,L] = size(H_Beam_Alignment_BS_angle); K = K_sub*L;
    Y_BS_ang_com = zeros(N_OFDM_BS_angle,K_sub,L);
    for ll_ang_bs = 1:L
        for kk_ang_ll = 1:K_sub
            for nn_BS_ll = 1:N_OFDM_BS_angle
                Y_BS_ang_com(nn_BS_ll,kk_ang_ll,ll_ang_bs) = alpha_init(ll_ang_bs)*...
                    H_Beam_Alignment_BS_angle(kk_ang_ll,nn_BS_ll,ll_ang_bs)*s_BS_angle_set(kk_ang_ll,ll_ang_bs)*Ptx_gain_ii...
                    + awgn_en*sigma_no*(normrnd(0,1) + 1i*normrnd(0,1))/sqrt(2);
            end
        end
    end
    
    for Ad_num = 1:Adjust_num
        Est_miu_BS_set_temp = zeros(L,1); Est_niu_BS_set_temp = zeros(L,1);
        Est_Azi_BS_set_temp = zeros(L,1); Est_Ele_BS_set_temp = zeros(L,1);
        for ll_ang_2 = 1:L
            if Ad_num == 1
                Y_BS_ang_ll = Y_BS_ang_com(:,:,ll_ang_2);
            else
                Mtx_a_vec_BD_Adjust_ll = zeros(N_OFDM_BS_angle,K_sub);
                for kk_adust = 1:K_sub
                    Mtx_a_vec_BD_Adjust_ll(:,kk_adust) = A0_Angle_Compensation(I_BS_h_bar,I_BS_v_bar,Azi_BS_set_init_error(ll_ang_2),...
                        Ele_BS_set_init_error(ll_ang_2),Azi_BS_set_init_xx(ll_ang_2),Ele_BS_set_init_xx(ll_ang_2),k_index_com(kk_adust,ll_ang_2),K,fs,fc);
                end
                Y_BS_ang_ll = conj(Mtx_a_vec_BD_Adjust_ll).*Y_BS_ang_com(:,:,ll_ang_2);
            end
            [miu_est_BS_ll_temp, niu_est_BS_ll_temp] = TDU_ESPRIT_Algorithm(Y_BS_ang_ll, K_miu_BS_Re, K_miu_BS_Im, K_niu_BS_Re, K_niu_BS_Im);
            miu_est_BS_ll = Find_real_value(miu_est_BS_ll_temp/Delta_BS, miu_BS_set_init_error(ll_ang_2), Delta_BS);
            niu_est_BS_ll = Find_real_value(niu_est_BS_ll_temp/Delta_BS, niu_BS_set_init_error(ll_ang_2), Delta_BS);
            Est_miu_BS_set_temp(ll_ang_2) = miu_est_BS_ll;
            Est_niu_BS_set_temp(ll_ang_2) = niu_est_BS_ll;
            Est_Ele_BS_set_temp(ll_ang_2) = asin(niu_est_BS_ll/pi);
            Est_Azi_BS_set_temp(ll_ang_2) = asin(miu_est_BS_ll/(pi*cos(Est_Ele_BS_set_temp(ll_ang_2))));
        end
        Azi_BS_set_init_xx = Est_Azi_BS_set_temp; Ele_BS_set_init_xx = Est_Ele_BS_set_temp;
    end
    Est_miu_BS_set = Est_miu_BS_set_temp; Est_niu_BS_set = Est_niu_BS_set_temp;
    Est_Azi_BS_set = Est_Azi_BS_set_temp; Est_Ele_BS_set = Est_Ele_BS_set_temp;
    
end