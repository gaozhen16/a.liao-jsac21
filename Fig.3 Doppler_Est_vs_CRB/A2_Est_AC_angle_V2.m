function [Est_miu_AC_set,Est_niu_AC_set,Est_Azi_AC_set,Est_Ele_AC_set] = A2_Est_AC_angle_V2(H_Beam_Alignment_AC_angle,G_ls_init,alpha_init,...
    tau_init,fs,fc,k_index_com,Power_Tx_ii,sigma_no,awgn_en,K_miu_AC_Re,K_miu_AC_Im,K_niu_AC_Re,K_niu_AC_Im,Ptx_gain_ii,...
    Delta_AC,miu_AC_set_init_error,niu_AC_set_init_error,I_AC_h_bar,I_AC_v_bar,Azi_AC_set_init_error,Ele_AC_set_init_error,Adjust_num,s_AC_angle_set)
    
    [K_sub,N_OFDM_AC_angle,L] = size(H_Beam_Alignment_AC_angle); K = K_sub*L;
    Y_AC_ang_com = zeros(N_OFDM_AC_angle,K_sub,L);
    for ll_ang_ac = 1:L
        for kk_ang_ll = 1:K_sub
            k_index_ll_ang = k_index_com(kk_ang_ll,ll_ang_ac);
            for nn_AC_ll = 1:N_OFDM_AC_angle
                Y_AC_ang_com(nn_AC_ll,kk_ang_ll,ll_ang_ac) = alpha_init(ll_ang_ac)*...
                    exp(-1i*2*pi*tau_init(ll_ang_ac)*(-1/2+(k_index_ll_ang-1)/K)*fs)*H_Beam_Alignment_AC_angle(kk_ang_ll,nn_AC_ll,ll_ang_ac)*...
                    s_AC_angle_set(kk_ang_ll,ll_ang_ac)*Ptx_gain_ii + awgn_en*sigma_no*(normrnd(0,1) + 1i*normrnd(0,1))/sqrt(2);
            end
        end
    end
    
    for Ad_num = 1:Adjust_num
        Est_miu_AC_set_temp = zeros(L,1); Est_niu_AC_set_temp = zeros(L,1);
        Est_Azi_AC_set_temp = zeros(L,1); Est_Ele_AC_set_temp = zeros(L,1);
        for ll_ang_2 = 1:L
            if Ad_num == 1
                Y_AC_ang_ll = Y_AC_ang_com(:,:,ll_ang_2);
            else
                Mtx_a_vec_BD_Adjust_ll = zeros(N_OFDM_AC_angle,K_sub);
                for kk_adust = 1:K_sub
                    Mtx_a_vec_BD_Adjust_ll(:,kk_adust) = A0_Angle_Compensation(I_AC_h_bar,I_AC_v_bar,Azi_AC_set_init_error(ll_ang_2),...
                        Ele_AC_set_init_error(ll_ang_2),Azi_AC_set_init_xx(ll_ang_2),Ele_AC_set_init_xx(ll_ang_2),k_index_com(kk_adust,ll_ang_2),K,fs,fc);
                end
                Y_AC_ang_ll = conj(Mtx_a_vec_BD_Adjust_ll).*Y_AC_ang_com(:,:,ll_ang_2);
            end
            [miu_est_AC_ll_temp, niu_est_AC_ll_temp] = TDU_ESPRIT_Algorithm(Y_AC_ang_ll, K_miu_AC_Re, K_miu_AC_Im, K_niu_AC_Re, K_niu_AC_Im);
            miu_est_AC_ll = Find_real_value(miu_est_AC_ll_temp/Delta_AC, miu_AC_set_init_error(ll_ang_2), Delta_AC);
            niu_est_AC_ll = Find_real_value(niu_est_AC_ll_temp/Delta_AC, niu_AC_set_init_error(ll_ang_2), Delta_AC);
            Est_miu_AC_set_temp(ll_ang_2) = miu_est_AC_ll;
            Est_niu_AC_set_temp(ll_ang_2) = niu_est_AC_ll;
            Est_Ele_AC_set_temp(ll_ang_2) = asin(niu_est_AC_ll/pi);
            Est_Azi_AC_set_temp(ll_ang_2) = asin(miu_est_AC_ll/(pi*cos(Est_Ele_AC_set_temp(ll_ang_2))));
        end
        Azi_AC_set_init_xx = Est_Azi_AC_set_temp; Ele_AC_set_init_xx = Est_Ele_AC_set_temp;
    end
    Est_miu_AC_set = Est_miu_AC_set_temp; Est_niu_AC_set = Est_niu_AC_set_temp;
    Est_Azi_AC_set = Est_Azi_AC_set_temp; Est_Ele_AC_set = Est_Ele_AC_set_temp;
    
end