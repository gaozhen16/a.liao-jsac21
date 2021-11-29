clc,clear

%% System parameters
fc = 100e9;
lambda_z = 3e8/fc;
d_ant = lambda_z/2;
BW = 1e9;
fs = BW;
Ts = 1/fs;
vt = 200;
FFT_len = 2048;
K = FFT_len;
CIR_max_L = 128;
No_sym = K+CIR_max_L;
T_sym = No_sym*Ts;
N_bits = 0;
sigma_2_alpha = 1;
awgn_en = 1;
L = 2;
N_BS_h = 200; N_BS_v = N_BS_h;
N_BS = N_BS_h*N_BS_v;
Ns_BS = 1;
M_AC_h = 200; M_AC_v = M_AC_h;
M_AC = M_AC_h*M_AC_v;
I_AC_h = 1; N_AC_h = I_AC_h*M_AC_h;
I_AC_v = L; N_AC_v = I_AC_v*M_AC_v;
I_AC = I_AC_h*I_AC_v;
N_AC = N_AC_h*N_AC_v;
I_BS_h_bar = 5;
I_BS_v_bar = I_BS_h_bar;
I_BS_bar = I_BS_h_bar*I_BS_v_bar;
Delta_BS = 1;
N_BS_h_bar = N_BS_h - Delta_BS*(I_BS_h_bar - 1);
N_BS_v_bar = N_BS_v - Delta_BS*(I_BS_v_bar - 1);
N_BS_bar = N_BS_h_bar*N_BS_v_bar;
Index_W_RF_BS_Angle = A0_Gene_Index_W_RF_1BSangle(N_BS_bar,I_BS_bar,N_BS_v_bar,N_BS_h_bar,I_BS_h_bar,N_BS_v,N_BS_h,Delta_BS);
I_AC_h_bar = 5;
I_AC_v_bar = I_AC_h_bar;
I_AC_bar = I_AC_h_bar*I_AC_v_bar;
Delta_AC = 1;
M_AC_h_bar = M_AC_h - Delta_AC*(I_AC_h_bar - 1);
M_AC_v_bar = M_AC_v - Delta_AC*(I_AC_v_bar - 1);
M_AC_bar = M_AC_h_bar*M_AC_v_bar;
Index_W_RF_AC_Angle = A0_Gene_Index_W_RF_1ACangle(M_AC_bar,I_AC_bar,I_AC,I_AC_h,I_AC_h_bar,M_AC_h_bar,M_AC_v_bar,M_AC_h,M_AC_v,Delta_AC,N_AC_h,N_AC_v);
Index_W_RF_Sub_Array = A0_Gene_Index_W_RF_2Sub_Array(N_AC_h,N_AC_v,M_AC,M_AC_h,M_AC_v,I_AC,I_AC_h);
K_index_set = 1:K; K_sub = K/L;
D_subc = 2;
if mod(K,L) == 0
	k_index_com = reshape(K_index_set,L,K_sub).';
else
    error('Error: K should be the integral number of L !');
end
[K_miu_BS_Re,K_miu_BS_Im,K_niu_BS_Re,K_niu_BS_Im] = A0_Gene_2D_ESPRIT_Para(I_BS_h_bar,I_BS_v_bar);
[K_miu_AC_Re,K_miu_AC_Im,K_niu_AC_Re,K_niu_AC_Im] = A0_Gene_2D_ESPRIT_Para(I_AC_h_bar,I_AC_v_bar);
N_BS_Comp_h = 5; N_BS_Comp_v = N_BS_Comp_h;
N_BS_Comp = N_BS_Comp_h*N_BS_Comp_v;
M_BS_Comp_h = N_BS_h/N_BS_Comp_h; M_BS_Comp_v = N_BS_v/N_BS_Comp_v;
M_BS_Comp = M_BS_Comp_h*M_BS_Comp_v;
Index_BS_Ant_Comp = A0_Gene_Index_Antenna_Comp(N_BS_Comp,M_BS_Comp,N_BS_Comp_h,N_BS_Comp_v,M_BS_Comp_h,N_BS_h,N_BS_v);
N_AC_Comp_h = 5; N_AC_Comp_v = N_AC_Comp_h;
N_AC_Comp = N_AC_Comp_h*N_AC_Comp_v;
M_AC_Comp_h = M_AC_h/N_AC_Comp_h; M_AC_Comp_v = M_AC_v/N_AC_Comp_v;
M_AC_Comp = M_AC_Comp_h*M_AC_Comp_v;
Index_AC_Ant_Comp = A0_Gene_Index_Antenna_Comp(N_AC_Comp,M_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp_h,N_AC_h,N_AC_v,I_AC,I_AC_h,M_AC_h,M_AC_v);

%% Preliminary
load('Z0_CoorEF_vt_Azi_Ele_BS_AC_tau_set.mat');
load('Z0_alpha_init_set.mat');
load('Z0_Other_Parameters_CRB_v0.mat');
D = 10e3;
E = 100e3; F = 100e3;
G = 200e3;
Radius_AC = 50e3;
G_AC = 1; G_BS = 1;
E_temp = E_temp_set(rand_num);
F_temp = F_temp_set(rand_num);
Coo_BS1_A = [0,0,D].';
Coo_BS2_B = [0,G,D].';
Coo_Air_C_init = [E_temp,F_temp,0].';
Distance_init(1) = norm(Coo_BS1_A-Coo_Air_C_init);
Distance_init(2) = norm(Coo_BS2_B-Coo_Air_C_init);
G_ls_init = G_AC*G_BS*lambda_z^2./(4*pi*Distance_init).^2;
unit_d = unit_d_set(:,rand_num);
vt_vec = vt*unit_d;
radial_vt_BS_set_init = Calculate_Radial_Velocity(vt_vec, Coo_Air_C_init, Coo_BS1_A, Coo_BS2_B, L);
Doppler_z_init = radial_vt_BS_set_init./lambda_z;
Doppler_init_k_set = zeros(K_sub,L);
for ll_dop_k = 1:L
    Doppler_init_k_set(:,ll_dop_k) = Doppler_z_init(ll_dop_k) + radial_vt_BS_set_init(ll_dop_k)*((k_index_com(:,ll_dop_k)-1)/K-0.5)*fs/3e8;
end
Azi_BS_set_init = Azi_BS_set_init_set(:,rand_num);
Ele_BS_set_init = Ele_BS_set_init_set(:,rand_num);
miu_BS_set_init = pi*sin(Azi_BS_set_init).*cos(Ele_BS_set_init);
niu_BS_set_init = pi*sin(Ele_BS_set_init);
Azi_AC_set_init = Azi_AC_set_init_set(:,rand_num);
Ele_AC_set_init = Ele_AC_set_init_set(:,rand_num);
miu_AC_set_init = pi*sin(Azi_AC_set_init).*cos(Ele_AC_set_init);
niu_AC_set_init = pi*sin(Ele_AC_set_init);
tau_init = tau_init_set(:,rand_num);
miu_tau_set_init = -2*pi*fs*tau_init/K;
alpha_init = alpha_init_set(:,rand_num);
N_OFDM_BS_angle = I_BS_bar;
N_OFDM_AC_angle = I_AC_bar;
N_OFDM_Angle_Est = N_OFDM_BS_angle + N_OFDM_AC_angle;
N_OFDM_Doppler = 6;
N_OFDM_Delay = 10;
N_OFDM = N_OFDM_Angle_Est + N_OFDM_Doppler + N_OFDM_Delay;
T_CCT = N_OFDM*T_sym;
niu_Doppler_init = 2*pi*Doppler_z_init*T_sym;
m_CCT_init = 1;
n_ofdm_doppler = 1:N_OFDM_Doppler;
Phase_Doppler_Dop_set = 2*pi*Doppler_z_init*T_sym*((n_ofdm_doppler-1)+(m_CCT_init-1)*N_OFDM);
Angle_diff = deg2rad(5);
miu_BS_set_init_error = pi*sin(Azi_BS_set_init_error).*cos(Ele_BS_set_init_error);
niu_BS_set_init_error = pi*sin(Ele_BS_set_init_error);
miu_AC_set_init_error = pi*sin(Azi_AC_set_init_error).*cos(Ele_AC_set_init_error);
niu_AC_set_init_error = pi*sin(Ele_AC_set_init_error);
Jacobian_Mat_Azi_Ele_BS = Gene_Jacobian_Mat_Azi_Ele(miu_BS_set_init, niu_BS_set_init, Delta_BS, lambda_z, d_ant, L);
Jacobian_Mat_Azi_Ele_AC = Gene_Jacobian_Mat_Azi_Ele(miu_AC_set_init, niu_AC_set_init, Delta_AC, lambda_z, d_ant, L);
Jacobian_Mat_Doppler = Gene_Jacobian_Mat_Doppler(niu_Doppler_init, T_sym, L);
Jacobian_Mat_Delay = Gene_Jacobian_Mat_tau(miu_tau_set_init, K, 1, D_subc, L);

%% set snr
iterMax = 2e3;
Power_Tx_dBm = 0:5:150;
Ptx_gain = ones(1,length(Power_Tx_dBm));
SNR_dBs = -120:10:20;
sigma2_NSD = -174;
sigma2_no = 1e-3*10^(sigma2_NSD/10)*BW;
sigma_no = sqrt(sigma2_no);
snr_0_set = [-20,20];
MSE_Delay_set_temp_0 = zeros(length(SNR_dBs),length(snr_0_set));
CRB_Delay_set = zeros(length(SNR_dBs),length(snr_0_set));
MSE_Delay_set_temp_1 = zeros(length(SNR_dBs),length(snr_0_set));

%% 
Power_Tx_ii = (1e-3*10^(Power_Tx_dBm(1)/10)); Ptx_gain_ii = Ptx_gain(1);
Est_Azi_BS_set_02 = zeros(L,length(snr_0_set)); Est_Ele_BS_set_02 = zeros(L,length(snr_0_set));
for snr_bs_0 = 1:length(snr_0_set)
    if snr_bs_0 == 1
        a_vec_AC_init_set = zeros(N_AC,L); a_vec_BS_init_set = zeros(N_BS,L);
    end
    H_Beam_Alignment_BS_angle_0 = zeros(N_OFDM_BS_angle,L);
    for bs_ll_0 = 1:L
        a_vec_AC_init_ll = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init(bs_ll_0),Ele_AC_set_init(bs_ll_0));
        a_vec_BS_init_ll = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init(bs_ll_0),Ele_BS_set_init(bs_ll_0));
        if snr_bs_0 == 1
            a_vec_AC_init_set(:,bs_ll_0) = a_vec_AC_init_ll; a_vec_BS_init_set(:,bs_ll_0) = a_vec_BS_init_ll;
        end
        a_vec_AC_init_ll_error = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init_error(bs_ll_0),Ele_AC_set_init_error(bs_ll_0));
        a_vec_BS_init_ll_error = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init_error(bs_ll_0),Ele_BS_set_init_error(bs_ll_0));
        f_RF_AC_init_ll = zeros(N_AC,1);
        f_RF_AC_init_ll(Index_W_RF_Sub_Array(:,bs_ll_0)) = Quantize(a_vec_AC_init_ll_error(Index_W_RF_Sub_Array(:,bs_ll_0)),N_bits)/sqrt(M_AC); % 
        Equ_a_AC_ll_init = a_vec_AC_init_ll'*f_RF_AC_init_ll;
        Subarray_BS_ll_init = Quantize(a_vec_BS_init_ll_error(Index_W_RF_BS_Angle(:,1)),N_bits)/sqrt(N_BS_bar); % 
        for nn_bs = 1:N_OFDM_BS_angle
            w_RF_BS_init_ll_nn = zeros(N_BS,1);
            w_RF_BS_init_ll_nn(Index_W_RF_BS_Angle(:,nn_bs)) = Subarray_BS_ll_init;
            Equ_a_BS_ll_init = w_RF_BS_init_ll_nn'*a_vec_BS_init_ll;
            H_Beam_Alignment_BS_angle_0(nn_bs,bs_ll_0) = Equ_a_BS_ll_init*Equ_a_AC_ll_init;
        end
    end
    sigma_no_0 = sqrt(10^(-(snr_0_set(snr_bs_0)/10)));
    [~,~,Est_Azi_BS_set_02(:,snr_bs_0),Est_Ele_BS_set_02(:,snr_bs_0)] = A1_Est_BS_angle_V1(H_Beam_Alignment_BS_angle_0,G_ls_init,alpha_init,...
        k_index_com,Power_Tx_ii,sigma_no_0,awgn_en,K_miu_BS_Re, K_miu_BS_Im, K_niu_BS_Re, K_niu_BS_Im, ...
        Ptx_gain_ii,Delta_BS,miu_BS_set_init_error,niu_BS_set_init_error,s_BS_angle_set);
end

Est_Azi_AC_set_02 = zeros(L,length(snr_0_set)); Est_Ele_AC_set_02 = zeros(L,length(snr_0_set));
f_RF_BS_init_set_0 = zeros(N_BS,L,length(snr_0_set));
for snr_ac_0 = 1:length(snr_0_set)
    H_Beam_Alignment_AC_angle_0 = zeros(N_OFDM_AC_angle,L);
    for ac_ll_0 = 1:L
        a_vec_BS_init_ll_est = A0_Array_Response_Vector(N_BS_h,N_BS_v,Est_Azi_BS_set_02(ac_ll_0,snr_ac_0),Est_Ele_BS_set_02(ac_ll_0,snr_ac_0));
        a_vec_AC_init_ll_error = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init_error(ac_ll_0),Ele_AC_set_init_error(ac_ll_0));
        f_RF_BS_init_ll = Quantize(a_vec_BS_init_ll_est,N_bits)/sqrt(N_BS);
        f_RF_BS_init_set_0(:,ac_ll_0,snr_ac_0) = f_RF_BS_init_ll;
        Equ_a_BS_ll_init = a_vec_BS_init_set(:,ac_ll_0)'*f_RF_BS_init_ll;
        Subarray_AC_ll_init = Quantize(a_vec_AC_init_ll_error(Index_W_RF_AC_Angle(:,1,ac_ll_0)),N_bits)/sqrt(M_AC_bar); % 
        for nn_ac_0 = 1:N_OFDM_AC_angle
            w_RF_AC_init_ll_nn = zeros(N_AC,1);
            w_RF_AC_init_ll_nn(Index_W_RF_AC_Angle(:,nn_ac_0,ac_ll_0)) = Subarray_AC_ll_init;
            Equ_a_AC_ll_init = w_RF_AC_init_ll_nn'*a_vec_AC_init_set(:,ac_ll_0);
            H_Beam_Alignment_AC_angle_0(nn_ac_0,ac_ll_0) = Equ_a_AC_ll_init*Equ_a_BS_ll_init;
        end
    end
    sigma_no_0 = sqrt(10^(-(snr_0_set(snr_ac_0)/10)));
    [~,~,Est_Azi_AC_set_02(:,snr_ac_0),Est_Ele_AC_set_02(:,snr_ac_0)] = A2_Est_AC_angle_V1(H_Beam_Alignment_AC_angle_0,G_ls_init,alpha_init,...
        tau_init,fs,k_index_com,Power_Tx_ii,sigma_no_0,awgn_en,K_miu_AC_Re,K_miu_AC_Im,K_niu_AC_Re,K_niu_AC_Im,...
        Ptx_gain_ii,Delta_AC,miu_AC_set_init_error,niu_AC_set_init_error,s_AC_angle_set);
end

Est_Doppler_set_02 = zeros(L,length(snr_0_set)); w_RF_AC_init_set_0 = zeros(N_AC,L,length(snr_0_set));
H_Beam_Alignment_Doppler_set_0 = zeros(1,L,length(snr_0_set));
for snr_dop_0 = 1:length(snr_0_set)
    for ll_dop_0 = 1:L
        a_vec_AC_init_ll_est = A0_Array_Response_Vector(N_AC_h,N_AC_v,Est_Azi_AC_set_02(ll_dop_0,snr_dop_0),Est_Ele_AC_set_02(ll_dop_0,snr_dop_0));
        w_RF_AC_init_ll = zeros(N_AC,1);
        w_RF_AC_init_ll(Index_W_RF_Sub_Array(:,ll_dop_0)) = Quantize(a_vec_AC_init_ll_est(Index_W_RF_Sub_Array(:,ll_dop_0)),N_bits)/sqrt(M_AC); % 
        w_RF_AC_init_set_0(:,ll_dop_0,snr_dop_0) = w_RF_AC_init_ll;
        Equ_a_BS_ll_init_dop = a_vec_BS_init_set(:,ll_dop_0)'*f_RF_BS_init_set_0(:,ll_dop_0,snr_dop_0);
        Equ_a_AC_ll_init_dop = w_RF_AC_init_ll'*a_vec_AC_init_set(:,ll_dop_0);
        H_Beam_Alignment_Doppler_set_0(1,ll_dop_0,snr_dop_0) = Equ_a_AC_ll_init_dop*Equ_a_BS_ll_init_dop;
    end
    sigma_no_0 = sqrt(10^(-(snr_0_set(snr_dop_0)/10)));
    [~,Est_Doppler_set_02(:,snr_dop_0)] = A3_Est_Dopple_V1(H_Beam_Alignment_Doppler_set_0(:,:,snr_dop_0),G_ls_init,alpha_init,Phase_Doppler_Dop_set,...
        tau_init,fs,k_index_com,N_OFDM_Doppler,Power_Tx_ii,sigma_no_0,awgn_en,Ptx_gain_ii,T_sym,s_dop_set);
end

%%
H_Beam_Alignment_Delay_set_1 = zeros(K_sub,L,length(snr_0_set));
for snr_del_1 = 1:length(snr_0_set)
    for ll_del1 = 1:L
        for kk_del1 = 1:K_sub
            a_vec_BS_init_lk = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init(ll_del1),Ele_BS_set_init(ll_del1),...
                k_index_com(kk_del1,ll_del1),K,fs,fc);
            a_vec_AC_init_lk = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init(ll_del1),Ele_AC_set_init(ll_del1),...
                k_index_com(kk_del1,ll_del1),K,fs,fc);
            a_vec_BS_Comp_lk = A0_Array_Response_Vector_Comp_BS(Index_BS_Ant_Comp,N_BS_h,N_BS_v,N_BS_Comp,N_BS_Comp_h,N_BS_Comp_v,M_BS_Comp,M_BS_Comp_h,...
                Est_Azi_BS_set_02(ll_del1,snr_del_1),Est_Ele_BS_set_02(ll_del1,snr_del_1),k_index_com(kk_del1,ll_del1),K,fs,fc); % 
            a_vec_AC_Comp_lk = A0_Array_Response_Vector_Comp_AC(ll_del1,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp(:,:,ll_del1),...
                N_AC_h,N_AC_v,M_AC_h,M_AC_v,N_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp,M_AC_Comp_h,...
                Est_Azi_AC_set_02(ll_del1,snr_del_1),Est_Ele_AC_set_02(ll_del1,snr_del_1),k_index_com(kk_del1,ll_del1),K,fs,fc); % 
            Equ_a_BS_ll_init_del = (a_vec_BS_init_lk.*a_vec_BS_Comp_lk)'*f_RF_BS_init_set_0(:,ll_del1,snr_del_1); % 
            Equ_a_AC_ll_init_del = w_RF_AC_init_set_0(:,ll_del1,snr_del_1)'*(a_vec_AC_init_lk.*a_vec_AC_Comp_lk); % 
            H_Beam_Alignment_Delay_set_1(kk_del1,ll_del1,snr_del_1) = Equ_a_AC_ll_init_del*Equ_a_BS_ll_init_del;
        end
    end
end

radial_vt_BS_set_est = lambda_z*Est_Doppler_set_02;
n_ofdm_delay = 1:N_OFDM_Delay;
Phase_Delay_Dop_Comp_set_0 = zeros(L,N_OFDM_Delay,length(snr_0_set));
Phase_Delay_Dop_Comp_set_1 = zeros(L,N_OFDM_Delay,K_sub,length(snr_0_set));
for snr_ddc = 1:length(snr_0_set)
    Phase_Delay_Dop_Comp_set_0(:,:,snr_ddc) = 2*pi*(Doppler_z_init-Est_Doppler_set_02(:,snr_ddc))*(n_ofdm_delay-1)*T_sym;
    Doppler_est_k_set = zeros(K_sub,L);
    for ll_ddc = 1:L
        Doppler_est_k_set(:,ll_ddc) = Est_Doppler_set_02(ll_ddc,snr_ddc) + radial_vt_BS_set_est(ll_ddc,snr_ddc)*((k_index_com(:,ll_ddc)-1)/K-0.5)*fs/3e8;
    end
    for kk_ddc = 1:K_sub
        Phase_Delay_Dop_Comp_set_1(:,:,kk_ddc,snr_ddc) = 2*pi*(Doppler_init_k_set(kk_ddc,:)-Doppler_est_k_set(kk_ddc,:)).'*(n_ofdm_delay-1)*T_sym;
    end
end

%% Start
for ii = 1:length(SNR_dBs)
    Power_Tx_ii = (1e-3*10^(Power_Tx_dBm(ii)/10));
    Ptx_gain_ii = Ptx_gain(ii);
    sigma_no = sqrt(10^(-(SNR_dBs(ii)/10)));
    tic
    for iter = 1:iterMax
        for snr_del = 1:length(snr_0_set)
            [Est_miu_tau_0,Est_tau_set_0,Y_tau_bar_set_0] = A4_Est_Delay_V1(H_Beam_Alignment_Doppler_set_0(:,:,snr_del),G_ls_init,alpha_init,...
                Phase_Delay_Dop_Comp_set_0(:,:,snr_del),tau_init,fs,k_index_com,N_OFDM_Delay,Power_Tx_ii,sigma_no,awgn_en,D_subc,Ptx_gain_ii,s_tau_set);
            MSE_Delay_set_temp_0(ii,snr_del) = MSE_Delay_set_temp_0(ii,snr_del) + norm((Est_tau_set_0 - tau_init)*fs)^2/L;
            
            if ii == 1 && iter == 1
                Mat_ZPZ_Del = zeros(1,1,L);
                for ll_crb_del1 = 1:L
                    Mat_Z_Del_ll = Diff_Para_Delay(miu_tau_set_init(ll_crb_del1), k_index_com(:,ll_crb_del1));
                    a_bar_Del = exp(1i*(k_index_com(:,ll_crb_del1)-1)*miu_tau_set_init(ll_crb_del1));
                    Mat_P_A_orth_Del = eye(size(a_bar_Del,1)) - a_bar_Del*inv(a_bar_Del'*a_bar_Del)*a_bar_Del';
                    Mat_ZPZ_Del(:,:,ll_crb_del1) = Mat_Z_Del_ll'*Mat_P_A_orth_Del*Mat_Z_Del_ll;
                end
            end
            if iter == 1
                CRB_miu_tau_set_temp = zeros(L,1); CRB_Delay_set_temp = zeros(L,1);
                for ll_crb_del2 = 1:L
                    gamma_del_ll = alpha_init(ll_crb_del2)*H_Beam_Alignment_Doppler_set_0(:,ll_crb_del2,snr_del)*exp(1i*pi*fs*tau_init(ll_crb_del2));
                    Tx_eff_Del_set_ll = exp(1i*Phase_Delay_Dop_Comp_set_0(ll_crb_del2,:,snr_del)).'.*s_tau_set(:,ll_crb_del2); % gamma_del_ll*
                    Mat_1D_dop_ll = 0;
                    for nn_Del = 1:N_OFDM_Delay
                        Mat_1D_dop_ll = Mat_1D_dop_ll + Tx_eff_Del_set_ll(nn_Del)'*Mat_ZPZ_Del(:,:,ll_crb_del2)*Tx_eff_Del_set_ll(nn_Del);
                    end
                    CRB_miu_tau_set_temp(ll_crb_del2) = sigma_no^2*inv(real(Mat_1D_dop_ll))/2/abs(gamma_del_ll/N_OFDM_Delay)^2; % 
                    CRB_Delay_set_temp(ll_crb_del2) = Jacobian_Mat_Delay(:,:,ll_crb_del2)*CRB_miu_tau_set_temp(ll_crb_del2)*Jacobian_Mat_Delay(:,:,ll_crb_del2).';
                end
                CRB_Delay_set(ii,snr_del) = sum(CRB_Delay_set_temp)/L;
            end
            
            [Est_miu_tau_1,Est_tau_set_1,Y_tau_bar_set_1] = A4_Est_Delay_V1(H_Beam_Alignment_Delay_set_1(:,:,snr_del),G_ls_init,alpha_init,...
                Phase_Delay_Dop_Comp_set_1(:,:,:,snr_del),tau_init,fs,k_index_com,N_OFDM_Delay,Power_Tx_ii,sigma_no,awgn_en,D_subc,Ptx_gain_ii,s_tau_set);
            MSE_Delay_set_temp_1(ii,snr_del) = MSE_Delay_set_temp_1(ii,snr_del) + norm((Est_tau_set_1 - tau_init)*fs)^2/L;
            
        end
        
        %%
        if mod(iter,100) == 0
            toc
            disp(['    SNR_dB = ' num2str(SNR_dBs(ii)) ',   iter = ' num2str(iter)])
        end
    end
    disp(['  Finish SNR_dB = ' num2str(SNR_dBs(ii))])
end
MSE_Delay_set_0 = MSE_Delay_set_temp_0/iterMax;
MSE_Delay_set_1 = MSE_Delay_set_temp_1/iterMax;
disp('Finish  All')

MarkerSize = 3;
LineWidth = 1.2;
Fontsize = 15;

RCRB_Delay_set = sqrt(CRB_Delay_set);
RMSE_Delay_set_0 = sqrt(MSE_Delay_set_0);
RMSE_Delay_set_1 = sqrt(MSE_Delay_set_1);

figure
semilogy(SNR_dBs,RMSE_Delay_set_1(:,1),'--bp','LineWidth',LineWidth); hold on; grid on;
semilogy(SNR_dBs,RMSE_Delay_set_0(:,1),'--ro','LineWidth',LineWidth);
semilogy(SNR_dBs,RCRB_Delay_set(:,1),'--k','LineWidth',LineWidth+0.5);
semilogy(SNR_dBs,RMSE_Delay_set_1(:,2),'-bp','LineWidth',LineWidth);
semilogy(SNR_dBs,RMSE_Delay_set_0(:,2),'-ro','LineWidth',LineWidth);
semilogy(SNR_dBs,RCRB_Delay_set(:,2),'-k','LineWidth',LineWidth+0.5);
xlabel('SNR [dB]','Fontsize',Fontsize),ylabel('RMSE','Fontsize',Fontsize);
title('Delay vs CRB','Fontsize',Fontsize);
set(gca, 'GridLineStyle', '-.','FontSize',Fontsize, 'linewidth',1.5,'Fontname','Times New Roman');
set(gcf, 'position', [700 300 650 550]); axis normal;
h1 = legend('Triple squint','No triple squint','CRLB, SNR = -20 dB',...
    'Triple squint','No triple squint','CRLB, SNR = 20 dB','Location','southwest');
set(h1,'Fontsize',11);
