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
K_index_set = (1:K).'; K_sub = K/L;
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
niu_Doppler_init = 2*pi*Doppler_z_init*T_sym;
Doppler_diff = Doppler_z_init*1e-2;
Doppler_z_init_error = Doppler_z_init+Doppler_diff;
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
Power_Tx_dBm = 0:5:150;
Ptx_gain = ones(1,length(Power_Tx_dBm));
Power_Tx_ii = (1e-3*10^(Power_Tx_dBm(1)/10));
Ptx_gain_ii = Ptx_gain(1);
load('Z0_Other_Est_BS_AC_angles_CRB_v1.mat');
Delta_BS_tr = 1;
N_BS_h_bar_tr = N_BS_h - Delta_BS_tr*(I_BS_h_bar - 1);
N_BS_v_bar_tr = N_BS_v - Delta_BS_tr*(I_BS_v_bar - 1);
N_BS_bar_tr = N_BS_h_bar_tr*N_BS_v_bar_tr;
Index_W_RF_BS_Angle_tr = A0_Gene_Index_W_RF_1BSangle(N_BS_bar_tr,I_BS_bar,N_BS_v_bar_tr,N_BS_h_bar_tr,I_BS_h_bar,N_BS_v,N_BS_h,Delta_BS_tr);
Delta_AC_tr = 4;
M_AC_h_bar_tr = M_AC_h - Delta_AC_tr*(I_AC_h_bar - 1);
M_AC_v_bar_tr = M_AC_v - Delta_AC_tr*(I_AC_v_bar - 1);
M_AC_bar_tr = M_AC_h_bar_tr*M_AC_v_bar_tr;
Index_W_RF_AC_Angle_tr = A0_Gene_Index_W_RF_1ACangle(M_AC_bar_tr,I_AC_bar,I_AC,I_AC_h,I_AC_h_bar,M_AC_h_bar_tr,M_AC_v_bar_tr,M_AC_h,M_AC_v,Delta_AC_tr,N_AC_h,N_AC_v);

%%
N_OFDM_C = 70;
T_CCT = N_OFDM_C*T_sym;
N_CCT = 100;
T_CCT_tot = N_CCT*T_CCT;
% rho_theta_phi_BS = randsrc*pi/12;
% rho_theta_phi_AC = randsrc*pi/4;
% rho_Doppler = -0.01*Doppler_z_init;
% rho_tau = -tau_init/2;
% rho_alpha = randsrc*alpha_init/2;
% save Z0_Other_Parameters_Change_Rate.mat rho_theta_phi_BS rho_theta_phi_AC rho_Doppler rho_tau rho_alpha
load('Z0_Other_Parameters_Change_Rate.mat');
Azi_BS_set_trac = Azi_BS_set_init + rho_theta_phi_BS*T_CCT_tot;
Ele_BS_set_trac = Ele_BS_set_init + rho_theta_phi_BS*T_CCT_tot;
miu_BS_set_trac = pi*sin(Azi_BS_set_trac).*cos(Ele_BS_set_trac);
niu_BS_set_trac = pi*sin(Ele_BS_set_trac);
Azi_AC_set_trac = Azi_AC_set_init + rho_theta_phi_AC*T_CCT_tot;
Ele_AC_set_trac = Ele_AC_set_init + rho_theta_phi_AC*T_CCT_tot;
miu_AC_set_trac = pi*sin(Azi_AC_set_trac).*cos(Ele_AC_set_trac);
niu_AC_set_trac = pi*sin(Ele_AC_set_trac);
Doppler_z_trac = Doppler_z_init + rho_Doppler*T_CCT_tot;
radial_vt_BS_set_trac = lambda_z*Doppler_z_trac;
tau_trac = tau_init + rho_tau*T_CCT_tot;
alpha_trac = alpha_init + rho_alpha*T_CCT_tot;
Jacobian_Mat_Azi_Ele_BS_tr = Gene_Jacobian_Mat_Azi_Ele(miu_BS_set_trac, niu_BS_set_trac, Delta_BS_tr, lambda_z, d_ant, L);
Jacobian_Mat_Azi_Ele_AC_tr = Gene_Jacobian_Mat_Azi_Ele(miu_AC_set_trac, niu_AC_set_trac, Delta_AC_tr, lambda_z, d_ant, L);

H_Beam_Alignment_BS_angle_trac = zeros(N_OFDM_BS_angle,L);
a_vec_AC_trac_set = zeros(N_AC,L); a_vec_BS_trac_set = zeros(N_BS,L);
for ll_bs_tr1 = 1:L
    a_vec_AC_trac_ll = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_trac(ll_bs_tr1),Ele_AC_set_trac(ll_bs_tr1));
    a_vec_BS_trac_ll = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_trac(ll_bs_tr1),Ele_BS_set_trac(ll_bs_tr1));
    a_vec_AC_trac_set(:,ll_bs_tr1) = a_vec_AC_trac_ll; a_vec_BS_trac_set(:,ll_bs_tr1) = a_vec_BS_trac_ll;
    a_vec_AC_trac_ll_est = A0_Array_Response_Vector(N_AC_h,N_AC_v,Est_Azi_AC_set(ll_bs_tr1),Est_Ele_AC_set(ll_bs_tr1));
    a_vec_BS_trac_ll_est = A0_Array_Response_Vector(N_BS_h,N_BS_v,Est_Azi_BS_set(ll_bs_tr1),Est_Ele_BS_set(ll_bs_tr1));
    f_RF_AC_trac_ll = zeros(N_AC,1);
    f_RF_AC_trac_ll(Index_W_RF_Sub_Array(:,ll_bs_tr1)) = Quantize(a_vec_AC_trac_ll_est(Index_W_RF_Sub_Array(:,ll_bs_tr1)),N_bits)/sqrt(M_AC);
    Equ_a_AC_ll_trac = a_vec_AC_trac_ll'*f_RF_AC_trac_ll;
    Subarray_BS_ll_trac = Quantize(a_vec_BS_trac_ll_est(Index_W_RF_BS_Angle_tr(:,1)),N_bits)/sqrt(N_BS_bar_tr);
    for nn_bs_tr1 = 1:N_OFDM_BS_angle
        w_RF_BS_trac_ll_nn = zeros(N_BS,1);
        w_RF_BS_trac_ll_nn(Index_W_RF_BS_Angle_tr(:,nn_bs_tr1)) = Subarray_BS_ll_trac;
        Equ_a_BS_ll_trac = w_RF_BS_trac_ll_nn'*a_vec_BS_trac_ll;
        H_Beam_Alignment_BS_angle_trac(nn_bs_tr1,ll_bs_tr1) = Equ_a_BS_ll_trac*Equ_a_AC_ll_trac;
    end
end
% sigma_no_trac_bs = sqrt(10^(-(-60/10)));
% [Est_miu_BS_set_trac,Est_niu_BS_set_trac,Est_Azi_BS_set_trac,Est_Ele_BS_set_trac] = A1_Est_BS_angle_V1(H_Beam_Alignment_BS_angle_trac,G_ls_init,alpha_trac,...
%     k_index_com,Power_Tx_ii,sigma_no_trac_bs,awgn_en,K_miu_BS_Re, K_miu_BS_Im, K_niu_BS_Re, K_niu_BS_Im, ...
%     Ptx_gain_ii,Delta_BS_tr,Est_miu_BS_set,Est_niu_BS_set,s_BS_angle_set);
% save Z0_Other_Est_Azi_Ele_BS_set_trac.mat Est_Azi_BS_set_trac Est_Ele_BS_set_trac
load('Z0_Other_Est_Azi_Ele_BS_set_trac.mat');

H_Beam_Alignment_AC_angle_trac = zeros(N_OFDM_AC_angle,L);
f_RF_BS_trac_AC_set = zeros(N_BS,L); w_RF_AC_trac_AC_1_set = zeros(N_AC,L);
for ll_ac_tr1 = 1:L
    a_vec_BS_trac_ll_est = A0_Array_Response_Vector(N_BS_h,N_BS_v,Est_Azi_BS_set_trac(ll_ac_tr1),Est_Ele_BS_set_trac(ll_ac_tr1));
    a_vec_AC_trac_ll_est = A0_Array_Response_Vector(N_AC_h,N_AC_v,Est_Azi_AC_set(ll_ac_tr1),Est_Ele_AC_set(ll_ac_tr1));
    f_RF_BS_trac_ll = Quantize(a_vec_BS_trac_ll_est,N_bits)/sqrt(N_BS);
    f_RF_BS_trac_AC_set(:,ll_ac_tr1) = f_RF_BS_trac_ll;
    Equ_a_BS_ll_trac = a_vec_BS_trac_set(:,ll_ac_tr1)'*f_RF_BS_trac_ll;
    Subarray_AC_ll_trac = Quantize(a_vec_AC_trac_ll_est(Index_W_RF_AC_Angle_tr(:,1,ll_ac_tr1)),N_bits)/sqrt(M_AC_bar_tr);
    for nn_ac_tr1 = 1:N_OFDM_AC_angle
        w_RF_AC_trac_ll_nn = zeros(N_AC,1);
        w_RF_AC_trac_ll_nn(Index_W_RF_AC_Angle_tr(:,nn_ac_tr1,ll_ac_tr1)) = Subarray_AC_ll_trac;
        if nn_ac_tr1 == 1
            w_RF_AC_trac_AC_1_set(:,ll_ac_tr1) = w_RF_AC_trac_ll_nn;
        end
        Equ_a_AC_ll_trac = w_RF_AC_trac_ll_nn'*a_vec_AC_trac_set(:,ll_ac_tr1);
        H_Beam_Alignment_AC_angle_trac(nn_ac_tr1,ll_ac_tr1) = Equ_a_AC_ll_trac*Equ_a_BS_ll_trac;
    end
end

H_Beam_Alignment_AC_angle_trac_BeSq = zeros(K_sub,N_OFDM_AC_angle,L);
for ac_ll_tr2 = 1:L
    a_vec_BS_trac_ll_est = A0_Array_Response_Vector(N_BS_h,N_BS_v,Est_Azi_BS_set_trac(ac_ll_tr2),Est_Ele_BS_set_trac(ac_ll_tr2));
    a_vec_AC_trac_ll_est = A0_Array_Response_Vector(N_AC_h,N_AC_v,Est_Azi_AC_set(ac_ll_tr2),Est_Ele_AC_set(ac_ll_tr2));
    f_RF_BS_trac_ll = Quantize(a_vec_BS_trac_ll_est,N_bits)/sqrt(N_BS);
    Subarray_AC_ll_trac = Quantize(a_vec_AC_trac_ll_est(Index_W_RF_AC_Angle_tr(:,1,ac_ll_tr2)),N_bits)/sqrt(M_AC_bar_tr);
    Index_AC_Ant_Comp_ll = Index_AC_Ant_Comp(:,:,ac_ll_tr2);
    for kk_ac_tr_ll = 1:K_sub
        a_vec_BS_trac_lk = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_trac(ac_ll_tr2),Ele_BS_set_trac(ac_ll_tr2),...
            k_index_com(kk_ac_tr_ll,ac_ll_tr2),K,fs,fc);
        a_vec_BS_Comp_lk = A0_Array_Response_Vector_Comp_BS(Index_BS_Ant_Comp,N_BS_h,N_BS_v,N_BS_Comp,N_BS_Comp_h,N_BS_Comp_v,M_BS_Comp,M_BS_Comp_h,...
            Est_Azi_BS_set_trac(ac_ll_tr2),Est_Ele_BS_set_trac(ac_ll_tr2),k_index_com(kk_ac_tr_ll,ac_ll_tr2),K,fs,fc);
        a_vec_AC_trac_lk = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_trac(ac_ll_tr2),Ele_AC_set_trac(ac_ll_tr2),...
            k_index_com(kk_ac_tr_ll,ac_ll_tr2),K,fs,fc);
        a_vec_AC_Comp_lk = A0_Array_Response_Vector_Comp_AC(ac_ll_tr2,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp_ll,N_AC_h,N_AC_v,M_AC_h,M_AC_v,...
            N_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp,M_AC_Comp_h,Est_Azi_AC_set(ac_ll_tr2),Est_Ele_AC_set(ac_ll_tr2),...
            k_index_com(kk_ac_tr_ll,ac_ll_tr2),K,fs,fc);
        Equ_a_BS_lk_trac = (a_vec_BS_trac_lk.*a_vec_BS_Comp_lk)'*f_RF_BS_trac_ll;
        for nn_ac_tr = 1:N_OFDM_AC_angle
            w_RF_AC_trac_ll_nn = zeros(N_AC,1);
            w_RF_AC_trac_ll_nn(Index_W_RF_AC_Angle_tr(:,nn_ac_tr,ac_ll_tr2)) = Subarray_AC_ll_trac;
            Equ_a_AC_lk_trac = w_RF_AC_trac_ll_nn'*(a_vec_AC_trac_lk.*a_vec_AC_Comp_lk);
            H_Beam_Alignment_AC_angle_trac_BeSq(kk_ac_tr_ll,nn_ac_tr,ac_ll_tr2) = Equ_a_AC_lk_trac*Equ_a_BS_lk_trac;
        end
    end
end

%% set snr
iterMax = 2e3;
SNR_dBs = -120:10:-30;
MSE_Azi_AC_set_tr_temp_0 = zeros(length(SNR_dBs),1);
MSE_Ele_AC_set_tr_temp_0 = zeros(length(SNR_dBs),1);
CRB_Azi_AC_set_trac = zeros(length(SNR_dBs),1);
CRB_Ele_AC_set_trac = zeros(length(SNR_dBs),1);
MSE_Azi_AC_set_tr_temp_1 = zeros(length(SNR_dBs),1);
MSE_Ele_AC_set_tr_temp_1 = zeros(length(SNR_dBs),1);

%% Start
for ii = 1:length(SNR_dBs)
    sigma_no_trac = sqrt(10^(-(SNR_dBs(ii)/10)));
    tic
    for iter = 1:iterMax
        [Est_miu_AC_set_trac_0,Est_niu_AC_set_trac_0,Est_Azi_AC_set_trac_0,Est_Ele_AC_set_trac_0] = A2_Est_AC_angle_V1(H_Beam_Alignment_AC_angle_trac,...
            G_ls_init,alpha_trac,tau_trac,fs,k_index_com,Power_Tx_ii,sigma_no_trac,awgn_en,K_miu_AC_Re,K_miu_AC_Im,...
            K_niu_AC_Re,K_niu_AC_Im,Ptx_gain_ii,Delta_AC_tr,Est_miu_AC_set,Est_niu_AC_set,s_AC_angle_set);
        MSE_Azi_AC_set_tr_temp_0(ii,1) = MSE_Azi_AC_set_tr_temp_0(ii,1) + norm(Est_Azi_AC_set_trac_0 - Azi_AC_set_trac)^2/L;
        MSE_Ele_AC_set_tr_temp_0(ii,1) = MSE_Ele_AC_set_tr_temp_0(ii,1) + norm(Est_Ele_AC_set_trac_0 - Ele_AC_set_trac)^2/L;
        
        if ii == 1 && iter == 1
            Mat_ZPZ_2D_AC = zeros(2,2,L);
            for ll_crb_ac1 = 1:L
                i_AC_h_bar = (0:I_BS_h_bar-1).'; a_AC_miu_eff_ll = exp(1i*i_AC_h_bar*miu_AC_set_trac(ll_crb_ac1));
                i_AC_v_bar = (0:I_BS_v_bar-1).'; a_AC_niu_eff_ll = exp(1i*i_AC_v_bar*niu_AC_set_trac(ll_crb_ac1));
                Diff_a_AC_miu = Diff_Para_Azi_Ele(I_AC_h_bar, miu_AC_set_trac(ll_crb_ac1));
                Com_Diff_a_AC_miu = Khatri_Rao(a_AC_niu_eff_ll,Diff_a_AC_miu);
                Diff_a_AC_niu = Diff_Para_Azi_Ele(I_AC_v_bar, niu_AC_set_trac(ll_crb_ac1));
                Com_Diff_a_AC_niu = Khatri_Rao(Diff_a_AC_niu,a_AC_miu_eff_ll);
                Mat_Z_2D_AC_ll = [Com_Diff_a_AC_miu, Com_Diff_a_AC_niu];
                a_AC_bar_2D = Khatri_Rao(a_AC_niu_eff_ll,a_AC_miu_eff_ll);
                Mat_P_A_orth_2D_AC_ll = eye(size(a_AC_bar_2D,1)) - a_AC_bar_2D*inv(a_AC_bar_2D'*a_AC_bar_2D)*a_AC_bar_2D';
                Mat_ZPZ_2D_AC(:,:,ll_crb_ac1) = Mat_Z_2D_AC_ll'*Mat_P_A_orth_2D_AC_ll*Mat_Z_2D_AC_ll;
            end
        end
        if iter == 1
            CRB_Azi_AC_temp = zeros(L,1); CRB_Ele_AC_temp = zeros(L,1);
            for ll_crb_ac2 = 1:L
                Equ_a_BS_ac_crb_ll = a_vec_BS_trac_set(:,ll_crb_ac2)'*f_RF_BS_trac_AC_set(:,ll_crb_ac2);
                Equ_a_AC_ac_crb_ll = w_RF_AC_trac_AC_1_set(:,ll_crb_ac2)'*a_vec_AC_trac_set(:,ll_crb_ac2);
                gamma_angle_ac_ll = alpha_trac(ll_crb_ac2)*Equ_a_AC_ac_crb_ll*Equ_a_BS_ac_crb_ll;
                Tx_eff_AC_set_ll = exp(-1i*2*pi*tau_trac(ll_crb_ac2)*(-1/2+(k_index_com(:,ll_crb_ac2)-1)/K)*fs).*s_AC_angle_set(:,ll_crb_ac2);
                Mat_2D_AC = zeros(2,2);
                for kk_ac = 1:K_sub
                    Mat_B_2D_AC_nn = kron(eye(2),diag(Tx_eff_AC_set_ll(kk_ac)));
                    Mat_2D_AC = Mat_2D_AC + Mat_B_2D_AC_nn'*Mat_ZPZ_2D_AC(:,:,ll_crb_ac2)*Mat_B_2D_AC_nn;
                end
                CRB_miu_niu_Mtx_AC_ll = sigma_no_trac^2*inv(real(Mat_2D_AC))/2/abs(gamma_angle_ac_ll)^2;
                CRB_Azi_Ele_Mtx_AC_ll = Jacobian_Mat_Azi_Ele_AC_tr(:,:,ll_crb_ac2)*CRB_miu_niu_Mtx_AC_ll*Jacobian_Mat_Azi_Ele_AC_tr(:,:,ll_crb_ac2).';
                CRB_Azi_AC_temp(ll_crb_ac2) = CRB_Azi_Ele_Mtx_AC_ll(2,2);
                CRB_Ele_AC_temp(ll_crb_ac2) = CRB_Azi_Ele_Mtx_AC_ll(1,1);
            end
            CRB_Azi_AC_set_trac(ii) = sum(CRB_Azi_AC_temp)/L;
            CRB_Ele_AC_set_trac(ii) = sum(CRB_Ele_AC_temp)/L;
        end
        Adjust_num = 2;
        [Est_miu_AC_set_trac_1,Est_niu_AC_set_trac_1,Est_Azi_AC_set_trac_1,Est_Ele_AC_set_trac_1] = A2_Est_AC_angle_V2(H_Beam_Alignment_AC_angle_trac_BeSq,...
            G_ls_init,alpha_trac,tau_trac,fs,fc,k_index_com,Power_Tx_ii,sigma_no_trac,awgn_en,K_miu_AC_Re,K_miu_AC_Im,...
            K_niu_AC_Re,K_niu_AC_Im,Ptx_gain_ii,Delta_AC_tr,Est_miu_AC_set,Est_niu_AC_set,I_AC_h_bar,I_AC_v_bar,Est_Azi_AC_set,Est_Ele_AC_set,...
            Adjust_num,s_AC_angle_set);
        MSE_Azi_AC_set_tr_temp_1(ii,1) = MSE_Azi_AC_set_tr_temp_1(ii,1) + norm(Est_Azi_AC_set_trac_1 - Azi_AC_set_trac)^2/L;
        MSE_Ele_AC_set_tr_temp_1(ii,1) = MSE_Ele_AC_set_tr_temp_1(ii,1) + norm(Est_Ele_AC_set_trac_1 - Ele_AC_set_trac)^2/L;
        
        if mod(iter,50) == 0
            toc
            disp(['    SNR_dB = ' num2str(SNR_dBs(ii)) ',   iter = ' num2str(iter)])
        end
    end
    disp(['  Finish SNR_dB = ' num2str(SNR_dBs(ii))])
end
MSE_Azi_AC_set_trac_0 = MSE_Azi_AC_set_tr_temp_0/iterMax;
MSE_Ele_AC_set_trac_0 = MSE_Ele_AC_set_tr_temp_0/iterMax;
MSE_Azi_AC_set_trac_1 = MSE_Azi_AC_set_tr_temp_1/iterMax;
MSE_Ele_AC_set_trac_1 = MSE_Ele_AC_set_tr_temp_1/iterMax;
disp('Finish  All')

MarkerSize = 3;
LineWidth = 1.2;
Fontsize = 15;
RMSE_Azi_AC_set_trac_0 = sqrt(MSE_Azi_AC_set_trac_0); RCRB_Azi_AC_set_trac = sqrt(CRB_Azi_AC_set_trac);
RMSE_Ele_AC_set_trac_0 = sqrt(MSE_Ele_AC_set_trac_0); RCRB_Ele_AC_set_trac = sqrt(CRB_Ele_AC_set_trac);
RMSE_Azi_AC_set_trac_1 = sqrt(MSE_Azi_AC_set_trac_1);
RMSE_Ele_AC_set_trac_1 = sqrt(MSE_Ele_AC_set_trac_1);

figure
semilogy(SNR_dBs,RMSE_Azi_AC_set_trac_1,'--bp','LineWidth',LineWidth); hold on; grid on;
semilogy(SNR_dBs,RMSE_Azi_AC_set_trac_0,'--ro','LineWidth',LineWidth);
semilogy(SNR_dBs,RCRB_Azi_AC_set_trac,'--k','LineWidth',LineWidth+0.5);
semilogy(SNR_dBs,RMSE_Ele_AC_set_trac_1,'-bp','LineWidth',LineWidth);
semilogy(SNR_dBs,RMSE_Ele_AC_set_trac_0,'-ro','LineWidth',LineWidth);
semilogy(SNR_dBs,RCRB_Ele_AC_set_trac,'-k','LineWidth',LineWidth+0.5);
xlabel('SNR [dB]','Fontsize',Fontsize),ylabel('RMSE','Fontsize',Fontsize);
title('Track AC Azi & Ele vs CRB {\it{\Omega}} = 4','Fontsize',Fontsize);
set(gca, 'GridLineStyle', '-.','FontSize',Fontsize, 'linewidth',1.5,'Fontname','Times New Roman');
set(gcf, 'position', [700 300 700 500]); axis normal;
h1 = legend('{\it{¦È}}^{AC} + BeSq','{\it{¦È}}^{AC} + no BeSq','CRB {\it{¦È}}^{AC}',...
    '{\it{¦Õ}}^{AC} + BeSq','{\it{¦Õ}}^{AC} + no BeSq','CRB {\it{¦Õ}}^{AC}','Location','northeast');
set(h1,'Fontsize',11);
