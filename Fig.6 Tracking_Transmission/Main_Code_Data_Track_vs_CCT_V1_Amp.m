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
SNR_dB = -20;
sigma_no_est = sqrt(10^(-(SNR_dB/10)));

%% 
a_vec_AC_init_be_sq_set = zeros(N_AC,K_sub,L); a_vec_BS_init_be_sq_set = zeros(N_BS,K_sub,L);
H_Beam_Alignment_BS_angle_sq = zeros(K_sub,N_OFDM_BS_angle,L);
for bs_ll_1 = 1:L
    a_vec_AC_init_ll_error = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init_error(bs_ll_1),Ele_AC_set_init_error(bs_ll_1));
    a_vec_BS_init_ll_error = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init_error(bs_ll_1),Ele_BS_set_init_error(bs_ll_1));
    f_RF_AC_init_ll = zeros(N_AC,1);
    f_RF_AC_init_ll(Index_W_RF_Sub_Array(:,bs_ll_1)) = Quantize(a_vec_AC_init_ll_error(Index_W_RF_Sub_Array(:,bs_ll_1)),N_bits)/sqrt(M_AC);
    Subarray_BS_ll_init = Quantize(a_vec_BS_init_ll_error(Index_W_RF_BS_Angle(:,1)),N_bits)/sqrt(N_BS_bar);
    Index_AC_Ant_Comp_ll = Index_AC_Ant_Comp(:,:,bs_ll_1);
    for kk_bs_ll = 1:K_sub
        a_vec_AC_init_lk = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init(bs_ll_1),Ele_AC_set_init(bs_ll_1),...
            k_index_com(kk_bs_ll,bs_ll_1),K,fs,fc);
        a_vec_AC_Comp_lk = A0_Array_Response_Vector_Comp_AC(bs_ll_1,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp_ll,N_AC_h,N_AC_v,M_AC_h,M_AC_v,...
            N_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp,M_AC_Comp_h,Azi_AC_set_init_error(bs_ll_1),Ele_AC_set_init_error(bs_ll_1),...
            k_index_com(kk_bs_ll,bs_ll_1),K,fs,fc);
        a_vec_BS_init_lk = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init(bs_ll_1),Ele_BS_set_init(bs_ll_1),...
            k_index_com(kk_bs_ll,bs_ll_1),K,fs,fc);
        a_vec_BS_Comp_lk = A0_Array_Response_Vector_Comp_BS(Index_BS_Ant_Comp,N_BS_h,N_BS_v,N_BS_Comp,N_BS_Comp_h,N_BS_Comp_v,M_BS_Comp,M_BS_Comp_h,...
            Azi_BS_set_init_error(bs_ll_1),Ele_BS_set_init_error(bs_ll_1),k_index_com(kk_bs_ll,bs_ll_1),K,fs,fc);
        a_vec_AC_init_be_sq_set(:,kk_bs_ll,bs_ll_1) = a_vec_AC_init_lk;
        a_vec_BS_init_be_sq_set(:,kk_bs_ll,bs_ll_1) = a_vec_BS_init_lk;
        Equ_a_AC_lk_init = (a_vec_AC_init_lk.*a_vec_AC_Comp_lk)'*f_RF_AC_init_ll;
        for nn_bs = 1:N_OFDM_BS_angle
            w_RF_BS_init_ll_nn = zeros(N_BS,1);
            w_RF_BS_init_ll_nn(Index_W_RF_BS_Angle(:,nn_bs)) = Subarray_BS_ll_init;
            Equ_a_BS_lk_init = w_RF_BS_init_ll_nn'*(a_vec_BS_init_lk.*a_vec_BS_Comp_lk);
            H_Beam_Alignment_BS_angle_sq(kk_bs_ll,nn_bs,bs_ll_1) = Equ_a_BS_lk_init*Equ_a_AC_lk_init;
        end
    end
end
Adjust_num = 2;
[Est_miu_BS_set,Est_niu_BS_set,Est_Azi_BS_set,Est_Ele_BS_set] = A1_Est_BS_angle_V2(H_Beam_Alignment_BS_angle_sq,G_ls_init,alpha_init,...
    k_index_com,fs,fc,Power_Tx_ii,sigma_no_est,awgn_en,K_miu_BS_Re,K_miu_BS_Im,K_niu_BS_Re,K_niu_BS_Im,Ptx_gain_ii,...
    Delta_BS,miu_BS_set_init_error,niu_BS_set_init_error,I_BS_h_bar,I_BS_v_bar,Azi_BS_set_init_error,Ele_BS_set_init_error,...
    Adjust_num,s_BS_angle_set);

H_Beam_Alignment_AC_angle_sq = zeros(K_sub,N_OFDM_AC_angle,L);
f_RF_BS_init_set_sq = zeros(N_BS,L); a_vec_BS_Comp_est_set = zeros(N_BS,K_sub,L);
for ac_ll_1 = 1:L
    a_vec_BS_init_ll_est = A0_Array_Response_Vector(N_BS_h,N_BS_v,Est_Azi_BS_set(ac_ll_1),Est_Ele_BS_set(ac_ll_1));
    a_vec_AC_init_ll_error = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init_error(ac_ll_1),Ele_AC_set_init_error(ac_ll_1));
    f_RF_BS_init_ll = Quantize(a_vec_BS_init_ll_est,N_bits)/sqrt(N_BS);
    f_RF_BS_init_set_sq(:,ac_ll_1) = f_RF_BS_init_ll;
    Subarray_AC_ll_init = Quantize(a_vec_AC_init_ll_error(Index_W_RF_AC_Angle(:,1,ac_ll_1)),N_bits)/sqrt(M_AC_bar);
    Index_AC_Ant_Comp_ll = Index_AC_Ant_Comp(:,:,ac_ll_1);
    for kk_ac_ll = 1:K_sub
        a_vec_BS_Comp_lk = A0_Array_Response_Vector_Comp_BS(Index_BS_Ant_Comp,N_BS_h,N_BS_v,N_BS_Comp,N_BS_Comp_h,N_BS_Comp_v,M_BS_Comp,M_BS_Comp_h,...
            Est_Azi_BS_set(ac_ll_1),Est_Ele_BS_set(ac_ll_1),k_index_com(kk_ac_ll,ac_ll_1),K,fs,fc);
        a_vec_BS_Comp_est_set(:,kk_ac_ll,ac_ll_1) = a_vec_BS_Comp_lk;
        a_vec_AC_Comp_lk = A0_Array_Response_Vector_Comp_AC(ac_ll_1,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp_ll,N_AC_h,N_AC_v,M_AC_h,M_AC_v,...
            N_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp,M_AC_Comp_h,Azi_AC_set_init_error(ac_ll_1),Ele_AC_set_init_error(ac_ll_1),...
            k_index_com(kk_ac_ll,ac_ll_1),K,fs,fc);
        Equ_a_BS_lk_init = (a_vec_BS_init_be_sq_set(:,kk_ac_ll,ac_ll_1).*a_vec_BS_Comp_lk)'*f_RF_BS_init_ll;
        for nn_ac = 1:N_OFDM_AC_angle
            w_RF_AC_init_ll_nn = zeros(N_AC,1);
            w_RF_AC_init_ll_nn(Index_W_RF_AC_Angle(:,nn_ac,ac_ll_1)) = Subarray_AC_ll_init;
            Equ_a_AC_lk_init = w_RF_AC_init_ll_nn'*(a_vec_AC_init_be_sq_set(:,kk_ac_ll,ac_ll_1).*a_vec_AC_Comp_lk);
            H_Beam_Alignment_AC_angle_sq(kk_ac_ll,nn_ac,ac_ll_1) = Equ_a_AC_lk_init*Equ_a_BS_lk_init;
        end
    end
end
Adjust_num = 2;
[Est_miu_AC_set,Est_niu_AC_set,Est_Azi_AC_set,Est_Ele_AC_set] = A2_Est_AC_angle_V2(H_Beam_Alignment_AC_angle_sq,G_ls_init,alpha_init,...
    tau_init,fs,fc,k_index_com,Power_Tx_ii,sigma_no_est,awgn_en,K_miu_AC_Re,K_miu_AC_Im,K_niu_AC_Re,K_niu_AC_Im,Ptx_gain_ii,...
    Delta_AC,miu_AC_set_init_error,niu_AC_set_init_error,I_AC_h_bar,I_AC_v_bar,Azi_AC_set_init_error,Ele_AC_set_init_error,Adjust_num,s_AC_angle_set);

H_Beam_Alignment_Doppler_set_sq = zeros(K_sub,L);
w_RF_AC_init_set_sq = zeros(N_AC,L); a_vec_AC_Comp_est_set = zeros(N_AC,K_sub,L);
for dop_ll_1 = 1:L
    a_vec_AC_init_ll_est = A0_Array_Response_Vector(N_AC_h,N_AC_v,Est_Azi_AC_set(dop_ll_1),Est_Ele_AC_set(dop_ll_1));
    w_RF_AC_init_ll = zeros(N_AC,1);
    w_RF_AC_init_ll(Index_W_RF_Sub_Array(:,dop_ll_1)) = Quantize(a_vec_AC_init_ll_est(Index_W_RF_Sub_Array(:,dop_ll_1)),N_bits)/sqrt(M_AC);
    w_RF_AC_init_set_sq(:,dop_ll_1) = w_RF_AC_init_ll;
    for kk_dop_ll = 1:K_sub
        a_vec_AC_Comp_lk = A0_Array_Response_Vector_Comp_AC(dop_ll_1,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp_ll,N_AC_h,N_AC_v,M_AC_h,M_AC_v,...
            N_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp,M_AC_Comp_h,Est_Azi_AC_set(dop_ll_1),Est_Ele_AC_set(dop_ll_1),...
            k_index_com(kk_dop_ll,dop_ll_1),K,fs,fc);
        a_vec_AC_Comp_est_set(:,kk_dop_ll,dop_ll_1) = a_vec_AC_Comp_lk;
        Equ_a_BS_ll_init_dop = (a_vec_BS_init_be_sq_set(:,kk_dop_ll,dop_ll_1).*a_vec_BS_Comp_est_set(:,kk_dop_ll,dop_ll_1))'*...
            f_RF_BS_init_set_sq(:,dop_ll_1);
        Equ_a_AC_ll_init_dop = w_RF_AC_init_ll'*(a_vec_AC_init_be_sq_set(:,kk_dop_ll,dop_ll_1).*a_vec_AC_Comp_lk);
        H_Beam_Alignment_Doppler_set_sq(kk_dop_ll,dop_ll_1) = Equ_a_AC_ll_init_dop*Equ_a_BS_ll_init_dop;
    end
end
Adjust_num_Dop = 3;
[~,Est_Doppler_set_sq] = A3_Est_Dopple_V3(H_Beam_Alignment_Doppler_set_sq,G_ls_init,alpha_init,Doppler_init_k_set,...
    tau_init,fs,k_index_com,N_OFDM_Doppler,Power_Tx_ii,sigma_no_est,awgn_en,Ptx_gain_ii,T_sym,s_dop_set,Adjust_num_Dop,Doppler_z_init_error,lambda_z);

N_BS_Comp_h_del = 2; N_BS_Comp_v_del = N_BS_Comp_h_del;
N_BS_Comp_del = N_BS_Comp_h_del*N_BS_Comp_v_del;
M_BS_Comp_h_del = N_BS_h/N_BS_Comp_h_del; M_BS_Comp_v_del = N_BS_v/N_BS_Comp_v_del;
M_BS_Comp_del = M_BS_Comp_h_del*M_BS_Comp_v_del;
Index_BS_Ant_Comp_del = A0_Gene_Index_Antenna_Comp(N_BS_Comp_del,M_BS_Comp_del,N_BS_Comp_h_del,N_BS_Comp_v_del,M_BS_Comp_h_del,N_BS_h,N_BS_v);
N_AC_Comp_h_del = 2; N_AC_Comp_v_del = N_AC_Comp_h_del;
N_AC_Comp_del = N_AC_Comp_h_del*N_AC_Comp_v_del;
M_AC_Comp_h_del = M_AC_h/N_AC_Comp_h_del; M_AC_Comp_v_del = M_AC_v/N_AC_Comp_v_del;
M_AC_Comp_del = M_AC_Comp_h_del*M_AC_Comp_v_del;
Index_AC_Ant_Comp_del = A0_Gene_Index_Antenna_Comp(N_AC_Comp_del,M_AC_Comp_del,N_AC_Comp_h_del,N_AC_Comp_v_del,M_AC_Comp_h_del,N_AC_h,N_AC_v,...
    I_AC,I_AC_h,M_AC_h,M_AC_v);
H_Beam_Alignment_Delay_set_sq = zeros(K_sub,L);
for del_ll_1 = 1:L
    for kk_del1 = 1:K_sub
        a_vec_BS_Comp_lk = A0_Array_Response_Vector_Comp_BS(Index_BS_Ant_Comp_del,N_BS_h,N_BS_v,N_BS_Comp_del,N_BS_Comp_h_del,N_BS_Comp_v_del,...
            M_BS_Comp_del,M_BS_Comp_h_del,Est_Azi_BS_set(del_ll_1),Est_Ele_BS_set(del_ll_1),k_index_com(kk_del1,del_ll_1),K,fs,fc);
        a_vec_AC_Comp_lk = A0_Array_Response_Vector_Comp_AC(del_ll_1,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp_del(:,:,del_ll_1),...
            N_AC_h,N_AC_v,M_AC_h,M_AC_v,N_AC_Comp_del,N_AC_Comp_h_del,N_AC_Comp_v_del,M_AC_Comp_del,M_AC_Comp_h_del,...
            Est_Azi_AC_set(del_ll_1),Est_Ele_AC_set(del_ll_1),k_index_com(kk_del1,del_ll_1),K,fs,fc);
        Equ_a_BS_ll_init_del = (a_vec_BS_init_be_sq_set(:,kk_del1,del_ll_1).*a_vec_BS_Comp_lk)'*f_RF_BS_init_set_sq(:,del_ll_1);
        Equ_a_AC_ll_init_del = w_RF_AC_init_set_sq(:,del_ll_1)'*(a_vec_AC_init_be_sq_set(:,kk_del1,del_ll_1).*a_vec_AC_Comp_lk);
        H_Beam_Alignment_Delay_set_sq(kk_del1,del_ll_1) = Equ_a_AC_ll_init_del*Equ_a_BS_ll_init_del;
    end
end
radial_vt_BS_set_est = lambda_z*Est_Doppler_set_sq;
n_ofdm_delay = 1:N_OFDM_Delay;
Doppler_est_k_set = zeros(K_sub,L);
for ll_ddc = 1:L
    Doppler_est_k_set(:,ll_ddc) = Est_Doppler_set_sq(ll_ddc) + radial_vt_BS_set_est(ll_ddc)*((k_index_com(:,ll_ddc)-1)/K-0.5)*fs/3e8;
end
Phase_Delay_Dop_Comp_set_1 = zeros(L,N_OFDM_Delay,K_sub);
for kk_ddc = 1:K_sub
    Phase_Delay_Dop_Comp_set_1(:,:,kk_ddc) = 2*pi*(Doppler_init_k_set(kk_ddc,:)-Doppler_est_k_set(kk_ddc,:)).'*(n_ofdm_delay-1)*T_sym;
end
[Est_miu_tau_sq,Est_tau_set_sq,Y_tau_bar_set_sq] = A4_Est_Delay_V1(H_Beam_Alignment_Delay_set_sq,G_ls_init,alpha_init,Phase_Delay_Dop_Comp_set_1,...
    tau_init,fs,k_index_com,N_OFDM_Delay,Power_Tx_ii,sigma_no_est,awgn_en,D_subc,Ptx_gain_ii,s_tau_set);

Y_tau_hat_set = zeros(K_sub,N_OFDM_Delay,L);
for alpha_ll_1 = 1:L
    gamma_del_hat_ll = exp(1i*pi*fs*Est_tau_set_sq(alpha_ll_1))*sqrt(N_BS)*sqrt(M_AC);
    a_bar_Del_hat_ll = exp(1i*(k_index_com(:,alpha_ll_1)-1)*Est_miu_tau_sq(alpha_ll_1));
    Y_tau_hat_set(:,:,alpha_ll_1) = gamma_del_hat_ll*a_bar_Del_hat_ll*s_tau_set(:,alpha_ll_1).';
end
Est_alpha_set_sq = zeros(L,1);
for alpha_ll_2 = 1:L
    Est_alpha_set_sq(alpha_ll_2) = sum(sum(Y_tau_bar_set_sq(:,:,alpha_ll_2)./Y_tau_hat_set(:,:,alpha_ll_2)))/K_sub/N_OFDM_Delay;
end

%% 
turboEnc = comm.TurboEncoder('InterleaverIndicesSource','Input port');
turboDec = comm.TurboDecoder('InterleaverIndicesSource','Input port','NumIterations',4);
Mod_set = [2,4];

K_tmp = 3*floor(K/3);
K_tmp_index_set = (1:K_tmp).';
Num_Mod_init = 2; Mod_init = Mod_set(Num_Mod_init);
bps_init = log2(Mod_init);
EnData_perOFDM_init = Ns_BS*K_tmp*bps_init;
pktLen_init = (EnData_perOFDM_init-12)/3;

N_CCT = 500;
N_OFDM_C = 50;
T_CCT = N_OFDM_C*T_sym;
N_OFDM_tot = N_OFDM_C*N_CCT;
n_ofdm_data = 1:N_OFDM_C;
Amp_set_EDD = zeros(L,N_OFDM_tot);
Amp_set_EDD_CCEDC = zeros(L,N_OFDM_tot);
Amp_set_Real_eff = zeros(L,N_OFDM_tot);

%% 
epsilon = 0.2;
K_tilde = floor(K_tmp/2);
rho_alpha = randsrc*alpha_init/2;
rho_tau = -randsrc*tau_init/2;
rho_Doppler = -0.01*Doppler_z_init;
rho_theta_phi_AC = randsrc*pi/4;
rho_theta_phi_BS = randsrc*pi/12;

%% 
Doppler_k_set_init = zeros(K_tmp,L);
for ll_data_k1 = 1:L
    Doppler_k_set_init(:,ll_data_k1) = Doppler_z_init(ll_data_k1) + radial_vt_BS_set_init(ll_data_k1)*((K_tmp_index_set-1)/K-0.5)*fs/3e8;
end
Doppler_k_set_est = zeros(K_tmp,L);
for ll_data_k2 = 1:L
    Doppler_k_set_est(:,ll_data_k2) = Est_Doppler_set_sq(ll_data_k2) + radial_vt_BS_set_est(ll_data_k2)*((K_tmp_index_set-1)/K-0.5)*fs/3e8;
end
Phase_Data_Dop_Comp_init_set = zeros(L,N_OFDM_C,K_tmp);
for kk_data_1 = 1:K_tmp
    Phase_Data_Dop_Comp_init_set(:,:,kk_data_1) = 2*pi*(Doppler_k_set_init(kk_data_1,:)-Doppler_k_set_est(kk_data_1,:)).'*...
        ((1-1)*N_OFDM_C + n_ofdm_data-1)*T_sym;
end
h_vec_init_est_set_0 = zeros(K_tmp,L); h_vec_init_real_set = zeros(K_tmp,L,N_OFDM_C);
a_vec_BS_init_set = zeros(N_BS,K_tmp,L); a_vec_AC_init_set = zeros(N_AC,K_tmp,L);
a_vec_BS_Comp_set = zeros(N_BS,K_tmp,L); a_vec_AC_Comp_set = zeros(N_AC,K_tmp,L);
H_Beam_Alignment_init_sq = zeros(K_tmp,L);
for hh_ll_1 = 1:L
    gamma_h_hat_ll = Est_alpha_set_sq(hh_ll_1)*sqrt(N_BS)*sqrt(M_AC);
    h_vec_init_est_set_0(:,hh_ll_1) = gamma_h_hat_ll*exp(-1i*2*pi*Est_tau_set_sq(hh_ll_1)*(-1/2+(K_tmp_index_set-1)/K)*fs);
    for hh_kk_1 = 1:K_tmp
        a_vec_BS_init_lk = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init(hh_ll_1),Ele_BS_set_init(hh_ll_1),...
            K_tmp_index_set(hh_kk_1),K,fs,fc);
        a_vec_BS_init_set(:,hh_kk_1,hh_ll_1) = a_vec_BS_init_lk;
        a_vec_BS_Comp_lk = A0_Array_Response_Vector_Comp_BS(Index_BS_Ant_Comp,N_BS_h,N_BS_v,N_BS_Comp,N_BS_Comp_h,N_BS_Comp_v,M_BS_Comp,M_BS_Comp_h,...
            Est_Azi_BS_set(hh_ll_1),Est_Ele_BS_set(hh_ll_1),K_tmp_index_set(hh_kk_1),K,fs,fc);
        a_vec_BS_Comp_set(:,hh_kk_1,hh_ll_1) = a_vec_BS_Comp_lk;
        a_vec_AC_init_lk = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init(hh_ll_1),Ele_AC_set_init(hh_ll_1),...
            K_tmp_index_set(hh_kk_1),K,fs,fc);
        a_vec_AC_init_set(:,hh_kk_1,hh_ll_1) = a_vec_AC_init_lk;
        a_vec_AC_Comp_lk = A0_Array_Response_Vector_Comp_AC(hh_ll_1,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp(:,:,hh_ll_1),...
            N_AC_h,N_AC_v,M_AC_h,M_AC_v,N_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp,M_AC_Comp_h,...
            Est_Azi_AC_set(hh_ll_1),Est_Ele_AC_set(hh_ll_1),K_tmp_index_set(hh_kk_1),K,fs,fc);
        a_vec_AC_Comp_set(:,hh_kk_1,hh_ll_1) = a_vec_AC_Comp_lk;
        Equ_a_BS_lk_init_data = (a_vec_BS_init_lk.*a_vec_BS_Comp_lk)'*f_RF_BS_init_set_sq(:,hh_ll_1);
        Equ_a_AC_lk_init_data = w_RF_AC_init_set_sq(:,hh_ll_1)'*(a_vec_AC_init_lk.*a_vec_AC_Comp_lk);
        H_Beam_Alignment_init_sq(hh_kk_1,hh_ll_1) = Equ_a_AC_lk_init_data*Equ_a_BS_lk_init_data;
    end
    h_vec_init_real_set_temp = alpha_init(hh_ll_1)*exp(-1i*2*pi*tau_init(hh_ll_1)*(-1/2+(K_tmp_index_set-1)/K)*fs).*H_Beam_Alignment_init_sq(:,hh_ll_1);
    for n_ofdm_c_1 = 1:N_OFDM_C
        h_vec_init_real_set(:,hh_ll_1,n_ofdm_c_1) = reshape(exp(1i*Phase_Data_Dop_Comp_init_set(hh_ll_1,n_ofdm_c_1,:)),K_tmp,1).*h_vec_init_real_set_temp;
    end
end
result_h_vec_init = h_vec_init_est_set_0 - h_vec_init_real_set(:,:,1);
result_h_vec_init_amp = abs(h_vec_init_est_set_0) - abs(h_vec_init_real_set(:,:,1));
NMSE_eff_chan_temp = zeros(L,1);
for lll = 1:L
    NMSE_eff_chan_temp(lll) = norm(result_h_vec_init(:,lll))^2/norm(h_vec_init_real_set(:,lll,1))^2;
end
NMSE_eff_chan = sum(NMSE_eff_chan_temp)/L;
H_Beam_Alignment_init_interf = zeros(K_tmp,L); h_vec_init_interf = zeros(K_tmp,L,N_OFDM_C);
for hh_ll_2 = 1:L
    L_set = 1:L; L_set(hh_ll_2) = []; ll_set_re = L_set;
    for hh_kk_2 = 1:K_tmp
        Equ_a_BS_lk_init_interf = (a_vec_BS_init_set(:,hh_kk_2,ll_set_re).*a_vec_BS_Comp_set(:,hh_kk_2,ll_set_re))'*f_RF_BS_init_set_sq(:,ll_set_re); % 
        Equ_a_AC_lk_init_interf = w_RF_AC_init_set_sq(:,hh_ll_2)'*(a_vec_AC_init_set(:,hh_kk_2,ll_set_re).*a_vec_AC_Comp_set(:,hh_kk_2,hh_ll_2)); % 
        H_Beam_Alignment_init_interf(hh_kk_2,hh_ll_2) = Equ_a_AC_lk_init_interf*Equ_a_BS_lk_init_interf;
    end
    h_vec_init_interf_temp = alpha_init(ll_set_re)*exp(-1i*2*pi*tau_init(ll_set_re)*(-1/2+(K_tmp_index_set-1)/K)*fs).*H_Beam_Alignment_init_interf(:,hh_ll_2);
    for n_ofdm_c_2 = 1:N_OFDM_C
        h_vec_init_interf(:,hh_ll_2,n_ofdm_c_2) = reshape(exp(1i*Phase_Data_Dop_Comp_init_set(ll_set_re,n_ofdm_c_2,:)),K_tmp,1).*h_vec_init_interf_temp;
    end
end
sigma_no = sqrt(10^(-(-20/10)));
%% Start
for ii = 1:N_CCT
    tic
    Est_eff_Chan_DD_EDD_CCEDC_ii = zeros(K_tmp,L,2,N_OFDM_C+1);
    %% 
    if ii == 1
        h_vec_real_set_temp = h_vec_init_real_set;
        h_vec_interf_temp = h_vec_init_interf;
        for dd_1 = 1:2
            Est_eff_Chan_DD_EDD_CCEDC_ii(:,:,dd_1,1) = h_vec_init_est_set_0;
        end
    else
        Est_eff_Chan_DD_EDD_CCEDC_ii(:,:,:,1) = Est_eff_Chan_ii_temp;
        Azi_BS_set_temp = Azi_BS_set_init + rho_theta_phi_BS*(ii-1)*T_CCT;
        Ele_BS_set_temp = Ele_BS_set_init + rho_theta_phi_BS*(ii-1)*T_CCT;
        Azi_AC_set_temp = Azi_AC_set_init + rho_theta_phi_AC*(ii-1)*T_CCT;
        Ele_AC_set_temp = Ele_AC_set_init + rho_theta_phi_AC*(ii-1)*T_CCT;
        Doppler_z_temp = Doppler_z_init + rho_Doppler*(ii-1)*T_CCT;
        radial_vt_BS_set_temp = lambda_z*Doppler_z_temp;
        tau_temp = tau_init + rho_tau*(ii-1)*T_CCT;
        alpha_temp = alpha_init + rho_alpha*(ii-1)*T_CCT;
        Doppler_k_set_temp = zeros(K_tmp,L);
        for ll_tr_k1 = 1:L
            Doppler_k_set_temp(:,ll_tr_k1) = Doppler_z_temp(ll_tr_k1) + radial_vt_BS_set_temp(ll_tr_k1)*((K_tmp_index_set-1)/K-0.5)*fs/3e8;
        end
        Phase_Data_Dop_Comp_set_temp = zeros(L,N_OFDM_C,K_tmp);
        for kk_tr_1 = 1:K_tmp
            Phase_Data_Dop_Comp_set_temp(:,:,kk_tr_1) = 2*pi*(Doppler_k_set_temp(kk_tr_1,:)-Doppler_k_set_est(kk_tr_1,:)).'*...
                ((ii-1)*N_OFDM_C + n_ofdm_data-1)*T_sym;
        end
        h_vec_real_set_temp = zeros(K_tmp,L,N_OFDM_C);
        a_vec_BS_set_temp = zeros(N_BS,K_tmp,L); a_vec_AC_set_temp = zeros(N_AC,K_tmp,L);
        H_Beam_Alignment_sq_temp = zeros(K_tmp,L);
        for hh_tr_ll_1 = 1:L
            for hh_tr_kk_1 = 1:K_tmp
                a_vec_BS_init_lk = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_temp(hh_tr_ll_1),Ele_BS_set_temp(hh_tr_ll_1),...
                    K_tmp_index_set(hh_tr_kk_1),K,fs,fc);
                a_vec_BS_set_temp(:,hh_tr_kk_1,hh_tr_ll_1) = a_vec_BS_init_lk;
                a_vec_AC_init_lk = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_temp(hh_tr_ll_1),Ele_AC_set_temp(hh_tr_ll_1),...
                    K_tmp_index_set(hh_tr_kk_1),K,fs,fc);
                a_vec_AC_set_temp(:,hh_tr_kk_1,hh_tr_ll_1) = a_vec_AC_init_lk;
                Equ_a_BS_lk_init_data = (a_vec_BS_init_lk.*a_vec_BS_Comp_set(:,hh_tr_kk_1,hh_tr_ll_1))'*f_RF_BS_init_set_sq(:,hh_tr_ll_1);
                Equ_a_AC_lk_init_data = w_RF_AC_init_set_sq(:,hh_tr_ll_1)'*(a_vec_AC_init_lk.*a_vec_AC_Comp_set(:,hh_tr_kk_1,hh_tr_ll_1));
                H_Beam_Alignment_sq_temp(hh_tr_kk_1,hh_tr_ll_1) = Equ_a_AC_lk_init_data*Equ_a_BS_lk_init_data;
            end
            h_vec_init_real_set_temp = alpha_temp(hh_tr_ll_1)*exp(-1i*2*pi*tau_temp(hh_tr_ll_1)*(-1/2+(K_tmp_index_set-1)/K)*fs)...
                .*H_Beam_Alignment_sq_temp(:,hh_tr_ll_1);
            for n_ofdm_c_tr_1 = 1:N_OFDM_C
                h_vec_real_set_temp(:,hh_tr_ll_1,n_ofdm_c_tr_1) = reshape(exp(1i*Phase_Data_Dop_Comp_set_temp(hh_tr_ll_1,n_ofdm_c_tr_1,:)),K_tmp,1)...
                    .*h_vec_init_real_set_temp;
            end
        end
        H_Beam_Alignment_interf_temp = zeros(K_tmp,L); h_vec_interf_temp = zeros(K_tmp,L,N_OFDM_C);
        for hh_tr_ll_2 = 1:L
            L_set_tr = 1:L; L_set_tr(hh_tr_ll_2) = []; ll_set_re_tr = L_set_tr;
            for hh_tr_kk_2 = 1:K_tmp
                Equ_a_BS_lk_interf_temp = (a_vec_BS_set_temp(:,hh_tr_kk_2,ll_set_re_tr).*a_vec_BS_Comp_set(:,hh_tr_kk_2,ll_set_re_tr))'...
                    *f_RF_BS_init_set_sq(:,ll_set_re_tr); % 
                Equ_a_AC_lk_interf_temp = w_RF_AC_init_set_sq(:,hh_tr_ll_2)'*(a_vec_AC_set_temp(:,hh_tr_kk_2,ll_set_re_tr)...
                    .*a_vec_AC_Comp_set(:,hh_tr_kk_2,hh_tr_ll_2)); % 
                H_Beam_Alignment_interf_temp(hh_tr_kk_2,hh_tr_ll_2) = Equ_a_AC_lk_interf_temp*Equ_a_BS_lk_interf_temp;
            end
            h_vec_init_interf_temp_ll = alpha_temp(ll_set_re_tr)*exp(-1i*2*pi*tau_temp(ll_set_re_tr)*(-1/2+(K_tmp_index_set-1)/K)*fs)...
                .*H_Beam_Alignment_interf_temp(:,hh_tr_ll_2);
            for n_ofdm_c_tr_2 = 1:N_OFDM_C
                h_vec_interf_temp(:,hh_tr_ll_2,n_ofdm_c_tr_2) = reshape(exp(1i*Phase_Data_Dop_Comp_set_temp(ll_set_re_tr,n_ofdm_c_tr_2,:)),K_tmp,1)...
                    .*h_vec_init_interf_temp_ll;
            end
        end
    end
    %%
    for mm = 1:N_OFDM_C
        for ll_data_tr = 1:L
            Tx_data_ll = randi([0 1],pktLen_init,1);
            intrlvrInd_ll = randperm(pktLen_init);
            EnData_ll = turboEnc(Tx_data_ll,intrlvrInd_ll);
            Tx_Mod_data_ll = qammod(EnData_ll,Mod_init,'InputType','bit','UnitAveragePower',true);
            Rx_Symbols_set = zeros(K_tmp,1);
            Rx_Mod_data_ll_set = zeros(K_tmp,3);
            for kk_data_ll_3 = 1:K_tmp
                Rx_Symbols_lk = h_vec_real_set_temp(kk_data_ll_3,ll_data_tr,mm)*Tx_Mod_data_ll(kk_data_ll_3) + h_vec_interf_temp(kk_data_ll_3,ll_data_tr,mm)...
                    + awgn_en*sigma_no*(normrnd(0,1) + 1i*normrnd(0,1))/sqrt(2);
                Rx_Symbols_set(kk_data_ll_3) = Rx_Symbols_lk;
                for dd_2 = 1:2
                    Rx_Mod_data_ll_set(kk_data_ll_3,dd_2) = Est_eff_Chan_DD_EDD_CCEDC_ii(kk_data_ll_3,ll_data_tr,dd_2,mm)\Rx_Symbols_lk;
                end
            end
            Rx_demod_data_ll = qamdemod(Rx_Mod_data_ll_set,Mod_init,'UnitAveragePower',true,'OutputType','llr','NoiseVariance',sigma_no^2);
            
            Rx_Bits_EDD_ll = turboDec(-Rx_demod_data_ll(:,1),intrlvrInd_ll);
            Err_bits_EDD_ll = biterr(Tx_data_ll,Rx_Bits_EDD_ll);
            Re_EnData_EDD_ll = turboEnc(Rx_Bits_EDD_ll,intrlvrInd_ll);
            Re_Tx_Mod_data_EDD_ll = qammod(Re_EnData_EDD_ll,Mod_init,'InputType','bit','UnitAveragePower',true);
            Est_eff_Chan_DD_EDD_CCEDC_ii(:,ll_data_tr,1,mm+1) = Rx_Symbols_set./Re_Tx_Mod_data_EDD_ll;
            
            Rx_Bits_EDD_CCEDC_ll = turboDec(-Rx_demod_data_ll(:,2),intrlvrInd_ll);
            Err_bits_EDD_CCEDC_ll = biterr(Tx_data_ll,Rx_Bits_EDD_CCEDC_ll);
            Re_EnData_EDD_CCEDC_ll = turboEnc(Rx_Bits_EDD_CCEDC_ll,intrlvrInd_ll);
            Re_Tx_Mod_data_EDD_CCEDC_ll = qammod(Re_EnData_EDD_CCEDC_ll,Mod_init,'InputType','bit','UnitAveragePower',true);
            Est_eff_Chan_EDD_CCEDC_temp_ll = Rx_Symbols_set./Re_Tx_Mod_data_EDD_CCEDC_ll;
            kk_error_set_ll = [];
            Est_eff_Chan_EDD_CCEDC_ii_mm = Est_eff_Chan_DD_EDD_CCEDC_ii(:,ll_data_tr,2,mm);
            for k_ccedc_ll = 1:K_tmp
                abs_kk_residual = abs(Est_eff_Chan_EDD_CCEDC_temp_ll(k_ccedc_ll) - Est_eff_Chan_EDD_CCEDC_ii_mm(k_ccedc_ll));
                kk_threshold = epsilon*mean(abs(Est_eff_Chan_EDD_CCEDC_ii_mm));
                if abs_kk_residual > kk_threshold
                    kk_error_set_ll = [kk_error_set_ll,k_ccedc_ll];
                end
            end
            if length(kk_error_set_ll) > K_tilde
                error('Error: Channel parameter tracking need to be implemented !');
            end
            Est_eff_Chan_EDD_CCEDC_corr_ll = Est_eff_Chan_EDD_CCEDC_temp_ll;
            if ~isempty(kk_error_set_ll)
                K_tmp_index_temp_ll = 1:K_tmp; K_tmp_index_temp_ll(kk_error_set_ll) = [];
                for error_lk = 1:length(kk_error_set_ll)
                    Est_chan_coeff_mod_lk = (sum(abs(Est_eff_Chan_EDD_CCEDC_ii_mm)) + sum(abs(Est_eff_Chan_EDD_CCEDC_temp_ll(K_tmp_index_temp_ll))))...
                        /(2*K_tmp-length(kk_error_set_ll));
                    Est_chan_coeff_ang_lk = angle(Est_eff_Chan_EDD_CCEDC_ii_mm(error_lk));
                    Est_eff_Chan_EDD_CCEDC_corr_ll(error_lk) = Est_chan_coeff_mod_lk*exp(1i*Est_chan_coeff_ang_lk);
                end
            end
            Est_eff_Chan_DD_EDD_CCEDC_ii(:,ll_data_tr,2,mm+1) = Est_eff_Chan_EDD_CCEDC_corr_ll;
        end
        if mm == N_OFDM_C
            Est_eff_Chan_ii_temp = Est_eff_Chan_DD_EDD_CCEDC_ii(:,:,:,mm+1);
        end
        %%
        if mod(mm,5) == 0
            toc
            disp(['  N_CCT = ' num2str(ii) ',  n_OFDM = ' num2str(mm)])
        end
    end
    for ll_data_amp = 1:L
        for mm_amp = 1:N_OFDM_C
            Amp_set_EDD(ll_data_amp,(ii-1)*N_OFDM_C+mm_amp) = mean(abs(Est_eff_Chan_DD_EDD_CCEDC_ii(:,ll_data_amp,1,mm_amp+1)));
            Amp_set_EDD_CCEDC(ll_data_amp,(ii-1)*N_OFDM_C+mm_amp) = mean(abs(Est_eff_Chan_DD_EDD_CCEDC_ii(:,ll_data_amp,2,mm_amp+1)));
            Amp_set_Real_eff(ll_data_amp,(ii-1)*N_OFDM_C+mm_amp) = mean(abs(h_vec_real_set_temp(:,ll_data_amp,mm_amp)));
        end
    end
    disp(['  Finish N_CCT = ' num2str(ii)])
end
h_vec_init_est_set_amp = mean(abs(h_vec_init_est_set_0)).';
Amp_set_notrack = repmat(h_vec_init_est_set_amp,1,N_OFDM_tot);

disp('Finish  All')

N_OFDM_tot = N_OFDM_C*N_CCT;
n_OFDM_tot = 1:N_OFDM_tot;
MarkerSize = 3;
LineWidth = 1.2;
Fontsize = 15;

figure
plot(n_OFDM_tot,Amp_set_notrack(1,:),'--b','LineWidth',LineWidth); hold on; grid on;
plot(n_OFDM_tot,Amp_set_EDD(1,:),'--r','LineWidth',LineWidth+1);
plot(n_OFDM_tot,Amp_set_EDD_CCEDC(1,:),'--c','LineWidth',LineWidth);
plot(n_OFDM_tot,Amp_set_Real_eff(1,:),'--k','LineWidth',LineWidth+0.5);
plot(n_OFDM_tot,Amp_set_notrack(2,:),'-b','LineWidth',LineWidth);
plot(n_OFDM_tot,Amp_set_EDD(2,:),'-r','LineWidth',LineWidth+1);
plot(n_OFDM_tot,Amp_set_EDD_CCEDC(2,:),'-c','LineWidth',LineWidth);
plot(n_OFDM_tot,Amp_set_Real_eff(2,:),'-k','LineWidth',LineWidth+0.5);
xlabel('Number of CCT','Fontsize',Fontsize),ylabel('Amplitude','Fontsize',Fontsize);
title('Amplitude vs CCT','Fontsize',Fontsize);
set(gca, 'GridLineStyle', '-.','FontSize',Fontsize, 'linewidth',1.5,'Fontname','Times New Roman');
set(gcf, 'position', [700 300 650 550]); axis normal;
h1 = legend('BS1-No tracking','BS1-DADD','BS1-DADD CCEDC','BS1-Real Chan',...
    'BS2-No tracking','BS2-DADD','BS2-DADD CCEDC','BS2-Real Chan','Location','southwest');
set(h1,'Fontsize',11);

N_CCT_used = ii - 1;
N_OFDM_tot_used = N_OFDM_C*N_CCT_used;

h_vec_init_est_set_amp_used = mean(abs(h_vec_init_est_set_0)).';
Amp_set_notrack_used = repmat(h_vec_init_est_set_amp_used,1,N_OFDM_tot_used);
Amp_set_EDD_used = Amp_set_EDD(:,1:N_OFDM_tot_used);
Amp_set_EDD_CCEDC_used = Amp_set_EDD_CCEDC(:,1:N_OFDM_tot_used);
Amp_set_Real_eff_used = Amp_set_Real_eff(:,1:N_OFDM_tot_used);

n_OFDM_tot_used = 1:N_OFDM_tot_used;
MarkerSize = 3;
LineWidth = 1.2;
Fontsize = 15;

figure
plot(n_OFDM_tot_used,Amp_set_notrack_used(1,:),'--b','LineWidth',LineWidth); hold on; grid on;
plot(n_OFDM_tot_used,Amp_set_EDD_used(1,:),'--r','LineWidth',LineWidth+1);
plot(n_OFDM_tot_used,Amp_set_EDD_CCEDC_used(1,:),'--c','LineWidth',LineWidth);
plot(n_OFDM_tot_used,Amp_set_Real_eff_used(1,:),'--k','LineWidth',LineWidth+0.5);
plot(n_OFDM_tot_used,Amp_set_notrack_used(2,:),'-b','LineWidth',LineWidth);
plot(n_OFDM_tot_used,Amp_set_EDD_used(2,:),'-r','LineWidth',LineWidth+1);
plot(n_OFDM_tot_used,Amp_set_EDD_CCEDC_used(2,:),'-c','LineWidth',LineWidth);
plot(n_OFDM_tot_used,Amp_set_Real_eff_used(2,:),'-k','LineWidth',LineWidth+0.5);
xlabel('Number of CCT','Fontsize',Fontsize),ylabel('Amplitude','Fontsize',Fontsize);
title('Amplitude vs CCT','Fontsize',Fontsize);
set(gca, 'GridLineStyle', '-.','FontSize',Fontsize, 'linewidth',1.5,'Fontname','Times New Roman');
set(gcf, 'position', [700 300 650 550]); axis normal;
h2 = legend('BS1-No tracking','BS1-DADD','BS1-DADD CCEDC','BS1-Real Chan',...
    'BS2-No tracking','BS2-DADD','BS2-DADD CCEDC','BS2-Real Chan','Location','southwest');
set(h2,'Fontsize',11);
