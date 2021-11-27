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
Distance_init(1) = norm(Coo_BS1_A-Coo_Air_C_init); Distance_init(2) = norm(Coo_BS2_B-Coo_Air_C_init);
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
m_CCT_init = 1;
Doppler_diff = Doppler_z_init*1e-2;
Doppler_z_init_error = Doppler_z_init+Doppler_diff;
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
Power_Tx_dBm = 0:5:150;
Ptx_gain = ones(1,length(Power_Tx_dBm));
Power_Tx_ii = (1e-3*10^(Power_Tx_dBm(1)/10));
Ptx_gain_ii = Ptx_gain(1);
H_Beam_Alignment_BS_angle_0 = zeros(N_OFDM_BS_angle,L);
a_vec_AC_init_set = zeros(N_AC,L); a_vec_BS_init_set = zeros(N_BS,L);
for bs_ll_1 = 1:L
    a_vec_AC_init_ll = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init(bs_ll_1),Ele_AC_set_init(bs_ll_1));
    a_vec_BS_init_ll = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init(bs_ll_1),Ele_BS_set_init(bs_ll_1));
    a_vec_AC_init_set(:,bs_ll_1) = a_vec_AC_init_ll; a_vec_BS_init_set(:,bs_ll_1) = a_vec_BS_init_ll;
    a_vec_AC_init_ll_error = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init_error(bs_ll_1),Ele_AC_set_init_error(bs_ll_1));
    a_vec_BS_init_ll_error = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init_error(bs_ll_1),Ele_BS_set_init_error(bs_ll_1));
    f_RF_AC_init_ll = zeros(N_AC,1);
    f_RF_AC_init_ll(Index_W_RF_Sub_Array(:,bs_ll_1)) = Quantize(a_vec_AC_init_ll_error(Index_W_RF_Sub_Array(:,bs_ll_1)),N_bits)/sqrt(M_AC);
    Equ_a_AC_ll_init = a_vec_AC_init_ll'*f_RF_AC_init_ll;
    Subarray_BS_ll_init = Quantize(a_vec_BS_init_ll_error(Index_W_RF_BS_Angle(:,1)),N_bits)/sqrt(N_BS_bar);
    for bs_nn_1 = 1:N_OFDM_BS_angle
        w_RF_BS_init_ll_nn = zeros(N_BS,1);
        w_RF_BS_init_ll_nn(Index_W_RF_BS_Angle(:,bs_nn_1)) = Subarray_BS_ll_init;
        Equ_a_BS_ll_init = w_RF_BS_init_ll_nn'*a_vec_BS_init_ll;
        H_Beam_Alignment_BS_angle_0(bs_nn_1,bs_ll_1) = Equ_a_BS_ll_init*Equ_a_AC_ll_init;
    end
end
H_Beam_Alignment_BS_angle_BeSq = zeros(K_sub,N_OFDM_BS_angle,L);
a_vec_AC_init_besq_set = zeros(N_AC,K_sub,L); a_vec_BS_init_besq_set = zeros(N_BS,K_sub,L);
for bs_ll_2 = 1:L
    a_vec_AC_init_ll_error = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init_error(bs_ll_2),Ele_AC_set_init_error(bs_ll_2));
    a_vec_BS_init_ll_error = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init_error(bs_ll_2),Ele_BS_set_init_error(bs_ll_2));
    f_RF_AC_init_ll = zeros(N_AC,1);
    f_RF_AC_init_ll(Index_W_RF_Sub_Array(:,bs_ll_2)) = Quantize(a_vec_AC_init_ll_error(Index_W_RF_Sub_Array(:,bs_ll_2)),N_bits)/sqrt(M_AC);
    Subarray_BS_ll_init = Quantize(a_vec_BS_init_ll_error(Index_W_RF_BS_Angle(:,1)),N_bits)/sqrt(N_BS_bar);
    Index_AC_Ant_Comp_ll = Index_AC_Ant_Comp(:,:,bs_ll_2);
    for bs_kk_ll = 1:K_sub
        a_vec_AC_init_lk = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init(bs_ll_2),Ele_AC_set_init(bs_ll_2),...
            k_index_com(bs_kk_ll,bs_ll_2),K,fs,fc);
        a_vec_AC_Comp_lk = A0_Array_Response_Vector_Comp_AC(bs_ll_2,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp_ll,N_AC_h,N_AC_v,M_AC_h,M_AC_v,...
            N_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp,M_AC_Comp_h,Azi_AC_set_init_error(bs_ll_2),Ele_AC_set_init_error(bs_ll_2),...
            k_index_com(bs_kk_ll,bs_ll_2),K,fs,fc);
        a_vec_BS_init_lk = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init(bs_ll_2),Ele_BS_set_init(bs_ll_2),...
            k_index_com(bs_kk_ll,bs_ll_2),K,fs,fc);
        a_vec_BS_Comp_lk = A0_Array_Response_Vector_Comp_BS(Index_BS_Ant_Comp,N_BS_h,N_BS_v,N_BS_Comp,N_BS_Comp_h,N_BS_Comp_v,M_BS_Comp,M_BS_Comp_h,...
            Azi_BS_set_init_error(bs_ll_2),Ele_BS_set_init_error(bs_ll_2),k_index_com(bs_kk_ll,bs_ll_2),K,fs,fc);
        a_vec_AC_init_besq_set(:,bs_kk_ll,bs_ll_2) = a_vec_AC_init_lk;
        a_vec_BS_init_besq_set(:,bs_kk_ll,bs_ll_2) = a_vec_BS_init_lk;
        Equ_a_AC_lk_init = (a_vec_AC_init_lk.*a_vec_AC_Comp_lk)'*f_RF_AC_init_ll;
        for bs_nn_2 = 1:N_OFDM_BS_angle
            w_RF_BS_init_ll_nn = zeros(N_BS,1);
            w_RF_BS_init_ll_nn(Index_W_RF_BS_Angle(:,bs_nn_2)) = Subarray_BS_ll_init;
            Equ_a_BS_lk_init = w_RF_BS_init_ll_nn'*(a_vec_BS_init_lk.*a_vec_BS_Comp_lk);
            H_Beam_Alignment_BS_angle_BeSq(bs_kk_ll,bs_nn_2,bs_ll_2) = Equ_a_BS_lk_init*Equ_a_AC_lk_init;
        end
    end
end

%% set snr
iterMax = 2e3;
SNR_dBs = -35:3:-5;
ASE_Per_CSI_No_TrSq = zeros(length(SNR_dBs),1);
ASE_Per_CSI_TrSq = zeros(length(SNR_dBs),1);
ASE_Est_CSI_No_TrSq_temp = zeros(length(SNR_dBs),1);
ASE_Est_CSI_TrSq_temp = zeros(length(SNR_dBs),1);

for ii = 1:length(SNR_dBs)
    sigma_no = sqrt(10^(-(SNR_dBs(ii)/10)));
    tic
    for iter = 1:iterMax
        [~,~,Est_Azi_BS_set_0,Est_Ele_BS_set_0] = A1_Est_BS_angle_V1(H_Beam_Alignment_BS_angle_0,G_ls_init,alpha_init,...
            k_index_com,Power_Tx_ii,sigma_no,awgn_en,K_miu_BS_Re, K_miu_BS_Im, K_niu_BS_Re, K_niu_BS_Im, ...
            Ptx_gain_ii,Delta_BS,miu_BS_set_init_error,niu_BS_set_init_error,s_BS_angle_set);
        
        Adjust_num = 2;
        [~,~,Est_Azi_BS_set_1,Est_Ele_BS_set_1] = A1_Est_BS_angle_V2(H_Beam_Alignment_BS_angle_BeSq,G_ls_init,alpha_init,...
            k_index_com,fs,fc,Power_Tx_ii,sigma_no,awgn_en,K_miu_BS_Re,K_miu_BS_Im,K_niu_BS_Re,K_niu_BS_Im,Ptx_gain_ii,...
            Delta_BS,miu_BS_set_init_error,niu_BS_set_init_error,I_BS_h_bar,I_BS_v_bar,Azi_BS_set_init_error,Ele_BS_set_init_error,...
            Adjust_num,s_BS_angle_set);
        
        H_Beam_Alignment_AC_angle_0 = zeros(N_OFDM_AC_angle,L); f_RF_BS_init_set_0 = zeros(N_BS,L); a_vec_BS_init_set_est = zeros(N_BS,L);
        for ac_ll_1 = 1:L
            a_vec_BS_init_ll_est = A0_Array_Response_Vector(N_BS_h,N_BS_v,Est_Azi_BS_set_0(ac_ll_1),Est_Ele_BS_set_0(ac_ll_1));
            a_vec_AC_init_ll_error = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init_error(ac_ll_1),Ele_AC_set_init_error(ac_ll_1));
            a_vec_BS_init_set_est(:,ac_ll_1) = a_vec_BS_init_ll_est;
            f_RF_BS_init_ll = Quantize(a_vec_BS_init_ll_est,N_bits)/sqrt(N_BS);
            f_RF_BS_init_set_0(:,ac_ll_1) = f_RF_BS_init_ll;
            Equ_a_BS_ll_init = a_vec_BS_init_set(:,ac_ll_1)'*f_RF_BS_init_ll;
            Subarray_AC_ll_init = Quantize(a_vec_AC_init_ll_error(Index_W_RF_AC_Angle(:,1,ac_ll_1)),N_bits)/sqrt(M_AC_bar); % 
            for ac_nn_1 = 1:N_OFDM_AC_angle
                w_RF_AC_init_ll_nn = zeros(N_AC,1);
                w_RF_AC_init_ll_nn(Index_W_RF_AC_Angle(:,ac_nn_1,ac_ll_1)) = Subarray_AC_ll_init;
                Equ_a_AC_ll_init = w_RF_AC_init_ll_nn'*a_vec_AC_init_set(:,ac_ll_1);
                H_Beam_Alignment_AC_angle_0(ac_nn_1,ac_ll_1) = Equ_a_AC_ll_init*Equ_a_BS_ll_init;
            end
        end
        [~,~,Est_Azi_AC_set_0,Est_Ele_AC_set_0] = A2_Est_AC_angle_V1(H_Beam_Alignment_AC_angle_0,G_ls_init,alpha_init,...
            tau_init,fs,k_index_com,Power_Tx_ii,sigma_no,awgn_en,K_miu_AC_Re,K_miu_AC_Im,K_niu_AC_Re,K_niu_AC_Im,...
            Ptx_gain_ii,Delta_AC,miu_AC_set_init_error,niu_AC_set_init_error,s_AC_angle_set);
        
        H_Beam_Alignment_AC_angle_BeSq = zeros(K_sub,N_OFDM_AC_angle,L);
        f_RF_BS_init_set_BeSq = zeros(N_BS,L); a_vec_BS_Comp_est_set = zeros(N_BS,K_sub,L);
        for ac_ll_2 = 1:L
            a_vec_BS_init_ll_est = A0_Array_Response_Vector(N_BS_h,N_BS_v,Est_Azi_BS_set_1(ac_ll_2),Est_Ele_BS_set_1(ac_ll_2));
            a_vec_AC_init_ll_error = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init_error(ac_ll_2),Ele_AC_set_init_error(ac_ll_2));
            f_RF_BS_init_ll = Quantize(a_vec_BS_init_ll_est,N_bits)/sqrt(N_BS);
            f_RF_BS_init_set_BeSq(:,ac_ll_2) = f_RF_BS_init_ll;
            Subarray_AC_ll_init = Quantize(a_vec_AC_init_ll_error(Index_W_RF_AC_Angle(:,1,ac_ll_2)),N_bits)/sqrt(M_AC_bar);
            Index_AC_Ant_Comp_ll = Index_AC_Ant_Comp(:,:,ac_ll_2);
            for ac_kk_ll = 1:K_sub
                a_vec_BS_Comp_lk = A0_Array_Response_Vector_Comp_BS(Index_BS_Ant_Comp,N_BS_h,N_BS_v,N_BS_Comp,N_BS_Comp_h,N_BS_Comp_v,M_BS_Comp,M_BS_Comp_h,...
                    Est_Azi_BS_set_1(ac_ll_2),Est_Ele_BS_set_1(ac_ll_2),k_index_com(ac_kk_ll,ac_ll_2),K,fs,fc);
                a_vec_BS_Comp_est_set(:,ac_kk_ll,ac_ll_2) = a_vec_BS_Comp_lk;
                a_vec_AC_Comp_lk = A0_Array_Response_Vector_Comp_AC(ac_ll_2,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp_ll,N_AC_h,N_AC_v,M_AC_h,M_AC_v,...
                    N_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp,M_AC_Comp_h,Azi_AC_set_init_error(ac_ll_2),Ele_AC_set_init_error(ac_ll_2),...
                    k_index_com(ac_kk_ll,ac_ll_2),K,fs,fc);
                Equ_a_BS_lk_init = (a_vec_BS_init_besq_set(:,ac_kk_ll,ac_ll_2).*a_vec_BS_Comp_lk)'*f_RF_BS_init_ll;
                for ac_nn_2 = 1:N_OFDM_AC_angle
                    w_RF_AC_init_ll_nn = zeros(N_AC,1);
                    w_RF_AC_init_ll_nn(Index_W_RF_AC_Angle(:,ac_nn_2,ac_ll_2)) = Subarray_AC_ll_init;
                    Equ_a_AC_lk_init = w_RF_AC_init_ll_nn'*(a_vec_AC_init_besq_set(:,ac_kk_ll,ac_ll_2).*a_vec_AC_Comp_lk);
                    H_Beam_Alignment_AC_angle_BeSq(ac_kk_ll,ac_nn_2,ac_ll_2) = Equ_a_AC_lk_init*Equ_a_BS_lk_init;
                end
            end
        end
        Adjust_num_ac = 2;
        [~,~,Est_Azi_AC_set_1,Est_Ele_AC_set_1] = A2_Est_AC_angle_V2(H_Beam_Alignment_AC_angle_BeSq,G_ls_init,alpha_init,...
            tau_init,fs,fc,k_index_com,Power_Tx_ii,sigma_no,awgn_en,K_miu_AC_Re,K_miu_AC_Im,K_niu_AC_Re,K_niu_AC_Im,Ptx_gain_ii,...
            Delta_AC,miu_AC_set_init_error,niu_AC_set_init_error,I_AC_h_bar,I_AC_v_bar,Azi_AC_set_init_error,Ele_AC_set_init_error,Adjust_num_ac,s_AC_angle_set);
        
        H_Beam_Alignment_Doppler_0 = zeros(1,L); w_RF_AC_init_set_0 = zeros(N_AC,L); a_vec_AC_init_set_est = zeros(N_AC,L);
        for dop_ll_1 = 1:L
            a_vec_AC_init_ll_est = A0_Array_Response_Vector(N_AC_h,N_AC_v,Est_Azi_AC_set_0(dop_ll_1),Est_Ele_AC_set_0(dop_ll_1));
            a_vec_AC_init_set_est(:,dop_ll_1) = a_vec_AC_init_ll_est;
            Equ_a_BS_ll_init = a_vec_BS_init_set(:,dop_ll_1)'*f_RF_BS_init_set_0(:,dop_ll_1);
            w_RF_AC_init_ll = zeros(N_AC,1);
            w_RF_AC_init_ll(Index_W_RF_Sub_Array(:,dop_ll_1)) = Quantize(a_vec_AC_init_ll_est(Index_W_RF_Sub_Array(:,dop_ll_1)),N_bits)/sqrt(M_AC);
            w_RF_AC_init_set_0(:,dop_ll_1) = w_RF_AC_init_ll;
            Equ_a_AC_ll_init = w_RF_AC_init_ll'*a_vec_AC_init_set(:,dop_ll_1);
            H_Beam_Alignment_Doppler_0(:,dop_ll_1) = Equ_a_AC_ll_init*Equ_a_BS_ll_init;
        end
        [~,Est_Doppler_set_0] = A3_Est_Dopple_V1(H_Beam_Alignment_Doppler_0,G_ls_init,alpha_init,Phase_Doppler_Dop_set,...
            tau_init,fs,k_index_com,N_OFDM_Doppler,Power_Tx_ii,sigma_no,awgn_en,Ptx_gain_ii,T_sym,s_dop_set);
        
        H_Beam_Alignment_Doppler_BeSq = zeros(K_sub,L);
        w_RF_AC_init_set_BeSq = zeros(N_AC,L);
        for dop_ll_2 = 1:L
            a_vec_AC_init_ll_est = A0_Array_Response_Vector(N_AC_h,N_AC_v,Est_Azi_AC_set_1(dop_ll_2),Est_Ele_AC_set_1(dop_ll_2));
            w_RF_AC_init_ll = zeros(N_AC,1);
            w_RF_AC_init_ll(Index_W_RF_Sub_Array(:,dop_ll_2)) = Quantize(a_vec_AC_init_ll_est(Index_W_RF_Sub_Array(:,dop_ll_2)),N_bits)/sqrt(M_AC);
            w_RF_AC_init_set_BeSq(:,dop_ll_2) = w_RF_AC_init_ll;
            for dop_ll_kk = 1:K_sub
                a_vec_AC_Comp_lk = A0_Array_Response_Vector_Comp_AC(dop_ll_2,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp_ll,N_AC_h,N_AC_v,M_AC_h,M_AC_v,...
                    N_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp,M_AC_Comp_h,Est_Azi_AC_set_1(dop_ll_2),Est_Ele_AC_set_1(dop_ll_2),...
                    k_index_com(dop_ll_kk,dop_ll_2),K,fs,fc);
                Equ_a_BS_ll_init_dop = (a_vec_BS_init_besq_set(:,dop_ll_kk,dop_ll_2).*a_vec_BS_Comp_est_set(:,dop_ll_kk,dop_ll_2))'*...
                    f_RF_BS_init_set_BeSq(:,dop_ll_2);
                Equ_a_AC_ll_init_dop = w_RF_AC_init_ll'*(a_vec_AC_init_besq_set(:,dop_ll_kk,dop_ll_2).*a_vec_AC_Comp_lk);
                H_Beam_Alignment_Doppler_BeSq(dop_ll_kk,dop_ll_2) = Equ_a_AC_ll_init_dop*Equ_a_BS_ll_init_dop;
            end
        end
        Adjust_num_dop = 3;
        [~,Est_Doppler_set_1] = A3_Est_Dopple_V3(H_Beam_Alignment_Doppler_BeSq,G_ls_init,alpha_init,Doppler_init_k_set,...
            tau_init,fs,k_index_com,N_OFDM_Doppler,Power_Tx_ii,sigma_no,awgn_en,Ptx_gain_ii,T_sym,s_dop_set,Adjust_num_dop,Doppler_z_init_error,lambda_z);
        
        %%
        if ii == 1 && iter == 1
            h_vec_init_real_0 = zeros(K,L);
            for hh_ll_r1 = 1:L
                gamma_h_real_0_ll = alpha_init(hh_ll_r1)*sqrt(N_BS)*sqrt(M_AC);
                h_vec_init_real_0(:,hh_ll_r1) = gamma_h_real_0_ll*exp(-1i*2*pi*tau_init(hh_ll_r1)*(-1/2+(K_index_set-1)/K)*fs);
            end
            h_vec_init_real_interf_0 = zeros(K,L); f_RF_BS_init_real_set = zeros(N_BS,L); w_RF_AC_init_real_set = zeros(N_AC,L);
            for hh_ll_r2 = 1:L
                L_set_r2 = 1:L; L_set_r2(hh_ll_r2) = []; ll_set_re_r2 = L_set_r2;
                f_RF_BS_init_real_ll_re = Quantize(a_vec_BS_init_set(:,ll_set_re_r2),N_bits)/sqrt(N_BS);
                f_RF_BS_init_real_set(:,ll_set_re_r2) = f_RF_BS_init_real_ll_re;
                a_vec_AC_init_real_ll = a_vec_AC_init_set(:,hh_ll_r2);
                w_RF_AC_init_real_ll = zeros(N_AC,1);
                w_RF_AC_init_real_ll(Index_W_RF_Sub_Array(:,hh_ll_r2)) = Quantize(a_vec_AC_init_real_ll(Index_W_RF_Sub_Array(:,hh_ll_r2)),N_bits)/sqrt(M_AC);
                w_RF_AC_init_real_set(:,hh_ll_r2) = w_RF_AC_init_real_ll;
                Equ_a_BS_real_0_ll_interf = a_vec_BS_init_set(:,ll_set_re_r2)'*f_RF_BS_init_real_ll_re;
                Equ_a_AC_real_0_ll_interf = w_RF_AC_init_real_ll'*a_vec_AC_init_set(:,ll_set_re_r2);
                gamma_h_real_0_ll_interf = alpha_init(ll_set_re_r2)*Equ_a_AC_real_0_ll_interf*Equ_a_BS_real_0_ll_interf;
                h_vec_init_real_interf_0(:,hh_ll_r2) = gamma_h_real_0_ll_interf*exp(-1i*2*pi*tau_init(ll_set_re_r2)*(-1/2+(K_index_set-1)/K)*fs);
            end
        end
        if iter == 1
            for ase_ll_r0 = 1:L
                ASE_Per_CSI_No_TrSq_temp_ll = zeros(K,1);
                for ase_lk_r0 = 1:K
                    ASE_Per_CSI_No_TrSq_temp_ll(ase_lk_r0) = log2(1 + (abs(h_vec_init_real_0(ase_lk_r0,ase_ll_r0))^2/...
                        (abs(h_vec_init_real_interf_0(ase_lk_r0,ase_ll_r0))^2 + sigma_no^2*norm(w_RF_AC_init_real_set(:,ase_ll_r0))^2)));
                end
                ASE_Per_CSI_No_TrSq(ii) = ASE_Per_CSI_No_TrSq(ii) + mean(ASE_Per_CSI_No_TrSq_temp_ll);
            end
        end
        h_vec_init_est_0 = zeros(K,L);
        for hh_ll_e1 = 1:L
            Equ_a_BS_est_ll = a_vec_BS_init_set(:,hh_ll_e1)'*f_RF_BS_init_set_0(:,hh_ll_e1);
            Equ_a_AC_est_ll = w_RF_AC_init_set_0(:,hh_ll_e1)'*a_vec_AC_init_set(:,hh_ll_e1);
            gamma_h_est_0_ll = alpha_init(hh_ll_e1)*exp(1i*2*pi*(Doppler_z_init(hh_ll_e1)-Est_Doppler_set_0(hh_ll_e1))*T_sym)*Equ_a_AC_est_ll*Equ_a_BS_est_ll;
            h_vec_init_est_0(:,hh_ll_e1) = gamma_h_est_0_ll*exp(-1i*2*pi*tau_init(hh_ll_e1)*(-1/2+(K_index_set-1)/K)*fs);
        end
        h_vec_init_est_interf_0 = zeros(K,L);
        for hh_ll_e2 = 1:L
            L_set_e2 = 1:L; L_set_e2(hh_ll_e2) = []; ll_set_re_e2 = L_set_e2;
            Equ_a_BS_est_ll_interf = a_vec_BS_init_set(:,ll_set_re_e2)'*f_RF_BS_init_set_0(:,ll_set_re_e2);
            Equ_a_AC_est_ll_interf = w_RF_AC_init_set_0(:,hh_ll_e2)'*a_vec_AC_init_set(:,ll_set_re_e2);
            gamma_h_est_0_ll_interf = alpha_init(ll_set_re_e2)*exp(1i*2*pi*(Doppler_z_init(ll_set_re_e2)-Est_Doppler_set_0(ll_set_re_e2))*T_sym)*...
                Equ_a_AC_est_ll_interf*Equ_a_BS_est_ll_interf;
            h_vec_init_est_interf_0(:,hh_ll_e2) = gamma_h_est_0_ll_interf*exp(-1i*2*pi*tau_init(ll_set_re_e2)*(-1/2+(K_index_set-1)/K)*fs);
        end
        for ase_ll_e0 = 1:L
            ASE_Est_CSI_No_TrSq_temp_ll = zeros(K,1);
            for ase_lk_e0 = 1:K
                ASE_Est_CSI_No_TrSq_temp_ll(ase_lk_e0) = log2(1 + (abs(h_vec_init_est_0(ase_lk_e0,ase_ll_e0))^2/...
                    (abs(h_vec_init_est_interf_0(ase_lk_e0,ase_ll_e0))^2 + sigma_no^2*norm(w_RF_AC_init_set_0(:,ase_ll_e0))^2)));
            end
            ASE_Est_CSI_No_TrSq_temp(ii) = ASE_Est_CSI_No_TrSq_temp(ii) + mean(ASE_Est_CSI_No_TrSq_temp_ll);
        end
        if ii == 1 && iter == 1
            h_vec_init_real_1 = zeros(K,L);
            a_vec_BS_init_set_TrSq = zeros(N_BS,K,L); a_vec_AC_init_set_TrSq = zeros(N_AC,K,L);
            a_vec_BS_Comp_real_set = zeros(N_BS,K,L); a_vec_AC_Comp_real_set = zeros(N_AC,K,L);
            for hh_ll_rts1 = 1:L
                h_vec_init_real_ll_temp = zeros(K,1);
                for hh_kk_rts1 = 1:K
                    a_vec_BS_init_lk = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init(hh_ll_rts1),Ele_BS_set_init(hh_ll_rts1),...
                        K_index_set(hh_kk_rts1),K,fs,fc);
                    a_vec_BS_init_set_TrSq(:,hh_kk_rts1,hh_ll_rts1) = a_vec_BS_init_lk;
                    a_vec_BS_Comp_lk = A0_Array_Response_Vector_Comp_BS(Index_BS_Ant_Comp,N_BS_h,N_BS_v,N_BS_Comp,N_BS_Comp_h,N_BS_Comp_v,...
                        M_BS_Comp,M_BS_Comp_h,Azi_BS_set_init(hh_ll_rts1),Ele_BS_set_init(hh_ll_rts1),K_index_set(hh_kk_rts1),K,fs,fc);
                    a_vec_BS_Comp_real_set(:,hh_kk_rts1,hh_ll_rts1) = a_vec_BS_Comp_lk;
                    a_vec_AC_init_lk = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init(hh_ll_rts1),Ele_AC_set_init(hh_ll_rts1),...
                        K_index_set(hh_kk_rts1),K,fs,fc);
                    a_vec_AC_init_set_TrSq(:,hh_kk_rts1,hh_ll_rts1) = a_vec_AC_init_lk;
                    a_vec_AC_Comp_lk = A0_Array_Response_Vector_Comp_AC(hh_ll_rts1,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp(:,:,hh_ll_rts1),...
                        N_AC_h,N_AC_v,M_AC_h,M_AC_v,N_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp,M_AC_Comp_h,...
                        Azi_AC_set_init(hh_ll_rts1),Ele_AC_set_init(hh_ll_rts1),K_index_set(hh_kk_rts1),K,fs,fc);
                    a_vec_AC_Comp_real_set(:,hh_kk_rts1,hh_ll_rts1) = a_vec_AC_Comp_lk;
                    Equ_a_BS_lk_real_ase = (a_vec_BS_init_lk.*a_vec_BS_Comp_lk)'*f_RF_BS_init_real_set(:,hh_ll_rts1);
                    Equ_a_AC_lk_real_ase = w_RF_AC_init_real_set(:,hh_ll_rts1)'*(a_vec_AC_init_lk.*a_vec_AC_Comp_lk);
                    h_vec_init_real_ll_temp(hh_kk_rts1) = Equ_a_AC_lk_real_ase*Equ_a_BS_lk_real_ase;
                end
                h_vec_init_real_1(:,hh_ll_rts1) = alpha_init(hh_ll_rts1)*exp(-1i*2*pi*tau_init(hh_ll_rts1)*(-1/2+(K_index_set-1)/K)*fs).*...
                    h_vec_init_real_ll_temp;
            end
            h_vec_init_real_interf_1 = zeros(K,L);
            for hh_ll_rts2 = 1:L
                L_set_rts2 = 1:L; L_set_rts2(hh_ll_rts2) = []; ll_set_rts2 = L_set_rts2;
                h_vec_init_real_interf_ll_temp = zeros(K,1);
                for hh_kk_rts2 = 1:K
                    Equ_a_BS_real_1_ll_interf = (a_vec_BS_init_set_TrSq(:,hh_kk_rts2,ll_set_rts2).*a_vec_BS_Comp_real_set(:,hh_kk_rts2,ll_set_rts2))'*...
                        f_RF_BS_init_real_set(:,ll_set_rts2);
                    Equ_a_AC_real_1_ll_interf = w_RF_AC_init_real_set(:,hh_ll_rts2)'*(a_vec_AC_init_set_TrSq(:,hh_kk_rts2,ll_set_rts2).*...
                        a_vec_AC_Comp_real_set(:,hh_kk_rts2,hh_ll_rts2));
                    h_vec_init_real_interf_ll_temp(hh_kk_rts2) = Equ_a_AC_real_1_ll_interf*Equ_a_BS_real_1_ll_interf;
                end
                h_vec_init_real_interf_1(:,hh_ll_rts2) = alpha_init(ll_set_rts2)*exp(-1i*2*pi*tau_init(ll_set_rts2)*(-1/2+(K_index_set-1)/K)*fs).*...
                    h_vec_init_real_interf_ll_temp;
            end
        end
        if iter == 1
            for ase_ll_r1 = 1:L
                ASE_Per_CSI_TrSq_temp_ll = zeros(K,1);
                for ase_lk_r1 = 1:K
                    ASE_Per_CSI_TrSq_temp_ll(ase_lk_r1) = log2(1 + (abs(h_vec_init_real_1(ase_lk_r1,ase_ll_r1))^2/...
                        (abs(h_vec_init_real_interf_1(ase_lk_r1,ase_ll_r1))^2 + sigma_no^2*norm(w_RF_AC_init_real_set(:,ase_ll_r1))^2)));
                end
                ASE_Per_CSI_TrSq(ii) = ASE_Per_CSI_TrSq(ii) + mean(ASE_Per_CSI_TrSq_temp_ll);
            end
        end
        Doppler_k_set_init = zeros(K,L); Doppler_k_set_est = zeros(K,L);
        radial_vt_BS_set_est = lambda_z*Est_Doppler_set_1;
        for ll_ase_kk = 1:L
            Doppler_k_set_init(:,ll_ase_kk) = Doppler_z_init(ll_ase_kk) + radial_vt_BS_set_init(ll_ase_kk)*((K_index_set-1)/K-0.5)*fs/3e8;
            Doppler_k_set_est(:,ll_ase_kk) = Est_Doppler_set_1(ll_ase_kk) + radial_vt_BS_set_est(ll_ase_kk)*((K_index_set-1)/K-0.5)*fs/3e8;
        end
        h_vec_init_est_1 = zeros(K,L); a_vec_BS_Comp_est_set = zeros(N_BS,K,L); a_vec_AC_Comp_est_set = zeros(N_AC,K,L);
        for hh_ll_ets1 = 1:L
            h_vec_init_est_ll_temp = zeros(K,1);
            for hh_kk_ets1 = 1:K
                a_vec_BS_Comp_lk = A0_Array_Response_Vector_Comp_BS(Index_BS_Ant_Comp,N_BS_h,N_BS_v,N_BS_Comp,N_BS_Comp_h,N_BS_Comp_v,...
                    M_BS_Comp,M_BS_Comp_h,Est_Azi_BS_set_1(hh_ll_ets1),Est_Ele_BS_set_1(hh_ll_ets1),K_index_set(hh_kk_ets1),K,fs,fc);
                a_vec_BS_Comp_est_set(:,hh_kk_ets1,hh_ll_ets1) = a_vec_BS_Comp_lk;
                a_vec_AC_Comp_lk = A0_Array_Response_Vector_Comp_AC(hh_ll_ets1,I_AC_h,Index_W_RF_Sub_Array,Index_AC_Ant_Comp(:,:,hh_ll_ets1),...
                    N_AC_h,N_AC_v,M_AC_h,M_AC_v,N_AC_Comp,N_AC_Comp_h,N_AC_Comp_v,M_AC_Comp,M_AC_Comp_h,...
                    Est_Azi_AC_set_1(hh_ll_ets1),Est_Ele_AC_set_1(hh_ll_ets1),K_index_set(hh_kk_ets1),K,fs,fc);
                a_vec_AC_Comp_est_set(:,hh_kk_ets1,hh_ll_ets1) = a_vec_AC_Comp_lk;
                Equ_a_BS_lk_est_ase = (a_vec_BS_init_set_TrSq(:,hh_kk_ets1,hh_ll_ets1).*a_vec_BS_Comp_lk)'*f_RF_BS_init_set_BeSq(:,hh_ll_ets1);
                Equ_a_AC_lk_est_ase = w_RF_AC_init_set_BeSq(:,hh_ll_ets1)'*(a_vec_AC_init_set_TrSq(:,hh_kk_ets1,hh_ll_ets1).*a_vec_AC_Comp_lk);
                h_vec_init_est_ll_temp(hh_kk_ets1) = Equ_a_AC_lk_est_ase*Equ_a_BS_lk_est_ase;
            end
            h_vec_init_est_1(:,hh_ll_ets1) = alpha_init(hh_ll_ets1)*exp(-1i*2*pi*tau_init(hh_ll_ets1)*(-1/2+(K_index_set-1)/K)*fs).*...
                exp(1i*2*pi*(Doppler_k_set_init(:,hh_ll_ets1)-Doppler_k_set_est(:,hh_ll_ets1))*T_sym).*h_vec_init_est_ll_temp;
        end
        h_vec_init_est_interf_1 = zeros(K,L);
        for hh_ll_ets2 = 1:L
            L_set_ets2 = 1:L; L_set_ets2(hh_ll_ets2) = []; ll_set_ets2 = L_set_ets2;
            h_vec_init_est_interf_ll_temp = zeros(K,1);
            for hh_kk_ets2 = 1:K
                Equ_a_BS_est_1_ll_interf = (a_vec_BS_init_set_TrSq(:,hh_kk_ets2,ll_set_ets2).*a_vec_BS_Comp_est_set(:,hh_kk_ets2,ll_set_ets2))'*...
                    f_RF_BS_init_set_BeSq(:,ll_set_ets2);
                Equ_a_AC_est_1_ll_interf = w_RF_AC_init_set_BeSq(:,hh_ll_ets2)'*(a_vec_AC_init_set_TrSq(:,hh_kk_ets2,ll_set_ets2).*...
                    a_vec_AC_Comp_est_set(:,hh_kk_ets2,hh_ll_ets2));
                h_vec_init_est_interf_ll_temp(hh_kk_ets2) = Equ_a_AC_est_1_ll_interf*Equ_a_BS_est_1_ll_interf;
            end
            h_vec_init_est_interf_1(:,hh_ll_ets2) = alpha_init(ll_set_ets2)*exp(-1i*2*pi*tau_init(ll_set_ets2)*(-1/2+(K_index_set-1)/K)*fs).*...
                exp(1i*2*pi*(Doppler_k_set_init(:,ll_set_ets2)-Doppler_k_set_est(:,ll_set_ets2))*T_sym).*h_vec_init_est_interf_ll_temp;
        end
        for ase_ll_e1 = 1:L
            ASE_Est_CSI_TrSq_temp_ll = zeros(K,1);
            for ase_lk_e1 = 1:K
                ASE_Est_CSI_TrSq_temp_ll(ase_lk_e1) = log2(1 + (abs(h_vec_init_est_1(ase_lk_e1,ase_ll_e1))^2/...
                    (abs(h_vec_init_est_interf_1(ase_lk_e1,ase_ll_e1))^2 + sigma_no^2*norm(w_RF_AC_init_set_BeSq(:,ase_ll_e1))^2)));
            end
            ASE_Est_CSI_TrSq_temp(ii) = ASE_Est_CSI_TrSq_temp(ii) + mean(ASE_Est_CSI_TrSq_temp_ll);
        end
        
        %%
        if mod(iter,100) == 0
            toc
            disp(['  SNR_dB = ' num2str(SNR_dBs(ii)) ', iter = ' num2str(iter)])
        end
    end
    disp(['  SNR_dB = ' num2str(SNR_dBs(ii))])
end
ASE_Est_CSI_No_TrSq = ASE_Est_CSI_No_TrSq_temp/iterMax;
ASE_Est_CSI_TrSq = ASE_Est_CSI_TrSq_temp/iterMax;
disp('Finish  All')

MarkerSize = 6;
LineWidth = 1.2;
Fontsize = 15;

figure
plot(SNR_dBs,ASE_Per_CSI_No_TrSq,'-bd','LineWidth',LineWidth,'MarkerSize',MarkerSize); hold on; grid on;
plot(SNR_dBs,ASE_Per_CSI_TrSq,'-r^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
plot(SNR_dBs,ASE_Est_CSI_No_TrSq,'-mp','LineWidth',LineWidth,'MarkerSize',MarkerSize);
plot(SNR_dBs,ASE_Est_CSI_TrSq,'-ko','LineWidth',LineWidth,'MarkerSize',MarkerSize);
xlabel('Number of CCT','Fontsize',Fontsize),ylabel('NMSE','Fontsize',Fontsize);
title('NMSE vs CCT at -20 dB','Fontsize',Fontsize);
set(gca, 'GridLineStyle', '-.','FontSize',Fontsize, 'linewidth',1.5,'Fontname','Times New Roman');
set(gcf, 'position', [700 300 700 500]); axis normal;
h2 = legend('Perfect CSI + no triple squint','Perfect CSI + triple squint',...
    'Estimated CSI + no triple squint','Estimated CSI + triple squint','Location','northwest');
set(h2,'Fontsize',11);
