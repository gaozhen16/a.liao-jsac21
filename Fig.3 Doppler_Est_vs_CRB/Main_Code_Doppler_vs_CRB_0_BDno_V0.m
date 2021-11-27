clc,clear

%% System parameters
fc = 100e9;
lambda = 3e8/fc;
d_ant = lambda/2;
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
G_ls_init = G_AC*G_BS*lambda^2./(4*pi*Distance_init).^2;
unit_d = unit_d_set(:,rand_num);
vt_vec = vt*unit_d;
radial_vt_BS_set_init = Calculate_Radial_Velocity(vt_vec, Coo_Air_C_init, Coo_BS1_A, Coo_BS2_B, L);
Doppler_z_init = radial_vt_BS_set_init./lambda;
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
Doppler_diff = Doppler_z_init*1e-3;
n_ofdm_delay_doppler = 1:N_OFDM_Delay+N_OFDM_Doppler;
Phase_Doppler_Del_Dop_set = 2*pi*Doppler_z_init*T_sym*((n_ofdm_delay_doppler-1)+(m_CCT_init-1)*N_OFDM);
Phase_Doppler_Dop_set = Phase_Doppler_Del_Dop_set(:,1:N_OFDM_Doppler);
Angle_diff = deg2rad(5);
miu_BS_set_init_error = pi*sin(Azi_BS_set_init_error).*cos(Ele_BS_set_init_error);
niu_BS_set_init_error = pi*sin(Ele_BS_set_init_error);
miu_AC_set_init_error = pi*sin(Azi_AC_set_init_error).*cos(Ele_AC_set_init_error);
niu_AC_set_init_error = pi*sin(Ele_AC_set_init_error);
Jacobian_Mat_Azi_Ele_BS = Gene_Jacobian_Mat_Azi_Ele(miu_BS_set_init, niu_BS_set_init, Delta_BS, lambda, d_ant, L);
Jacobian_Mat_Azi_Ele_AC = Gene_Jacobian_Mat_Azi_Ele(miu_AC_set_init, niu_AC_set_init, Delta_AC, lambda, d_ant, L);
Jacobian_Mat_Doppler = Gene_Jacobian_Mat_Doppler(niu_Doppler_init, T_sym, L);
Jacobian_Mat_tau = Gene_Jacobian_Mat_tau(miu_tau_set_init, K, 1, D_subc, L);

%% set snr
iterMax = 1e3;
Power_Tx_dBm = 80:5:150;
Ptx_gain = ones(1,length(Power_Tx_dBm));
SNR_dBs = -120:10:-30;
sigma2_NSD = -174;
sigma2_no = 1e-3*10^(sigma2_NSD/10)*BW;
sigma_no = sqrt(sigma2_no);
MSE_Doppler_set_temp_0 = zeros(length(SNR_dBs),1);
CRB_Doppler_set = zeros(length(SNR_dBs),1);
MSE_Doppler_set_temp_1 = zeros(length(SNR_dBs),1);
MSE_Doppler_set_temp_2 = zeros(length(SNR_dBs),1);
MSE_Doppler_set_temp_3 = zeros(length(SNR_dBs),1);

%% 
H_Beam_Alignment_BS_angle = zeros(N_OFDM_BS_angle,L);
a_vec_AC_init_set = zeros(N_AC,L); a_vec_BS_init_set = zeros(N_BS,L);
f_RF_AC_init_BS_set = zeros(N_AC,L); w_RF_BS_init_BS_1_set = zeros(N_BS,L);
for ll_bs_ang1 = 1:L
    a_vec_AC_init_ll = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init(ll_bs_ang1),Ele_AC_set_init(ll_bs_ang1));
    a_vec_BS_init_ll = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init(ll_bs_ang1),Ele_BS_set_init(ll_bs_ang1));
    a_vec_AC_init_set(:,ll_bs_ang1) = a_vec_AC_init_ll; a_vec_BS_init_set(:,ll_bs_ang1) = a_vec_BS_init_ll;
    a_vec_AC_init_ll_est = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init_error(ll_bs_ang1),Ele_AC_set_init_error(ll_bs_ang1));
    a_vec_BS_init_ll_error = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init_error(ll_bs_ang1),Ele_BS_set_init_error(ll_bs_ang1));
    f_RF_AC_init_ll = zeros(N_AC,1);
    f_RF_AC_init_ll(Index_W_RF_Sub_Array(:,ll_bs_ang1)) = Quantize(a_vec_AC_init_ll_est(Index_W_RF_Sub_Array(:,ll_bs_ang1)),N_bits)/sqrt(M_AC);
    f_RF_AC_init_BS_set(:,ll_bs_ang1) = f_RF_AC_init_ll;
    Equ_a_AC_ll_init = a_vec_AC_init_ll'*f_RF_AC_init_ll;
    Subarray_BS_ll_init = Quantize(a_vec_BS_init_ll_error(Index_W_RF_BS_Angle(:,1)),N_bits)/sqrt(N_BS_bar);
    for nn_bs1 = 1:N_OFDM_BS_angle
        w_RF_BS_init_ll_nn = zeros(N_BS,1);
        w_RF_BS_init_ll_nn(Index_W_RF_BS_Angle(:,nn_bs1)) = Subarray_BS_ll_init;
        if nn_bs1 == 1
            w_RF_BS_init_BS_1_set(:,ll_bs_ang1) = w_RF_BS_init_ll_nn;
        end
        Equ_a_BS_ll_init = w_RF_BS_init_ll_nn'*a_vec_BS_init_ll;
        H_Beam_Alignment_BS_angle(nn_bs1,ll_bs_ang1) = Equ_a_BS_ll_init*Equ_a_AC_ll_init;
    end
end
sigma_no = sqrt(10^(-(-20/10)));
Power_Tx_ii = (1e-3*10^(Power_Tx_dBm(1)/10)); Ptx_gain_ii = Ptx_gain(1);
[Est_miu_BS_set,Est_niu_BS_set,Est_Azi_BS_set,Est_Ele_BS_set] = A1_Est_BS_angle_V1(H_Beam_Alignment_BS_angle,G_ls_init,alpha_init,...
    k_index_com,Power_Tx_ii,sigma_no,awgn_en,K_miu_BS_Re, K_miu_BS_Im, K_niu_BS_Re, K_niu_BS_Im, ...
    Ptx_gain_ii,Delta_BS,miu_BS_set_init_error,niu_BS_set_init_error,s_BS_angle_set);

H_Beam_Alignment_AC_angle = zeros(N_OFDM_AC_angle,L);
f_RF_BS_init_AC_set = zeros(N_BS,L); w_RF_AC_init_AC_1_set = zeros(N_AC,L);
for ll_ac_ang1 = 1:L
    a_vec_BS_init_ll_est = A0_Array_Response_Vector(N_BS_h,N_BS_v,Est_Azi_BS_set(ll_ac_ang1),Est_Ele_BS_set(ll_ac_ang1));
    a_vec_AC_init_ll_error = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init_error(ll_ac_ang1),Ele_AC_set_init_error(ll_ac_ang1));
    f_RF_BS_init_ll = Quantize(a_vec_BS_init_ll_est,N_bits)/sqrt(N_BS);
    f_RF_BS_init_AC_set(:,ll_ac_ang1) = f_RF_BS_init_ll;
    Equ_a_BS_ll_init = a_vec_BS_init_set(:,ll_ac_ang1)'*f_RF_BS_init_ll;
    Subarray_AC_ll_init = Quantize(a_vec_AC_init_ll_error(Index_W_RF_AC_Angle(:,1,ll_ac_ang1)),N_bits)/sqrt(M_AC_bar);
    for nn_ac1 = 1:N_OFDM_AC_angle
        w_RF_AC_init_ll_nn = zeros(N_AC,1);	% w_BB_BS_n = 1;
        w_RF_AC_init_ll_nn(Index_W_RF_AC_Angle(:,nn_ac1,ll_ac_ang1)) = Subarray_AC_ll_init;
        if nn_ac1 == 1
            w_RF_AC_init_AC_1_set(:,ll_ac_ang1) = w_RF_AC_init_ll_nn;
        end
        Equ_a_AC_ll_init = w_RF_AC_init_ll_nn'*a_vec_AC_init_set(:,ll_ac_ang1);
        H_Beam_Alignment_AC_angle(nn_ac1,ll_ac_ang1) = Equ_a_AC_ll_init*Equ_a_BS_ll_init;
    end
end
[Est_miu_AC_set,Est_niu_AC_set,Est_Azi_AC_set,Est_Ele_AC_set] = A2_Est_AC_angle_V1(H_Beam_Alignment_AC_angle,G_ls_init,alpha_init,...
    tau_init,fs,k_index_com,Power_Tx_ii,sigma_no,awgn_en,K_miu_AC_Re,K_miu_AC_Im,K_niu_AC_Re,K_niu_AC_Im,...
    Ptx_gain_ii,Delta_AC,miu_AC_set_init_error,niu_AC_set_init_error,s_AC_angle_set);

H_Beam_Alignment_Doppler = zeros(1,L);
f_RF_BS_init_Dop_set = zeros(N_BS,L); w_RF_AC_init_Dop_set = zeros(N_AC,L);
for ll_dop1 = 1:L
    a_vec_BS_init_ll_est = A0_Array_Response_Vector(N_BS_h,N_BS_v,Est_Azi_BS_set(ll_dop1),Est_Ele_BS_set(ll_dop1));
    a_vec_AC_init_ll_est = A0_Array_Response_Vector(N_AC_h,N_AC_v,Est_Azi_AC_set(ll_dop1),Est_Ele_AC_set(ll_dop1));
    f_RF_BS_init_ll = Quantize(a_vec_BS_init_ll_est,N_bits)/sqrt(N_BS);
    f_RF_BS_init_Dop_set(:,ll_dop1) = f_RF_BS_init_ll;
    Equ_a_BS_ll_init = a_vec_BS_init_set(:,ll_dop1)'*f_RF_BS_init_ll;
    w_RF_AC_init_ll = zeros(N_AC,1);
    w_RF_AC_init_ll(Index_W_RF_Sub_Array(:,ll_dop1)) = Quantize(a_vec_AC_init_ll_est(Index_W_RF_Sub_Array(:,ll_dop1)),N_bits)/sqrt(M_AC);
    w_RF_AC_init_Dop_set(:,ll_dop1) = w_RF_AC_init_ll;
    Equ_a_AC_ll_init = w_RF_AC_init_ll'*a_vec_AC_init_set(:,ll_dop1);
    H_Beam_Alignment_Doppler(:,ll_dop1) = Equ_a_AC_ll_init*Equ_a_BS_ll_init;
end

%% Start
for ii = 1:length(SNR_dBs)
    Power_Tx_ii = (1e-3*10^(Power_Tx_dBm(ii)/10));
    Ptx_gain_ii = Ptx_gain(ii);
    sigma_no = sqrt(10^(-(SNR_dBs(ii)/10)));
    tic
    for iter = 1:iterMax
        [Est_niu_Doppler,Est_Doppler_set,Y_dop_bar_set] = A3_Est_Dopple_V1(H_Beam_Alignment_Doppler,G_ls_init,alpha_init,Phase_Doppler_Dop_set,...
            tau_init,fs,k_index_com,N_OFDM_Doppler,Power_Tx_ii,sigma_no,awgn_en,Ptx_gain_ii,T_sym,s_dop_set);
        MSE_Doppler_set_temp_0(ii) = MSE_Doppler_set_temp_0(ii) + norm(Est_Doppler_set - Doppler_z_init)^2/L;
        
        if ii == 1 && iter == 1
            nn_dop = (1:N_OFDM_Doppler).';
            Mat_ZPZ_Dop = zeros(1,1,L);
            for ll_crb_dop1 = 1:L
                Mat_Z_Dop_ll = Diff_Para_Doppler(niu_Doppler_init(ll_crb_dop1), nn_dop);
                a_bar_Dop = exp(1i*(nn_dop-1).*niu_Doppler_init(ll_crb_dop1));
                Mat_P_A_orth_Dop = eye(size(a_bar_Dop,1)) - a_bar_Dop*inv(a_bar_Dop'*a_bar_Dop)*a_bar_Dop';
                Mat_ZPZ_Dop(:,:,ll_crb_dop1) = Mat_Z_Dop_ll'*Mat_P_A_orth_Dop*Mat_Z_Dop_ll;
            end
        end
        if iter == 1
            CRB_niu_Doppler_set_temp = zeros(L,1); CRB_Doppler_set_temp = zeros(L,1);
            for ll_crb_dop2 = 1:L
                Tx_eff_Dop_set_ll = alpha_init(ll_crb_dop2)*H_Beam_Alignment_Doppler(ll_crb_dop2)*...
                    exp(-1i*2*pi*tau_init(ll_crb_dop2)*(-1/2+(k_index_com(:,ll_crb_dop2)-1)/K)*fs).*s_dop_set(:,ll_crb_dop2); % 
                Mat_1D_dop_ll = 0;
                for kk_Dop = 1:K_sub
                    Mat_1D_dop_ll = Mat_1D_dop_ll + Tx_eff_Dop_set_ll(kk_Dop)'*Mat_ZPZ_Dop(:,:,ll_crb_dop2)*Tx_eff_Dop_set_ll(kk_Dop);
                end
                CRB_niu_Doppler_set_temp(ll_crb_dop2) = sigma_no^2*inv(real(Mat_1D_dop_ll))/2;
                CRB_Doppler_set_temp(ll_crb_dop2) = Jacobian_Mat_Doppler(:,:,ll_crb_dop2)*CRB_niu_Doppler_set_temp(ll_crb_dop2)*Jacobian_Mat_Doppler(:,:,ll_crb_dop2).';
            end
            CRB_Doppler_set(ii) = sum(CRB_Doppler_set_temp)/L;
        end
        
        %%
        if mod(iter,100) == 0
            toc
            disp(['    SNR_dB = ' num2str(SNR_dBs(ii)) ',   iter = ' num2str(iter)])
        end
    end
    disp(['  Finish SNR_dB = ' num2str(SNR_dBs(ii))])
end
MSE_Doppler_set_0 = MSE_Doppler_set_temp_0/iterMax;
disp('Finish  All')

MarkerSize = 3;
LineWidth = 1.2;
Fontsize = 15;
figure
semilogy(SNR_dBs,MSE_Doppler_set_0,'-ro','LineWidth',LineWidth); hold on; grid on;
semilogy(SNR_dBs,CRB_Doppler_set,'-k','LineWidth',LineWidth+0.5);
xlabel('SNR [dB]','Fontsize',Fontsize),ylabel('MSE','Fontsize',Fontsize);
title('Doppler vs CRB','Fontsize',Fontsize);
set(gca, 'GridLineStyle', '-.','FontSize',Fontsize, 'linewidth',1.5,'Fontname','Times New Roman');
set(gcf, 'position', [700 300 650 550]); axis normal;
h1 = legend('\psi+no DSS','CRB \psi','Location','southwest');
set(h1,'Fontsize',11);

RCRB_Doppler_set = sqrt(CRB_Doppler_set);
RMSE_Doppler_set_0 = sqrt(MSE_Doppler_set_0);
figure
semilogy(SNR_dBs,RMSE_Doppler_set_0,'-ro','LineWidth',LineWidth); hold on; grid on;
semilogy(SNR_dBs,RCRB_Doppler_set,'-k','LineWidth',LineWidth+0.5);
axis([-120 -30 1e-1 2.5e5])
xlabel('SNR [dB]','Fontsize',Fontsize),ylabel('RMSE','Fontsize',Fontsize);
title('Doppler vs CRB','Fontsize',Fontsize);
set(gca, 'GridLineStyle', '-.','FontSize',Fontsize, 'linewidth',1.5,'Fontname','Times New Roman');
set(gcf, 'position', [700 300 650 550]); axis normal;
h2 = legend('\psi+no DSS','CRB \psi','Location','southwest');
set(h2,'Fontsize',11);
