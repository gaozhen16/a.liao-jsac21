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
N_BS_h = 200;
N_BS_v = N_BS_h;
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
G_ls_init = G_AC*G_BS*lambda^2./(4*pi*Distance_init).^2;
unit_d = unit_d_set(:,rand_num);
vt_vec = vt*unit_d;
radial_vt_BS_set_init = Calculate_Radial_Velocity(vt_vec, Coo_Air_C_init, Coo_BS1_A, Coo_BS2_B, L);
Doppler_init = radial_vt_BS_set_init./lambda;
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
N_OFDM_Delay = 10;
N_OFDM_Doppler = 6;
N_OFDM = N_OFDM_Angle_Est + N_OFDM_Delay + N_OFDM_Doppler;
T_CCT = N_OFDM*T_sym;
niu_Doppler_init = 2*pi*Doppler_init*T_sym;
m_CCT_init = 1;
Angle_diff = deg2rad(10);
miu_BS_set_init_error = pi*sin(Azi_BS_set_init_error).*cos(Ele_BS_set_init_error);
niu_BS_set_init_error = pi*sin(Ele_BS_set_init_error);
miu_AC_set_init_error = pi*sin(Azi_AC_set_init_error).*cos(Ele_AC_set_init_error);
niu_AC_set_init_error = pi*sin(Ele_AC_set_init_error);
Jacobian_Mat_Azi_Ele_BS = Gene_Jacobian_Mat_Azi_Ele(miu_BS_set_init, niu_BS_set_init, Delta_BS, lambda, d_ant, L);
Jacobian_Mat_Azi_Ele_AC = Gene_Jacobian_Mat_Azi_Ele(miu_AC_set_init, niu_AC_set_init, Delta_AC, lambda, d_ant, L);
Jacobian_Mat_tau = Gene_Jacobian_Mat_tau(miu_tau_set_init, K, 1, D_subc, L);

%% set snr
iterMax = 1e3;
Power_Tx_dBm = 80:5:150;
Ptx_gain = ones(1,length(Power_Tx_dBm));
SNR_dBs = -50:10:30;
sigma2_NSD = -174;
sigma2_no = 1e-3*10^(sigma2_NSD/10)*BW;
sigma_no = sqrt(sigma2_no);
MSE_Azi_BS_set_temp = zeros(length(SNR_dBs),1);
MSE_Ele_BS_set_temp = zeros(length(SNR_dBs),1);
CRB_Azi_BS_set = zeros(length(SNR_dBs),1);
CRB_Ele_BS_set = zeros(length(SNR_dBs),1);

%%
H_Beam_Alignment_BS_angle = zeros(N_OFDM_BS_angle,L);
a_vec_AC_init_set = zeros(N_AC,L); a_vec_BS_init_set = zeros(N_BS,L);
f_RF_AC_init_BS_set = zeros(N_AC,L); w_RF_BS_init_BS_1_set = zeros(N_BS,L);
for ll_bs_ang1 = 1:L
    a_vec_AC_init_ll = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init(ll_bs_ang1),Ele_AC_set_init(ll_bs_ang1));
    a_vec_BS_init_ll = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init(ll_bs_ang1),Ele_BS_set_init(ll_bs_ang1));
    a_vec_AC_init_set(:,ll_bs_ang1) = a_vec_AC_init_ll; a_vec_BS_init_set(:,ll_bs_ang1) = a_vec_BS_init_ll;
    a_vec_AC_init_ll_error = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init_error(ll_bs_ang1),Ele_AC_set_init_error(ll_bs_ang1));
    a_vec_BS_init_ll_error = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init_error(ll_bs_ang1),Ele_BS_set_init_error(ll_bs_ang1));
    f_RF_AC_init_ll = zeros(N_AC,1);
    f_RF_AC_init_ll(Index_W_RF_Sub_Array(:,ll_bs_ang1)) = Quantize(a_vec_AC_init_ll_error(Index_W_RF_Sub_Array(:,ll_bs_ang1)),N_bits)/sqrt(M_AC);
    f_RF_AC_init_BS_set(:,ll_bs_ang1) = f_RF_AC_init_ll;
    Equ_a_AC_ll_init = a_vec_AC_init_ll'*f_RF_AC_init_ll;
    Subarray_BS_ll_init = Quantize(a_vec_BS_init_ll_error(Index_W_RF_BS_Angle(:,1)),N_bits)/sqrt(N_BS_bar); % 
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
%% Start
for ii = 1:length(SNR_dBs)
    Power_Tx_ii = (1e-3*10^(Power_Tx_dBm(ii)/10));
    Ptx_gain_ii = Ptx_gain(ii);
    sigma_no = sqrt(10^(-(SNR_dBs(ii)/10)));
    tic
    for iter = 1:iterMax
        [Est_miu_BS_set,Est_niu_BS_set,Est_Azi_BS_set,Est_Ele_BS_set] = A1_Est_BS_angle_V1(H_Beam_Alignment_BS_angle,G_ls_init,alpha_init,...
            k_index_com,Power_Tx_ii,sigma_no,awgn_en,K_miu_BS_Re, K_miu_BS_Im, K_niu_BS_Re, K_niu_BS_Im, ...
            Ptx_gain_ii,Delta_BS,miu_BS_set_init_error,niu_BS_set_init_error,s_BS_angle_set);
        MSE_Azi_BS_set_temp(ii,1) = MSE_Azi_BS_set_temp(ii,1) + norm(Est_Azi_BS_set - Azi_BS_set_init)^2/L;
        MSE_Ele_BS_set_temp(ii,1) = MSE_Ele_BS_set_temp(ii,1) + norm(Est_Ele_BS_set - Ele_BS_set_init)^2/L;
        
        if ii == 1 && iter == 1
            Mat_ZPZ_2D_BS = zeros(2,2,L);
            for ll_crb_bs1 = 1:L
                i_BS_h_bar = (0:I_BS_h_bar-1).'; a_BS_miu_eff_ll = exp(1i*i_BS_h_bar*miu_BS_set_init(ll_crb_bs1));
                i_BS_v_bar = (0:I_BS_v_bar-1).'; a_BS_niu_eff_ll = exp(1i*i_BS_v_bar*niu_BS_set_init(ll_crb_bs1));
                Diff_a_BS_miu = Diff_Para_Azi_Ele(I_BS_h_bar, miu_BS_set_init(ll_crb_bs1));
                Com_Diff_a_BS_miu = Khatri_Rao(a_BS_niu_eff_ll,Diff_a_BS_miu);
                Diff_a_AC_niu = Diff_Para_Azi_Ele(I_BS_v_bar, niu_BS_set_init(ll_crb_bs1));
                Com_Diff_a_BS_niu = Khatri_Rao(Diff_a_AC_niu,a_BS_miu_eff_ll);
                Mat_Z_2D_BS_ll = [Com_Diff_a_BS_miu, Com_Diff_a_BS_niu];
                a_BS_bar_2D = Khatri_Rao(a_BS_niu_eff_ll,a_BS_miu_eff_ll);
                Mat_P_A_orth_2D_BS_ll = eye(size(a_BS_bar_2D,1)) - a_BS_bar_2D*inv(a_BS_bar_2D'*a_BS_bar_2D)*a_BS_bar_2D';
                Mat_ZPZ_2D_BS(:,:,ll_crb_bs1) = Mat_Z_2D_BS_ll'*Mat_P_A_orth_2D_BS_ll*Mat_Z_2D_BS_ll;
            end
        end
        if iter == 1
            CRB_Azi_BS_temp = zeros(L,1); CRB_Ele_BS_temp = zeros(L,1);
            for ll_crb_bs2 = 1:L
                Equ_a_AC_bs_crb_ll = a_vec_AC_init_set(:,ll_crb_bs2)'*f_RF_AC_init_BS_set(:,ll_crb_bs2);
                Equ_a_BS_bs_crb_ll = w_RF_BS_init_BS_1_set(:,ll_crb_bs2)'*a_vec_BS_init_set(:,ll_crb_bs2);
                Tx_eff_BS_set_ll = alpha_init(ll_crb_bs2)*Equ_a_BS_bs_crb_ll*Equ_a_AC_bs_crb_ll*s_BS_angle_set(:,ll_crb_bs2);
                Mat_2D_BS = zeros(2,2);
                for kk_bs = 1:K_sub
                    Mat_B_2D_BS_nn = kron(eye(2),diag(Tx_eff_BS_set_ll(kk_bs)));
                    Mat_2D_BS = Mat_2D_BS + Mat_B_2D_BS_nn'*Mat_ZPZ_2D_BS(:,:,ll_crb_bs2)*Mat_B_2D_BS_nn;
                end
                CRB_miu_niu_Mtx_BS_ll = sigma_no^2*inv(real(Mat_2D_BS))/2;
                CRB_Azi_Ele_Mtx_BS_ll = Jacobian_Mat_Azi_Ele_BS(:,:,ll_crb_bs2)*CRB_miu_niu_Mtx_BS_ll*Jacobian_Mat_Azi_Ele_BS(:,:,ll_crb_bs2).';
                CRB_Azi_BS_temp(ll_crb_bs2) = CRB_Azi_Ele_Mtx_BS_ll(2,2);
                CRB_Ele_BS_temp(ll_crb_bs2) = CRB_Azi_Ele_Mtx_BS_ll(1,1);
            end
            CRB_Azi_BS_set(ii) = sum(CRB_Azi_BS_temp)/L;
            CRB_Ele_BS_set(ii) = sum(CRB_Ele_BS_temp)/L;
        end
        
        if mod(iter,100) == 0
            disp(['    SNR_dB = ' num2str(SNR_dBs(ii)) ',   iter = ' num2str(iter)])
        end
    end
    disp(['  Finish SNR_dB = ' num2str(SNR_dBs(ii))])
end
MSE_Azi_BS_set = MSE_Azi_BS_set_temp/iterMax;
MSE_Ele_BS_set = MSE_Ele_BS_set_temp/iterMax;
disp('Finish  All')

MarkerSize = 3;
LineWidth = 1.2;
Fontsize = 15;
figure
semilogy(SNR_dBs,MSE_Azi_BS_set(:,1),':ro','LineWidth',LineWidth); hold on; grid on;
semilogy(SNR_dBs,CRB_Azi_BS_set,'--k','LineWidth',LineWidth+0.5);
semilogy(SNR_dBs,MSE_Ele_BS_set(:,1),'-ro','LineWidth',LineWidth);
semilogy(SNR_dBs,CRB_Ele_BS_set,'-k','LineWidth',LineWidth+0.5);
xlabel('SNR [dB]','Fontsize',Fontsize),ylabel('MSE','Fontsize',Fontsize);
title('BS Azi & Ele vs CRB','Fontsize',Fontsize);
set(gca, 'GridLineStyle', '-.','FontSize',Fontsize, 'linewidth',1.5,'Fontname','Times New Roman');
set(gcf, 'position', [700 300 650 550]); axis normal;
h1 = legend('\theta_{BS} + no BD','CRB \theta_{BS}','\phi_{BS} + no BD','CRB \phi_{BS}','Location','northeast');
set(h1,'Fontsize',11);

RMSE_Azi_BS_set = sqrt(MSE_Azi_BS_set); RCRB_Azi_BS_set = sqrt(CRB_Azi_BS_set);
RMSE_Ele_BS_set = sqrt(MSE_Ele_BS_set); RCRB_Ele_BS_set = sqrt(CRB_Ele_BS_set);

figure
semilogy(SNR_dBs,RMSE_Azi_BS_set(:,1),':ro','LineWidth',LineWidth); hold on; grid on;
semilogy(SNR_dBs,RCRB_Azi_BS_set,'--k','LineWidth',LineWidth+0.5);
semilogy(SNR_dBs,RMSE_Ele_BS_set(:,1),'-ro','LineWidth',LineWidth);
semilogy(SNR_dBs,RCRB_Ele_BS_set,'-k','LineWidth',LineWidth+0.5);
xlabel('SNR [dB]','Fontsize',Fontsize),ylabel('RMSE','Fontsize',Fontsize);
title('BS Azi & Ele vs CRB','Fontsize',Fontsize);
set(gca, 'GridLineStyle', '-.','FontSize',Fontsize, 'linewidth',1.5,'Fontname','Times New Roman');
set(gcf, 'position', [700 300 650 550]); axis normal;
h2 = legend('\theta_{BS} + no BD','CRB \theta_{BS}','\phi_{BS} + no BD','CRB \phi_{BS}','Location','northeast');
set(h2,'Fontsize',11);
