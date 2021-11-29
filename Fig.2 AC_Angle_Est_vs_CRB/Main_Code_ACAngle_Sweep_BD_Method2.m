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
I_AC_h_bar = 5;
I_AC_v_bar = I_AC_h_bar;
I_AC_bar = I_AC_h_bar*I_AC_v_bar;
Index_W_RF_Sub_Array = A0_Gene_Index_W_RF_2Sub_Array(N_AC_h,N_AC_v,M_AC,M_AC_h,M_AC_v,I_AC,I_AC_h);
K_index_set = 1:K; K_sub = K/L;
D_subc = 2;
if mod(K,L) == 0
	k_index_com = reshape(K_index_set,L,K_sub).';
else
    error('Error: K should be the integral number of L !');
end

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
Angle_diff = deg2rad(5);
miu_BS_set_init_error = pi*sin(Azi_BS_set_init_error).*cos(Ele_BS_set_init_error);
niu_BS_set_init_error = pi*sin(Ele_BS_set_init_error);
miu_AC_set_init_error = pi*sin(Azi_AC_set_init_error).*cos(Ele_AC_set_init_error);
niu_AC_set_init_error = pi*sin(Ele_AC_set_init_error);
Angle_diff_sw = deg2rad(5);
Azi_BS_set_sw = zeros(I_BS_h_bar,L); Ele_BS_set_sw = zeros(I_BS_v_bar,L);
Azi_AC_set_sw = zeros(I_AC_h_bar,L); Ele_AC_set_sw = zeros(I_AC_v_bar,L);
for ll_sw = 1:L
    Azi_BS_set_sw(:,ll_sw) = Azi_BS_set_init_error(ll_sw)-Angle_diff_sw : 2*Angle_diff_sw/(I_BS_h_bar-1) : Azi_BS_set_init_error(ll_sw)+Angle_diff_sw;
    Ele_BS_set_sw(:,ll_sw) = Ele_BS_set_init_error(ll_sw)-Angle_diff_sw : 2*Angle_diff_sw/(I_BS_v_bar-1) : Ele_BS_set_init_error(ll_sw)+Angle_diff_sw;
    Azi_AC_set_sw(:,ll_sw) = Azi_AC_set_init_error(ll_sw)-Angle_diff_sw : 2*Angle_diff_sw/(I_AC_h_bar-1) : Azi_AC_set_init_error(ll_sw)+Angle_diff_sw;
    Ele_AC_set_sw(:,ll_sw) = Ele_AC_set_init_error(ll_sw)-Angle_diff_sw : 2*Angle_diff_sw/(I_AC_v_bar-1) : Ele_AC_set_init_error(ll_sw)+Angle_diff_sw;
end

%% 
H_BS_angle_Sw = zeros(K_sub,N_OFDM_BS_angle,L);
a_vec_AC_init_be_sq_set = zeros(N_AC,K_sub,L); a_vec_BS_init_be_sq_set = zeros(N_BS,K_sub,L);
for ll_bs_sw = 1:L
    a_vec_AC_init_ll_error = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init_error(ll_bs_sw),Ele_AC_set_init_error(ll_bs_sw));
    f_RF_AC_init_ll = zeros(N_AC,1);
    f_RF_AC_init_ll(Index_W_RF_Sub_Array(:,ll_bs_sw)) = Quantize(a_vec_AC_init_ll_error(Index_W_RF_Sub_Array(:,ll_bs_sw)),N_bits)/sqrt(M_AC); % 
    for ii_bs_h = 1:I_BS_h_bar
        for ii_bs_v = 1:I_BS_v_bar
            a_vec_BS_init_ll_sw = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_sw(ii_bs_h,ll_bs_sw),Ele_BS_set_sw(ii_bs_v,ll_bs_sw));
            w_RF_AC_init_ll_ii = Quantize(a_vec_BS_init_ll_sw,N_bits)/sqrt(N_BS); % 
            nn_bs_sw = (ii_bs_h-1)*I_BS_v_bar + ii_bs_v;
            for kk_BS_angle_ll = 1:K_sub
                if ii_bs_h == 1 && ii_bs_v == 1
                    a_vec_AC_init_lk = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_init(ll_bs_sw),Ele_AC_set_init(ll_bs_sw),...
                        k_index_com(kk_BS_angle_ll,ll_bs_sw),K,fs,fc);
                    a_vec_AC_init_be_sq_set(:,kk_BS_angle_ll,ll_bs_sw) = a_vec_AC_init_lk;
                    a_vec_BS_init_lk = A0_Array_Response_Vector(N_BS_h,N_BS_v,Azi_BS_set_init(ll_bs_sw),Ele_BS_set_init(ll_bs_sw),...
                        k_index_com(kk_BS_angle_ll,ll_bs_sw),K,fs,fc);
                    a_vec_BS_init_be_sq_set(:,kk_BS_angle_ll,ll_bs_sw) = a_vec_BS_init_lk;
                else
                    a_vec_AC_init_lk = a_vec_AC_init_be_sq_set(:,kk_BS_angle_ll,ll_bs_sw);
                    a_vec_BS_init_lk = a_vec_BS_init_be_sq_set(:,kk_BS_angle_ll,ll_bs_sw);
                end
                Equ_a_AC_lk_init = a_vec_AC_init_lk'*f_RF_AC_init_ll;
                Equ_a_BS_lk_init = w_RF_AC_init_ll_ii'*a_vec_BS_init_lk;
                H_BS_angle_Sw(kk_BS_angle_ll,nn_bs_sw,ll_bs_sw) = Equ_a_BS_lk_init*Equ_a_AC_lk_init;
            end
        end
    end
end

%% set snr
iterMax = 1e1;
SNR_dBs = -50:10:30;
MSE_Azi_AC_set_temp_1_Sw = zeros(length(SNR_dBs),1);
MSE_Ele_AC_set_temp_1_Sw = zeros(length(SNR_dBs),1);

%% Start
for ii = 1:length(SNR_dBs)
    sigma_no = sqrt(10^(-(SNR_dBs(ii)/10)));
    tic
    for iter = 1:iterMax
        Y_BS_ang_com = zeros(N_OFDM_BS_angle,K_sub,L);
        for ll_ang_bs = 1:L
            for kk_ang_ll = 1:K_sub
                for nn_BS_ll = 1:N_OFDM_BS_angle
                    Y_BS_ang_com(nn_BS_ll,kk_ang_ll,ll_ang_bs) = alpha_init(ll_ang_bs)*H_BS_angle_Sw(kk_ang_ll,nn_BS_ll,ll_ang_bs)*...
                        s_BS_angle_set(kk_ang_ll,ll_ang_bs) + awgn_en*sigma_no*(normrnd(0,1) + 1i*normrnd(0,1))/sqrt(2);
                end
            end
        end
        Est_Azi_BS_set = zeros(L,1); Est_Ele_BS_set = zeros(L,1);
        for ll_bs_pow = 1:L
            Y_BS_ang_ll = Y_BS_ang_com(:,:,ll_bs_pow);
            [maximum,index_max] = max(abs(sum(Y_BS_ang_ll,2).^2));
            ii_bs_h_hat = ceil(index_max/I_BS_v_bar);
            ii_bs_v_hat = index_max - (ii_bs_h_hat-1)*I_BS_v_bar;
            Est_Azi_BS_set(ll_bs_pow) = Azi_BS_set_sw(ii_bs_h_hat);
            Est_Ele_BS_set(ll_bs_pow) = Ele_BS_set_sw(ii_bs_v_hat);
        end
        
        H_AC_angle_Sw = zeros(K_sub,N_OFDM_AC_angle,L);
        for ll_ac_sw = 1:L
            a_vec_BS_init_ll_est = A0_Array_Response_Vector(N_BS_h,N_BS_v,Est_Azi_BS_set(ll_ac_sw),Est_Ele_BS_set(ll_ac_sw));
            f_RF_BS_init_ll = Quantize(a_vec_BS_init_ll_est,N_bits)/sqrt(N_BS); % f_BB_AC_n = 1;
            for ii_ac_h = 1:I_AC_h_bar
                for ii_ac_v = 1:I_AC_v_bar
                    a_vec_AC_init_ll_sw = A0_Array_Response_Vector(N_AC_h,N_AC_v,Azi_AC_set_sw(ii_ac_h,ll_ac_sw),Ele_AC_set_sw(ii_ac_v,ll_ac_sw));
                    w_RF_AC_init_ll_ii = zeros(N_AC,1);	% w_BB_BS_n = 1;
                    w_RF_AC_init_ll_ii(Index_W_RF_Sub_Array(:,ll_ac_sw)) = Quantize(a_vec_AC_init_ll_sw(Index_W_RF_Sub_Array(:,ll_ac_sw)),N_bits)/sqrt(M_AC);
                    nn_ac_sw = (ii_ac_h-1)*I_AC_v_bar + ii_ac_v;
                    for kk_AC_angle_ll = 1:K_sub
                        Equ_a_BS_lk_init = a_vec_BS_init_be_sq_set(:,kk_AC_angle_ll,ll_ac_sw)'*f_RF_BS_init_ll;
                        Equ_a_AC_lk_init = w_RF_AC_init_ll_ii'*a_vec_AC_init_be_sq_set(:,kk_AC_angle_ll,ll_ac_sw);
                        H_AC_angle_Sw(kk_AC_angle_ll,nn_ac_sw,ll_ac_sw) = Equ_a_AC_lk_init*Equ_a_BS_lk_init;
                    end
                end
            end
        end
        
        Y_AC_ang_com = zeros(N_OFDM_AC_angle,K_sub,L);
        for ll_ang_ac = 1:L
            for kk_ac_ll = 1:K_sub
                for nn_ac_ll = 1:N_OFDM_AC_angle
                    Y_AC_ang_com(nn_ac_ll,kk_ac_ll,ll_ang_ac) = alpha_init(ll_ang_ac)*H_AC_angle_Sw(kk_ac_ll,nn_ac_ll,ll_ang_ac)*...
                        s_AC_angle_set(kk_ac_ll,ll_ang_ac) + awgn_en*sigma_no*(normrnd(0,1) + 1i*normrnd(0,1))/sqrt(2);
                end
            end
        end
        Est_Azi_AC_set = zeros(L,1); Est_Ele_AC_set = zeros(L,1);
        for ll_ac_pow = 1:L
            Y_AC_ang_ll = Y_AC_ang_com(:,:,ll_ac_pow);
            [maximum,index_max] = max(abs(sum(Y_AC_ang_ll,2).^2));
            ii_ac_h_hat = ceil(index_max/I_AC_v_bar);
            ii_ac_v_hat = index_max - (ii_ac_h_hat-1)*I_AC_v_bar;
            Est_Azi_AC_set(ll_ac_pow) = Azi_AC_set_sw(ii_ac_h_hat);
            Est_Ele_AC_set(ll_ac_pow) = Ele_AC_set_sw(ii_ac_v_hat);
        end
        
        MSE_Azi_AC_set_temp_1_Sw(ii) = MSE_Azi_AC_set_temp_1_Sw(ii) + norm(Est_Azi_AC_set - Azi_AC_set_init)^2/L;
        MSE_Ele_AC_set_temp_1_Sw(ii) = MSE_Ele_AC_set_temp_1_Sw(ii) + norm(Est_Ele_AC_set - Ele_AC_set_init)^2/L;
        
        %%
        if mod(iter,1) == 0
            toc
            disp(['    SNR_dB = ' num2str(SNR_dBs(ii)) ',   iter = ' num2str(iter)])
        end
    end
    disp(['  Finish SNR_dB = ' num2str(SNR_dBs(ii))])
end
MSE_Azi_AC_set_1_Sw_2 = MSE_Azi_AC_set_temp_1_Sw/iterMax;
MSE_Ele_AC_set_1_Sw_2 = MSE_Ele_AC_set_temp_1_Sw/iterMax;
disp('Finish  All')

MarkerSize = 3;
LineWidth = 1.2;
Fontsize = 15;
RMSE_Azi_AC_set_1_Sw_2 = sqrt(MSE_Azi_AC_set_1_Sw_2); RMSE_Ele_AC_set_1_Sw_2 = sqrt(MSE_Ele_AC_set_1_Sw_2);

figure
semilogy(SNR_dBs,RMSE_Azi_AC_set_1_Sw_2,':ro','LineWidth',LineWidth); hold on; grid on;
semilogy(SNR_dBs,RMSE_Ele_AC_set_1_Sw_2,'-ro','LineWidth',LineWidth);
axis([-50 30 1e-5 3e0])
xlabel('SNR [dB]','Fontsize',Fontsize),ylabel('RMSE','Fontsize',Fontsize);
title('AC Azi & Ele RMSE','Fontsize',Fontsize);
set(gca, 'GridLineStyle', '-.','FontSize',Fontsize, 'linewidth',1.5,'Fontname','Times New Roman');
set(gcf, 'position', [700 300 650 550]); axis normal;
h2 = legend('Beam Sweeping method + \theta_{AC}','Beam Sweeping method + \phi_{AC}','Location','northeast');
set(h2,'Fontsize',11);
