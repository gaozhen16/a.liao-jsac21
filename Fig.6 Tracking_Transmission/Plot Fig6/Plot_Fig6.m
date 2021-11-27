clear, clc

MarkerSize = 8;
LineWidth = 3;
Fontsize = 14;

%% 
load('Z1_Save_Data_Track_v0_AMP_V1_70_2');
N_OFDM_C_used = N_OFDM_C;
N_OFDM_tot_used = N_OFDM_C_used*floor(N_CCT_used/10)*10;
n_OFDM_tot_used = 1:N_OFDM_tot_used;
figure
c11 = plot(n_OFDM_tot_used,Amp_set_notrack_used(1,1:N_OFDM_tot_used),'c-','LineWidth',LineWidth); hold on; grid on;
c14 = plot(n_OFDM_tot_used,Amp_set_Real_eff_used(1,1:N_OFDM_tot_used),'-b','LineWidth',LineWidth+1);
c12 = plot(n_OFDM_tot_used,Amp_set_EDD_used(1,1:N_OFDM_tot_used),'-r','LineWidth',LineWidth);
c21 = plot(n_OFDM_tot_used,Amp_set_notrack_used(2,1:N_OFDM_tot_used),'g-','LineWidth',LineWidth);
c24 = plot(n_OFDM_tot_used,Amp_set_Real_eff_used(2,1:N_OFDM_tot_used),'-k','LineWidth',LineWidth+1);
c22 = plot(n_OFDM_tot_used,Amp_set_EDD_used(2,1:N_OFDM_tot_used),'-m','LineWidth',LineWidth);
axis([1 N_OFDM_tot_used 0 3.8e4])
set(gca,'xtick',(10:10:100)*N_OFDM_C_used);
set(gca,'xticklabels',{'10','20','30','40','50','60','70','80','90','100'});
set(gca, 'GridLineStyle', '--','FontSize',Fontsize, 'linewidth',1.0,'Fontname','Times New Roman');
xlabel({'Number of TI';'(a)'},'Fontsize',18); ylabel('Amplitude of Effective Channels','Fontsize',18)
set(gcf, 'position', [250 200 700 500]);
set(gca,'position',[0.089 0.175 0.885 0.775])
axis normal;
h11 = legend([c11,c12,c14], 'Initial channel estimation, BS #1','DADD-based channel tracking, BS #1',...
    'Real-time effective channel, BS #1','Location','southwest');
ah = axes('position', get(gca,'position'), 'visible','off');
h12 = legend(ah, [c21,c22,c24], 'Initial channel estimation, BS #2','DADD-based channel tracking, BS #2',...
    'Real-time effective channel, BS #2');
set(h12, 'Position',[0.49 0.659 0.452856 0.158], 'FontSize',16);
set(h11,'Fontsize',Fontsize); set(h12,'Fontsize',Fontsize,'Fontname','Times New Roman'); set(gca, 'linewidth',1.0);
h_small_11 = axes('position',[0.16 0.395 0.245 0.175]);
axis(h_small_11);
plot(n_OFDM_tot_used(29*N_OFDM_C_used+N_OFDM_C_used/10:30*N_OFDM_C_used-N_OFDM_C_used/10),...
    Amp_set_Real_eff_used(1,29*N_OFDM_C_used+N_OFDM_C_used/10:30*N_OFDM_C_used-N_OFDM_C_used/10),...
    '-b','LineWidth',LineWidth); hold on; grid on;
plot(n_OFDM_tot_used(29*N_OFDM_C_used+N_OFDM_C_used/10:30*N_OFDM_C_used-N_OFDM_C_used/10),...
    Amp_set_EDD_used(1,29*N_OFDM_C_used+N_OFDM_C_used/10:30*N_OFDM_C_used-N_OFDM_C_used/10),'-r','LineWidth',LineWidth-1.5);
set(gca,'xtick',[29*N_OFDM_C_used+N_OFDM_C_used/10 29*N_OFDM_C_used+N_OFDM_C_used/2 30*N_OFDM_C_used-N_OFDM_C_used/10]);
xlim([29*N_OFDM_C_used+N_OFDM_C_used/10,30*N_OFDM_C_used-N_OFDM_C_used/10]);
set(gca,'xticklabels',{'29.1','29.5','29.9'}, 'linewidth',1.0);
annotation('line',[0.347143 0.2843],[0.66 0.568],'LineWidth',1.0,'LineStyle','-');
h_small_12 = axes('position',[0.670 0.447 0.245 0.175]);
axis(h_small_12);
plot(n_OFDM_tot_used(69*N_OFDM_C_used+N_OFDM_C_used/10:70*N_OFDM_C_used-N_OFDM_C_used/10),...
    Amp_set_Real_eff_used(2,69*N_OFDM_C_used+N_OFDM_C_used/10:70*N_OFDM_C_used-N_OFDM_C_used/10),...
    '-k','LineWidth',LineWidth); hold on; grid on;
plot(n_OFDM_tot_used(69*N_OFDM_C_used+N_OFDM_C_used/10:70*N_OFDM_C_used-N_OFDM_C_used/10),...
    Amp_set_EDD_used(2,69*N_OFDM_C_used+N_OFDM_C_used/10:70*N_OFDM_C_used-N_OFDM_C_used/10),'-m','LineWidth',LineWidth-1.5);
set(gca,'xtick',[69*N_OFDM_C_used+N_OFDM_C_used/10 69*N_OFDM_C_used+N_OFDM_C_used/2 70*N_OFDM_C_used-N_OFDM_C_used/10]);
xlim([69*N_OFDM_C_used+N_OFDM_C_used/10,70*N_OFDM_C_used-N_OFDM_C_used/10]);
set(gca,'xticklabels',{'69.1','69.5','69.9'}, 'linewidth',1.0);
annotation('line',[0.793 0.7043],[0.45 0.326],'LineWidth',1.0,'LineStyle','-');

%% 
load('Z1_Save_Data_Track_v0_NMSE_V1__20dB_70');
N_OFDM_C_used = N_OFDM_C;
N_OFDM_tot_used = N_OFDM_C_used*floor(N_CCT_used/10)*10;
n_OFDM_tot_used = 1:N_OFDM_tot_used;
NMSE_set_Notr_used__20dB = NMSE_set_Notr_used_dB(1:N_OFDM_tot_used);
NMSE_set_DADD__20dB = NMSE_set_EDD_used_dB(1:N_OFDM_tot_used);
load('Z1_Save_Data_Track_v0_NMSE_V1__10dB_70_2');
NMSE_set_Notr_used__10dB = NMSE_set_Notr_used_dB(1:N_OFDM_tot_used);
NMSE_set_DADD__10dB = NMSE_set_EDD_used_dB(1:N_OFDM_tot_used);
load('Z1_Save_Data_Track_v0_NMSE_V1_0dB_70_1');
NMSE_set_Notr_used_0dB = NMSE_set_Notr_used_dB(1:N_OFDM_tot_used);
load('Z1_Save_Data_Track_v0_NMSE_V1_0dB_70_4');
NMSE_set_DADD_0dB = NMSE_set_EDD_used_dB(1:N_OFDM_tot_used);
figure
plot(n_OFDM_tot_used,NMSE_set_Notr_used__20dB,'-g','LineWidth',LineWidth); hold on; grid on;
plot(n_OFDM_tot_used,NMSE_set_Notr_used__10dB,'-c','LineWidth',LineWidth);
plot(n_OFDM_tot_used,NMSE_set_Notr_used_0dB,'-m','LineWidth',LineWidth);
plot(n_OFDM_tot_used,NMSE_set_DADD__20dB,'-b','LineWidth',LineWidth);
plot(n_OFDM_tot_used,NMSE_set_DADD__10dB,'-r','LineWidth',LineWidth);
plot(n_OFDM_tot_used,NMSE_set_DADD_0dB,'-k','LineWidth',LineWidth);
axis([1 N_OFDM_tot_used -85 35])
set(gca,'xtick',(10:10:100)*N_OFDM_C_used);
set(gca,'xticklabels',{'10','20','30','40','50','60','70','80','90','100'});
set(gca,'ytick',-70:20:30);
set(gca, 'GridLineStyle', '--','FontSize',Fontsize, 'linewidth',1.0,'Fontname','Times New Roman');
xlabel({'Number of TI';'(b)'},'Fontsize',18); ylabel('NMSE of Effective Channels [dB]','Fontsize',18)
set(gcf, 'position', [960 200 700 500]);
set(gca,'position',[0.089 0.175 0.885 0.775])
axis normal;
h1 = legend('Initial channel estimation, SNR = -20 dB','Initial channel estimation, SNR = -10 dB',...
    'Initial channel estimation, SNR = 0 dB','DADD-based channel tracking, SNR = -20 dB',...
    'DADD-based channel tracking, SNR = -10 dB','DADD-based channel tracking, SNR = 0 dB');
set(h1, 'Position',[0.445 0.45 0.472 0.309], 'Fontsize',Fontsize-1);
h_small_21 = axes('position',[0.162 0.54 0.245 0.1]);
axis(h_small_21);
plot(n_OFDM_tot_used(19*50:21*50),NMSE_set_DADD__20dB(19*50:21*50),'-b','LineWidth',LineWidth-1.5);
set(gca, 'linewidth',1.3,'XTickLabel','');
h_small_22 = axes('position',[0.162 0.43 0.245 0.1]);
axis(h_small_22);
plot(n_OFDM_tot_used(19*50:21*50),NMSE_set_DADD__10dB(19*50:21*50),'-r','LineWidth',LineWidth-1.5);
set(gca, 'linewidth',1.0,'XTickLabel','');
h_small_23 = axes('position',[0.162 0.32 0.245 0.1]);
axis(h_small_23);
plot(n_OFDM_tot_used(19*50:21*50),NMSE_set_DADD_0dB(19*50:21*50),'-k','LineWidth',LineWidth-1.5);
set(gca,'xtick',[950 1000 1050]);
set(gca,'xticklabels',{'19','20','21'}, 'linewidth',1.0);
annotation('line',[0.2857 0.2643],[0.322 0.288],'LineWidth',1.0,'LineStyle','-');
annotation('ellipse',[0.2543 0.222 0.02 0.068],'LineWidth',1.8,'LineStyle','-');
