clear, clc

MarkerSize = 8;
LineWidth = 2;
Fontsize = 14;

%% 
load('Z1_Save_BSangle_CRB_Track_v0_Omega1_V1');
use_plot_bs = length(SNR_dBs)-1;
SNR_dBs_BS_used = SNR_dBs(1:use_plot_bs);
RMSE_Azi_BS_set_trac_0_Omega1 = sqrt(MSE_Azi_BS_set_trac_0(1:use_plot_bs));
RMSE_Azi_BS_set_trac_1_Omega1 = sqrt(MSE_Azi_BS_set_trac_1(1:use_plot_bs));
RCRB_Azi_BS_set_trac_Omega1 = sqrt(CRB_Azi_BS_set_trac(1:use_plot_bs));
load('Z1_Save_BSangle_CRB_Track_v0_Omega4_V1');
RMSE_Azi_BS_set_trac_0_Omega4 = sqrt(MSE_Azi_BS_set_trac_0(1:use_plot_bs));
RMSE_Azi_BS_set_trac_1_Omega4 = sqrt(MSE_Azi_BS_set_trac_1(1:use_plot_bs));
RCRB_Azi_BS_set_trac_Omega4 = sqrt(CRB_Azi_BS_set_trac(1:use_plot_bs));
load('Z1_Save_BSangle_CRB_Track_v0_Omega16_V1');
RMSE_Azi_BS_set_trac_0_Omega16 = sqrt(MSE_Ele_BS_set_trac_0(1:use_plot_bs));
RMSE_Azi_BS_set_trac_1_Omega16 = sqrt(MSE_Ele_BS_set_trac_1(1:use_plot_bs));
RCRB_Azi_BS_set_trac_Omega16 = sqrt(CRB_Azi_BS_set_trac(1:use_plot_bs));

figure
semilogy(SNR_dBs_BS_used,RMSE_Azi_BS_set_trac_1_Omega1,'-bp','LineWidth',LineWidth+1,'MarkerSize',MarkerSize); hold on; grid on;
semilogy(SNR_dBs_BS_used,RMSE_Azi_BS_set_trac_0_Omega1,'-rs','LineWidth',LineWidth-0.5,'MarkerSize',MarkerSize);
semilogy(SNR_dBs_BS_used,RCRB_Azi_BS_set_trac_Omega1,'-k','LineWidth',LineWidth + 0.5);
semilogy(SNR_dBs_BS_used,RMSE_Azi_BS_set_trac_1_Omega4,'-bp','LineWidth',LineWidth+1,'MarkerSize',MarkerSize);
semilogy(SNR_dBs_BS_used,RMSE_Azi_BS_set_trac_0_Omega4,'-rs','LineWidth',LineWidth-0.5,'MarkerSize',MarkerSize);
semilogy(SNR_dBs_BS_used,RCRB_Azi_BS_set_trac_Omega4,'-k','LineWidth',LineWidth + 0.5);
semilogy(SNR_dBs_BS_used,RMSE_Azi_BS_set_trac_1_Omega16,'-bp','LineWidth',LineWidth+1,'MarkerSize',MarkerSize);
semilogy(SNR_dBs_BS_used,RMSE_Azi_BS_set_trac_0_Omega16,'-rs','LineWidth',LineWidth-0.5,'MarkerSize',MarkerSize);
semilogy(SNR_dBs_BS_used,RCRB_Azi_BS_set_trac_Omega16,'-k','LineWidth',LineWidth + 0.5);
axis([-100 -20 5e-7 2e0])
set(gca,'ytick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0]);
set(gca, 'GridLineStyle', '--','FontSize',Fontsize, 'linewidth',1.0,'Fontname','Times New Roman');
xlabel({'SNR [dB]';'(a) azimuth angle at BS: {\it{¦È}}^{BS}'},'Fontsize',18); ylabel('RMSE','Fontsize',18)
set(gcf, 'position', [250 200 700 500]);
set(gca,'position',[0.095 0.19 0.885 0.80])
axis normal;
h1 = legend('Proposed Algorithm 1, {{\it{i}}_{BS}^{max}} = 2', 'Conventional scheme, ideal TTDU module', 'CRLB','Location','southwest'); % {\fontname{Times New Roman}}
set(h1,'Fontsize',Fontsize);
text(-54.883, 0.0504,'{\it{\Omega}} = 1','color','k','FontSize',16,'Fontname','Times New Roman')
annotation('arrow',[0.55 0.6043],[0.666 0.77],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
annotation('ellipse',[0.522 0.608 0.03043 0.072],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
text(-42.032, 0.01201,'{\it{\Omega}} = 4','color','k','FontSize',16,'Fontname','Times New Roman')
annotation('arrow',[0.6586 0.74143],[0.524 0.7],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
annotation('ellipse',[0.632 0.46 0.03043 0.072],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
text(-30.03, 0.003062,'{\it{\Omega}} = 16','color','k','FontSize',16,'Fontname','Times New Roman')
annotation('arrow',[0.767 0.8786],[0.396 0.624],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
annotation('ellipse',[0.7436 0.328 0.0305 0.072],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');

%%
load('Z1_Save_ACangle_CRB_Track_v0_Omega1_V1');
use_plot_ac = length(SNR_dBs);
SNR_dBs_AC_used = SNR_dBs(1:use_plot_ac);
RMSE_Azi_AC_set_trac_0_Omega1 = sqrt(MSE_Azi_AC_set_trac_0(1:use_plot_ac));
RMSE_Azi_AC_set_trac_1_Omega1 = sqrt(MSE_Azi_AC_set_trac_1(1:use_plot_ac));
RCRB_Azi_AC_set_trac_Omega1 = sqrt(CRB_Azi_AC_set_trac(1:use_plot_ac));
load('Z1_Save_ACangle_CRB_Track_v0_Omega4_V1');
RMSE_Azi_AC_set_trac_0_Omega4 = sqrt(MSE_Azi_AC_set_trac_0(1:use_plot_ac));
RMSE_Azi_AC_set_trac_1_Omega4 = sqrt(MSE_Azi_AC_set_trac_1(1:use_plot_ac));
RCRB_Azi_AC_set_trac_Omega4 = sqrt(CRB_Azi_AC_set_trac(1:use_plot_ac));
load('Z1_Save_ACangle_CRB_Track_v0_Omega16_V1');
RMSE_Azi_AC_set_trac_0_Omega16 = sqrt(MSE_Azi_AC_set_trac_0(1:use_plot_ac));
RMSE_Azi_AC_set_trac_1_Omega16 = sqrt(MSE_Azi_AC_set_trac_1(1:use_plot_ac));
RCRB_Azi_AC_set_trac_Omega16 = sqrt(CRB_Azi_AC_set_trac(1:use_plot_ac));
figure
semilogy(SNR_dBs_AC_used,RMSE_Azi_AC_set_trac_1_Omega1,'-bp','LineWidth',LineWidth+1,'MarkerSize',MarkerSize); hold on; grid on;
semilogy(SNR_dBs_AC_used,RMSE_Azi_AC_set_trac_0_Omega1,'-rs','LineWidth',LineWidth-0.5,'MarkerSize',MarkerSize);
semilogy(SNR_dBs_AC_used,RCRB_Azi_AC_set_trac_Omega1,'-k','LineWidth',LineWidth + 0.5);
semilogy(SNR_dBs_AC_used,RMSE_Azi_AC_set_trac_1_Omega4,'-bp','LineWidth',LineWidth+1,'MarkerSize',MarkerSize);
semilogy(SNR_dBs_AC_used,RMSE_Azi_AC_set_trac_0_Omega4,'-rs','LineWidth',LineWidth-0.5,'MarkerSize',MarkerSize);
semilogy(SNR_dBs_AC_used,RCRB_Azi_AC_set_trac_Omega4,'-k','LineWidth',LineWidth + 0.5);
semilogy(SNR_dBs_AC_used,RMSE_Azi_AC_set_trac_1_Omega16,'-bp','LineWidth',LineWidth+1,'MarkerSize',MarkerSize);
semilogy(SNR_dBs_AC_used,RMSE_Azi_AC_set_trac_0_Omega16,'-rs','LineWidth',LineWidth-0.5,'MarkerSize',MarkerSize);
semilogy(SNR_dBs_AC_used,RCRB_Azi_AC_set_trac_Omega16,'-k','LineWidth',LineWidth + 0.5);
axis([-120 -40 5e-7 2e0])
set(gca,'ytick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0]);
set(gca, 'GridLineStyle', '--','FontSize',Fontsize, 'linewidth',1.0,'Fontname','Times New Roman');
xlabel({'SNR [dB]';'(b) azimuth angle at aircraft: {\it{¦È}}^{AC}'},'Fontsize',18); ylabel('RMSE','Fontsize',18)
set(gcf, 'position', [960 200 700 500]);
set(gca,'position',[0.095 0.19 0.885 0.80])
axis normal;
h1 = legend('Proposed Algorithm 1, {{\it{i}}_{AC}^{max}} = 2', 'Conventional scheme, ideal TTDU module', 'CRLB','Location','southwest'); % {\fontname{Times New Roman}}
set(h1,'Fontsize',Fontsize);
text(-75.0123, 0.055,'{\it{\Omega}} = 1','color','k','FontSize',16,'Fontname','Times New Roman')
annotation('arrow',[0.55143 0.60143],[0.742 0.776],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
annotation('ellipse',[0.522 0.686 0.03043 0.08],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
text(-62.29, 0.01255,'{\it{\Omega}} = 4','color','k','FontSize',16,'Fontname','Times New Roman')
annotation('arrow',[0.6557 0.74143],[0.57 0.7],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
annotation('ellipse',[0.632 0.496 0.03043 0.08],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
text(-50.054, 0.0031,'{\it{\Omega}} = 16','color','k','FontSize',16,'Fontname','Times New Roman')
annotation('arrow',[0.767 0.8786],[0.376 0.622],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
annotation('ellipse',[0.74343 0.304 0.03043 0.08],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','-');
