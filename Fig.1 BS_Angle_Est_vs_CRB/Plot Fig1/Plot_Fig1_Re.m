clear, clc

load('Z1_Save_BSangle_CRB_v0_BD_V1');
load('Z1_Save_BSangle_CRB_v0_noBD_V1');
load('Z1_Save_BSangle_BeSq_Sw_V1');

MarkerSize = 8;
LineWidth = 2;
Fontsize = 14;

RMSE_Azi_BS_set_1_Sw = sqrt(MSE_Azi_BS_set_1_Sw);
RMSE_Azi_BS_set_BD_1 = sqrt(MSE_Azi_BS_set_BD_1);
RMSE_Azi_BS_set_BD_2 = sqrt(MSE_Azi_BS_set_BD_2);
RMSE_Azi_BS_set_BD_3 = sqrt(MSE_Azi_BS_set_BD_3);
RMSE_Azi_BS_set = sqrt(MSE_Azi_BS_set);
RCRB_Azi_BS_set = sqrt(CRB_Azi_BS_set);

figure
c16 = semilogy(SNR_dBs,RMSE_Azi_BS_set_1_Sw,'-x','Color',[1 0.5 0],'LineWidth',LineWidth,'MarkerSize',MarkerSize); hold on; grid on;
c15 = semilogy(SNR_dBs,RCRB_Azi_BS_set,'-k','LineWidth',LineWidth + 0.5);
c14 = semilogy(SNR_dBs,RMSE_Azi_BS_set,'-cp','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c13 = semilogy(SNR_dBs,RMSE_Azi_BS_set_BD_1,'-bd','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c12 = semilogy(SNR_dBs,RMSE_Azi_BS_set_BD_3,'-m^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c11 = semilogy(SNR_dBs,RMSE_Azi_BS_set_BD_2,'-rs','LineWidth',LineWidth,'MarkerSize',MarkerSize);
axis([-50 30 3e-6 3e0])
set(gca,'ytick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0]);
set(gca, 'GridLineStyle', '--','FontSize',Fontsize, 'linewidth',1.0,'Fontname','Times New Roman');
xlabel({'SNR [dB]';'(a) azimuth angle {\it{¦È}}^{BS}'},'Fontsize',18); ylabel('RMSE','Fontsize',18);
set(gcf, 'position', [250 200 700 500]);
set(gca,'position',[0.098 0.19 0.885 0.80])
axis normal;
h1 = legend([c16,c11,c12,c13,c14,c15],'Beam sweeping method [56]','Proposed Algorithm 1, {{\it{i}}_{BS}^{max}} = 1',...
    'Proposed Algorithm 1, {{\it{i}}_{BS}^{max}} = 2','Conventional scheme, no TTDU module',...
    'Conventional scheme, ideal TTDU module','CRLB','Location','southwest');
set(h1,'Fontsize',Fontsize-1);

RMSE_Ele_BS_set_1_Sw = sqrt(MSE_Ele_BS_set_1_Sw);
RMSE_Ele_BS_set_BD_1 = sqrt(MSE_Ele_BS_set_BD_1);
RMSE_Ele_BS_set_BD_2 = sqrt(MSE_Ele_BS_set_BD_2);
RMSE_Ele_BS_set_BD_3 = sqrt(MSE_Ele_BS_set_BD_3);
RMSE_Ele_BS_set = sqrt(MSE_Ele_BS_set);
RCRB_Ele_BS_set = sqrt(CRB_Ele_BS_set);

figure
c26 = semilogy(SNR_dBs,RMSE_Ele_BS_set_1_Sw,'-x','Color',[1 0.5 0],'LineWidth',LineWidth,'MarkerSize',MarkerSize); hold on; grid on;
c25 = semilogy(SNR_dBs,RCRB_Ele_BS_set,'-k','LineWidth',LineWidth + 0.5);
c24 = semilogy(SNR_dBs,RMSE_Ele_BS_set,'-cp','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c23 = semilogy(SNR_dBs,RMSE_Ele_BS_set_BD_1,'-bd','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c22 = semilogy(SNR_dBs,RMSE_Ele_BS_set_BD_3,'-m^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c21 = semilogy(SNR_dBs,RMSE_Ele_BS_set_BD_2,'-rs','LineWidth',LineWidth,'MarkerSize',MarkerSize);
axis([-50 30 3e-6 3e0])
set(gca,'ytick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0]);
set(gca, 'GridLineStyle', '--','FontSize',Fontsize, 'linewidth',1.0,'Fontname','Times New Roman');
xlabel({'SNR [dB]';'(b) elevation angle {\it{¦Õ}}^{BS}'},'Fontsize',18); ylabel('RMSE','Fontsize',18)
set(gcf, 'position', [960 200 700 500]);
set(gca,'position',[0.098 0.19 0.885 0.80])
axis normal;
h2 = legend([c26,c21,c22,c23,c24,c25],'Beam sweeping method [56]','Proposed Algorithm 1, {{\it{i}}_{BS}^{max}} = 1',...
    'Proposed Algorithm 1, {{\it{i}}_{BS}^{max}} = 2','Conventional scheme, no TTDU module',...
    'Conventional scheme, ideal TTDU module','CRLB','Location','southwest');
set(h2,'Fontsize',Fontsize-1);
