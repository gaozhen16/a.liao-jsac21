clear, clc

load('Z1_Save_ACangle_CRB_v0_noBD_V2_20');
RCRB_Azi_AC_set_20 = sqrt(CRB_Azi_AC_set_20);
RCRB_Ele_AC_set_20 = sqrt(CRB_Ele_AC_set_20);
num_dB = length(SNR_dBs);
load('Z1_Save_ACangle_CRB_v0_BD_V1_2');
RCRB_Azi_AC_set_v1 = sqrt(CRB_Azi_AC_set);
RMSE_Azi_AC_set_0_v1 = sqrt(MSE_Azi_AC_set_0);
RMSE_Azi_AC_set_1_v1 = sqrt(MSE_Azi_AC_set_1);
RMSE_Azi_AC_set_2_v1 = sqrt(MSE_Azi_AC_set_2);
RMSE_Azi_AC_set_3_v1 = sqrt(MSE_Azi_AC_set_3);
RCRB_Ele_AC_set_v1 = sqrt(CRB_Ele_AC_set);
RMSE_Ele_AC_set_0_v1 = sqrt(MSE_Ele_AC_set_0);
RMSE_Ele_AC_set_1_v1 = sqrt(MSE_Ele_AC_set_1);
RMSE_Ele_AC_set_2_v1 = sqrt(MSE_Ele_AC_set_2);
RMSE_Ele_AC_set_3_v1 = sqrt(MSE_Ele_AC_set_3);
load('Z1_Save_ACangle_CRB_v0_BD_V2_4');
RCRB_Azi_AC_set_v2 = sqrt(CRB_Azi_AC_set);
RCRB_Azi_AC_set_v2(num_dB) = RCRB_Azi_AC_set_20(num_dB);
RCRB_Ele_AC_set_v2 = sqrt(CRB_Ele_AC_set);
RCRB_Ele_AC_set_v2(num_dB) = RCRB_Ele_AC_set_20(num_dB);
RMSE_Azi_AC_set_0_v2 = sqrt(MSE_Azi_AC_set_0);
RMSE_Azi_AC_set_1_v2 = sqrt(MSE_Azi_AC_set_1);
RMSE_Azi_AC_set_2_v2 = sqrt(MSE_Azi_AC_set_2);
RMSE_Azi_AC_set_3_v2 = sqrt(MSE_Azi_AC_set_3);
RMSE_Ele_AC_set_0_v2 = sqrt(MSE_Ele_AC_set_0);
RMSE_Ele_AC_set_1_v2 = sqrt(MSE_Ele_AC_set_1);
RMSE_Ele_AC_set_2_v2 = sqrt(MSE_Ele_AC_set_2);
RMSE_Ele_AC_set_3_v2 = sqrt(MSE_Ele_AC_set_3);
load('Z1_Save_ACangle_BeSq_Sw_V1');
RMSE_Azi_AC_set_1_Sw_1 = sqrt(MSE_Azi_AC_set_1_Sw_1);
RMSE_Ele_AC_set_1_Sw_1 = sqrt(MSE_Ele_AC_set_1_Sw_1);
load('Z1_Save_ACangle_BeSq_Sw_V2');
RMSE_Azi_AC_set_1_Sw_2 = sqrt(MSE_Azi_AC_set_1_Sw_2);
RMSE_Ele_AC_set_1_Sw_2 = sqrt(MSE_Ele_AC_set_1_Sw_2);

MarkerSize = 8;
LineWidth = 2;
Fontsize = 14;

%% 
figure
c15 = semilogy(SNR_dBs,RMSE_Azi_AC_set_1_Sw_1,'x','Color',[1 0.5 0],'LineWidth',LineWidth,'MarkerSize',MarkerSize); hold on; grid on;
c11 = semilogy(SNR_dBs,RMSE_Azi_AC_set_1_v1,'bd','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c12 = semilogy(SNR_dBs,RMSE_Azi_AC_set_2_v1,'rs','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c13 = semilogy(SNR_dBs,RMSE_Azi_AC_set_3_v1,'m^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c14 = semilogy(SNR_dBs,RMSE_Azi_AC_set_0_v1,'cp','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c21 = semilogy(SNR_dBs,RCRB_Azi_AC_set_v1,'-k','LineWidth',LineWidth + 0.5);
semilogy(SNR_dBs,RMSE_Azi_AC_set_1_Sw_1,'-x','Color',[1 0.5 0],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Azi_AC_set_0_v1,'-cp','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Azi_AC_set_3_v1,'-m^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Azi_AC_set_2_v1,'-rs','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Azi_AC_set_1_v1,'-bd','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c22 = semilogy(SNR_dBs(1:num_dB),RCRB_Azi_AC_set_20,'--k','LineWidth',LineWidth + 0.5);
semilogy(SNR_dBs,RMSE_Azi_AC_set_1_Sw_2,'--x','Color',[1 0.5 0],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Azi_AC_set_0_v2,'--cp','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Azi_AC_set_3_v2,'--m^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Azi_AC_set_2_v2,'--rs','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Azi_AC_set_1_v2,'--bd','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs(num_dB:end),RCRB_Azi_AC_set_v2(num_dB:end),'--k','LineWidth',LineWidth + 0.5);
axis([-50 30 2e-8 4e0])
set(gca,'ytick',[1e-6 1e-4 1e-2 1e0 1e2]);
set(gca, 'GridLineStyle', '--','FontSize',Fontsize, 'linewidth',1.0,'Fontname','Times New Roman');
xlabel({'SNR [dB]';'(a) azimuth angle {\it{¦È}}^{AC}'},'Fontsize',18); ylabel('RMSE','Fontsize',18)
set(gca,'position',[0.098 0.19 0.885 0.80])
set(gcf, 'position', [250 200 700 500]);
text(5.0565, 0.01549,'Method 1','color','k','FontSize',16,'Fontname','Times New Roman')
annotation('arrow',[0.5957143 0.70143],[0.936 0.784],'Color',[0.34 0.83 0.0745],'LineWidth',1.5);
annotation('arrow',[0.5528573 0.7],[0.676 0.75],'Color',[0.34 0.83 0.0745],'LineWidth',1.5);
annotation('ellipse',[0.525 0.62 0.03043 0.07],'Color',[0.34 0.83 0.0745],'LineWidth',1.5);
text(-37.44, 5.302686e-05,'Method 2','color','k','FontSize',16,'Fontname','Times New Roman')
annotation('arrow',[0.492857 0.34143],[0.938 0.564],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','--');
annotation('arrow',[0.52143 0.3657],[0.564 0.532],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','--');
annotation('ellipse',[0.525 0.496 0.03043 0.12],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','--');
axis normal;
h11 = legend([c15,c12,c13,c11,c14],'Beam sweeping method [56]','Proposed Algorithm 1, {{\it{i}}_{AC}^{max}} = 1',...
    'Proposed Algorithm 1, {{\it{i}}_{AC}^{max}} = 2','Conventional scheme, no TTDU module',...
    'Conventional scheme, ideal TTDU module','Location','southwest');
ah1 = axes('position', get(gca,'position'), 'visible','off'); set(h11, 'FontSize',Fontsize-1);
h12 = legend(ah1, [c21,c22], 'CRLB, method 1','CRLB, method 2');
set(h12,'Fontsize',Fontsize-1,'Fontname','Times New Roman');
set(gca, 'Position',[0.64 0.224 0.243 0.09], 'linewidth',1.0);

%% 
figure
c35 = semilogy(SNR_dBs,RMSE_Ele_AC_set_1_Sw_1,'x','Color',[1 0.5 0],'LineWidth',LineWidth,'MarkerSize',MarkerSize); hold on; grid on;
c31 = semilogy(SNR_dBs,RMSE_Ele_AC_set_1_v1,'bd','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c32 = semilogy(SNR_dBs,RMSE_Ele_AC_set_2_v1,'rs','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c33 = semilogy(SNR_dBs,RMSE_Ele_AC_set_3_v1,'m^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c34 = semilogy(SNR_dBs,RMSE_Ele_AC_set_0_v1,'cp','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c41 = semilogy(SNR_dBs,RCRB_Ele_AC_set_v1,'-k','LineWidth',LineWidth + 0.5);
semilogy(SNR_dBs,RMSE_Ele_AC_set_1_Sw_1,'-x','Color',[1 0.5 0],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Ele_AC_set_0_v1,'-cp','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Ele_AC_set_3_v1,'-m^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Ele_AC_set_2_v1,'-rs','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Ele_AC_set_1_v1,'-bd','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c42 = semilogy(SNR_dBs(1:num_dB),RCRB_Ele_AC_set_20,'--k','LineWidth',LineWidth + 0.5);
semilogy(SNR_dBs,RMSE_Ele_AC_set_1_Sw_2,'--x','Color',[1 0.5 0],'LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Ele_AC_set_0_v2,'--cp','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Ele_AC_set_3_v2,'--m^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Ele_AC_set_2_v2,'--rs','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs,RMSE_Ele_AC_set_1_v2,'--bd','LineWidth',LineWidth,'MarkerSize',MarkerSize);
semilogy(SNR_dBs(num_dB:end),RCRB_Ele_AC_set_v2(num_dB:end),'--k','LineWidth',LineWidth + 0.5);
axis([-50 30 2e-8 4e0])
set(gca,'ytick',[1e-6 1e-4 1e-2 1e0 1e2]);
set(gca, 'GridLineStyle', '--','FontSize',Fontsize, 'linewidth',1.0,'Fontname','Times New Roman');
xlabel({'SNR [dB]';'(b) elevation angle {\it{¦Õ}}^{AC}'},'Fontsize',18); ylabel('RMSE','Fontsize',18)
set(gca,'position',[0.098 0.19 0.885 0.80])
set(gcf, 'position', [960 200 700 500]);
text(5.0565, 0.01549,'Method 1','color','k','FontSize',16,'Fontname','Times New Roman')
annotation('arrow',[0.5957143 0.70143],[0.936 0.784],'Color',[0.34 0.83 0.0745],'LineWidth',1.5);
annotation('arrow',[0.55143 0.7],[0.688 0.75],'Color',[0.34 0.83 0.0745],'LineWidth',1.5);
annotation('ellipse',[0.525 0.632 0.03043 0.068],'Color',[0.34 0.83 0.0745],'LineWidth',1.5);
text(-37.44, 5.302686e-05,'Method 2','color','k','FontSize',16,'Fontname','Times New Roman')
annotation('arrow',[0.492857 0.34143],[0.938 0.564],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','--');
annotation('arrow',[0.52143 0.3657],[0.578 0.532],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','--');
annotation('ellipse',[0.525 0.496 0.03043 0.134],'Color',[0.34 0.83 0.0745],'LineWidth',1.5,'LineStyle','--');
axis normal;
h21 = legend([c35,c32,c33,c31,c34],'Beam sweeping method [56]','Proposed Algorithm 1, {{\it{i}}_{AC}^{max}} = 1',...
    'Proposed Algorithm 1, {{\it{i}}_{AC}^{max}} = 2','Conventional scheme, no TTDU module',...
    'Conventional scheme, ideal TTDU module','Location','southwest');
ah2 = axes('position', get(gca,'position'), 'visible','off'); set(h21, 'FontSize',Fontsize-1);
h22 = legend(ah2, [c41,c42], 'CRLB, method 1','CRLB, method 2');
set(h22,'Fontsize',Fontsize-1,'Fontname','Times New Roman');
set(gca, 'Position',[0.64 0.224 0.243 0.09], 'linewidth',1.0);
