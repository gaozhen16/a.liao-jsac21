clear, clc

load('Z1_Save_Doppler_CRB_v0_BD_V1_1');

MarkerSize = 8;
LineWidth = 2;
Fontsize = 14;

%% 
RMSE_Doppler_set_1 = sqrt(MSE_Doppler_set_1);
RMSE_Doppler_set_2 = sqrt(MSE_Doppler_set_2);
RMSE_Doppler_set_3 = sqrt(MSE_Doppler_set_3);
RMSE_Doppler_set_0 = sqrt(MSE_Doppler_set_0);
RCRB_Doppler_set = sqrt(CRB_Doppler_set);

figure
c15 = semilogy(SNR_dBs,RCRB_Doppler_set,'-k','LineWidth',LineWidth + 0.5); hold on; grid on;
c14 = semilogy(SNR_dBs,RMSE_Doppler_set_0,'-cp','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c13 = semilogy(SNR_dBs,RMSE_Doppler_set_1,'-bd','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c12 = semilogy(SNR_dBs,RMSE_Doppler_set_3,'-m^','LineWidth',LineWidth,'MarkerSize',MarkerSize);
c11 = semilogy(SNR_dBs,RMSE_Doppler_set_2,'-rs','LineWidth',LineWidth,'MarkerSize',MarkerSize);
axis([-120 20 6e-4 2e5])
set(gca,'ytick',[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5]);
set(gca, 'GridLineStyle', '--','FontSize',Fontsize, 'linewidth',1.0,'Fontname','Times New Roman');
xlabel('SNR [dB]','Fontsize',18); ylabel('RMSE','Fontsize',18)
set(gcf, 'position', [600 250 700 500]);
set(gca,'position',[0.098 0.12 0.883 0.87])
axis normal;
h1 = legend([c11,c12,c13,c14,c15],'Proposed Algorithm 2, {{\it{i}}_{do}^{max}} = 1',...
    'Proposed Algorithm 2, {{\it{i}}_{do}^{max}} = 2','Conventional scheme, {{\it{i}}_{do}^{max}} = 0',...
    'Conventional scheme, no Doppler squint','CRLB','Location','northeast');
set(h1,'Fontsize',Fontsize);
