
clear all;close all;clc;


%%
JMP=0.04;
JMP=0.025; %  This defines the grid resolution (smaller number higher resolution)
NOISE_TYPE_GAUSSIAN=1; NOISE_TYPE_STRECH=2;NOISE_TYPE_TRIMODAL=3;NOISE_TYPE_FIXED=4;NOISE_TYPE_FIXED_UNIFORM=5;

IS_DENOISE_ONLY=true; % NOTE !!! IMPORTANT (skip adding noise in the synthesis step, compatability with a version of DDPM that do not include noise)

% DECIDE WHAT TO DO:
QQ0=[]; % parameters shared between coditions.
QQ0.JMP=JMP;
QQ0.IS_PLOT=true;
QQ0.IS_SWISS_ROLL=true;


todoS=cell(1,1); %initialize conditions
tcnt=0;


QQ=QQ0;
T=10;


QQ.NOISE_TYPE=NOISE_TYPE_GAUSSIAN;%S02 A
QQ.IS_DENOISE_ONLY=false;% NOTE !!! IMPORTANT (do not skip adding noise in the synthesis step, compatability with a version of DDPM that do not include noise)
sigma0=0.01;sigma1=0.01; sigma_vec=linspace(sigma0,sigma1,T);
tcnt=tcnt+1;
todoS{tcnt,1}.QQ=QQ;
todoS{tcnt,1}.sigma_vec=sigma_vec;
todoS{tcnt,1}.T=T;
todoS{tcnt,1}.name='Gaussian (T=10)';
% %
QQ.NOISE_TYPE=NOISE_TYPE_STRECH; % S02 B
QQ.IS_DENOISE_ONLY=false;% NOTE !!! IMPORTANT (do not skip adding noise in the synthesis step, compatability with a version of DDPM that do not include noise)
sigma0=0.0001;sigma1=0.0001;sigma_vec=linspace(sigma0,sigma1,T);
tcnt=tcnt+1;
todoS{tcnt,1}.QQ=QQ;
todoS{tcnt,1}.sigma_vec=sigma_vec;
todoS{tcnt,1}.T=T;
todoS{tcnt,1}.name='Bimodal (T=10)';

QQ.NOISE_TYPE=NOISE_TYPE_FIXED; % S02 C
QQ.IS_DENOISE_ONLY=false;% NOTE !!! IMPORTANT (do not skip adding noise in the synthesis step, compatability with a version of DDPM that do not include noise)
sigma0=0.0001;sigma1=0.0001;sigma_vec=linspace(sigma0,sigma1,T);

tcnt=tcnt+1;
todoS{tcnt,1}.QQ=QQ;
todoS{tcnt,1}.sigma_vec=sigma_vec;
todoS{tcnt,1}.T=T;
todoS{tcnt,1}.name='Fade (T=10)';


T=100;
QQ.NOISE_TYPE=NOISE_TYPE_GAUSSIAN;%S02 A
QQ.IS_DENOISE_ONLY=false;% NOTE !!! IMPORTANT (do not skip adding noise in the synthesis step, compatability with a version of DDPM that do not include noise)
sigma0=0.01;sigma1=0.01; sigma_vec=linspace(sigma0,sigma1,T);
tcnt=tcnt+1;
todoS{tcnt,1}.QQ=QQ;
todoS{tcnt,1}.sigma_vec=sigma_vec;
todoS{tcnt,1}.T=T;
todoS{tcnt,1}.name='Gaussian noise (T=100)';
% % 


% RUN SIMULATIONS:
for II=1:length(todoS)
    QQ=todoS{II}.QQ;
    sigma_vec=todoS{II}.sigma_vec;
    T=todoS{II}.T;

    figure(101+II);clf;
    set(gcf,'Units','normalized');
    set(gcf,'Position',[ 0,         0  ,  0.5363  ,  0.9143]);
    RES=DIF_simulated_once(sigma_vec,QQ);

    todoS{II}.RES=RES;


end

%%
% PLOT FIGURES

figure(200);clf;
for II=1:length(todoS)

    

     figure(200+II);clf;   
    QQ=todoS{II}.QQ; %load data
    RES=todoS{II}.RES; % load results

    % unpack results:
    pR=RES.pR;
    pU=RES.pU;
    pF=RES.pF;
    Z=RES.Z;
    xx=RES.xx;
    yy=RES.yy;
    N1=RES.N1;
    N2=RES.N2;
    M=RES.M;
    p_t_tp_s=RES.p_t_tp_s;%denoiser conditional at each stage
    q_tp_s=RES.q_tp_s;% forward process marginal that gets noisier
    q_tp_t_s=RES.q_tp_t_s;%noiser conditional at each stage
    p_tp_s=RES.p_tp_s; % generatd marginals across T iterations
    ps_t_t_s=RES.ps_t_t_s;% kernel of the sampling process
    pjoint_t_t_s=RES.pjoint_t_t_s; %joint distribution kernel of the sampling process
    T=RES.T;
    sigma_vec=RES.sigma_vec;
    stat=RES.stat; % further stats

      if T<50
        NDISP1=4;
        NDISP2=13;
        
    else
         NDISP1=11;
        NDISP2=11;
      end

    subplot(NDISP1,NDISP2,1+2*NDISP2*(1-1));
    imagesc(xx,yy,reshape(pR,N1,N2),[0 max(pR(:))]);axis xy;axis off;title('x_0');
    axis square;


  
    subplot(NDISP1,NDISP2,NDISP2*1+1+2*NDISP2*(1-1)+T);
      
    imagesc(xx,yy,reshape(pF,N1,N2),[0 max(pR(:))]);axis xy;axis off;title('x_T');
    axis square;

    for t=1:T
        subplot(NDISP1,NDISP2,t+1+2*NDISP2*(1-1));
        imagesc(xx,yy,reshape(q_tp_s{t},N1,N2),[0 max(pR(:))]);axis xy;title(sprintf('x_{%d}',t));
        axis off
        axis square
    end

    for t=1:T
        subplot(NDISP1,NDISP2,NDISP2*1+t+2*NDISP2*(1-1));
        imagesc(xx,yy,reshape(p_tp_s{t},N1,N2),[0 max(pR(:))]);axis xy;title(sprintf('x_{%d}',t));axis off
        axis square
    end

end
%%

figure(300);clf;

subplot(1,3,3);
legs=cell(size(todoS)); %adding legend
for II=1:length(todoS)
    legs{II}=todoS{II}.name;
end

fprintf('\nReport result of stationarity\n')
for II=1:length(todoS)
    
    QQ=todoS{II}.QQ; %load data
    RES=todoS{II}.RES; % load results
    stat=RES.stat; % further stats
    
    plot(II,stat.my_int_score,'o','MarkerFaceColor','k','MarkerSize',15);hold all;
    fprintf('%s score=%3.3g\n',todoS{II}.name,stat.my_int_score);
end
title('Stationarity')
set(gca,'Xtick',1:length(todoS))
set(gca,'XtickLabels',legs)
ylabel('Stationarity score')
%xlabel('Noise type')
set(gca,'FontSize',14)


subplot(1,3,2);

fprintf('\nReport result of complexity\n')
for II=1:length(todoS)
    
    QQ=todoS{II}.QQ; %load data
    RES=todoS{II}.RES; % load results
    stat=RES.stat; % further stats
    
    plot(II,max(stat.mdkl_dif_vec),'o','MarkerFaceColor','k','MarkerSize',15);hold all;
    fprintf('%s score=%3.3g\n',todoS{II}.name,max(stat.mdkl_dif_vec));
end
title('Complexity')
set(gca,'Xtick',1:length(todoS))
set(gca,'XtickLabels',legs)
ylabel('Complexity (bits)')
%xlabel('Noise type')
set(gca,'FontSize',14)


subplot(1,3,1);

fprintf('\nReport result of Performance (error)\n')
for II=1:length(todoS)
    
    QQ=todoS{II}.QQ; %load data
    RES=todoS{II}.RES; % load results
    stat=RES.stat; % further stats
    
    plot(II,stat.mdkl_score,'o','MarkerFaceColor','k','MarkerSize',15);hold all;
    fprintf('%s score=%3.3g\n',todoS{II}.name,stat.mdkl_score);
end
title('Performance (error)')
set(gca,'Xtick',1:length(todoS))
set(gca,'XtickLabels',legs)
ylabel('Distance to distribution (bits)')
%xlabel('Noise type')
set(gca,'FontSize',14)



%%
% plot(res_all(:,1),res_all(:,3),'k-+','MarkerFaceColor','y','MarkerSize',8, 'LineWidth',2);hold on; set(gca,'FontSize',14);
% xlabel('Number of steps');ylabel('Distance to distribution (bits)');title('Performance (error)');
% 
% subplot(4,1,3);
% 
% plot(res_all(:,1),res_all(:,2),'k-+','MarkerFaceColor','y','MarkerSize',8, 'LineWidth',2);hold on;set(gca,'FontSize',14)
% xlabel('Number of steps');ylabel('Complexity (bits)');title('Distance between Number of steps (complexity) ');
% 
% 
% res_all=[res_all; T,max(stat.mdkl_dif_vec),stat.mdkl_score,stat.my_int_score];