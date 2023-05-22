
clear all;close all;clc;
addpath('~/ResearchMIT/toolboxes/nUTIL');


%%
JMP=0.04;
JMP=0.025; %  This defines the grid resolution (smaller number higher resolution)
NOISE_TYPE_GAUSSIAN=1; NOISE_TYPE_STRECH=2;NOISE_TYPE_TRIMODAL=3;NOISE_TYPE_FIXED=4;NOISE_TYPE_FIXED_UNIFORM=5;


% DECIDE WHAT TO DO:
QQ0=[]; % parameters shared between coditions.
QQ0.JMP=JMP;
QQ0.IS_PLOT=true;
QQ0.IS_SWISS_ROLL=true;


todoS=cell(1,1); %initialize conditions
tcnt=0;

tcnt=tcnt+1; % 2D: fixed noise uniform condition 
QQ=QQ0;
QQ.NOISE_TYPE=NOISE_TYPE_FIXED_UNIFORM;
T=10;
sigma0=0.01;sigma1=1.0; sigma_vec=linspace(sigma0,sigma1,T);
todoS{tcnt,1}.QQ=QQ;
todoS{tcnt,1}.sigma_vec=sigma_vec;
todoS{tcnt,1}.T=T;

tcnt=tcnt+1; %2E: fixed noise non uniform condition
QQ=QQ0;
QQ.NOISE_TYPE=NOISE_TYPE_FIXED;
T=10;
sigma0=0.01;sigma1=1.0; sigma_vec=linspace(sigma0,sigma1,T);
todoS{tcnt,1}.QQ=QQ;
todoS{tcnt,1}.sigma_vec=sigma_vec;
todoS{tcnt,1}.T=T;


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
NDISP1=4;
NDISP2=13;
figure(200);clf;
for II=1:length(todoS)

        
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

    subplot(NDISP1,NDISP2,1+2*NDISP2*(II-1));
    imagesc(xx,yy,reshape(pR,N1,N2),[0 max(pR(:))]);axis xy;axis off;title('x_0');
    axis square;


    subplot(NDISP1,NDISP2,NDISP2*1+1+2*NDISP2*(II-1)+T);
    imagesc(xx,yy,reshape(pF,N1,N2),[0 max(pR(:))]);axis xy;axis off;title('x_T');
    axis square;

    for t=1:T
        subplot(NDISP1,NDISP2,t+1+2*NDISP2*(II-1));
        imagesc(xx,yy,reshape(q_tp_s{t},N1,N2),[0 max(pR(:))]);axis xy;title(sprintf('x_{%d}',t));
        axis off
        axis square
    end

    for t=1:T
        subplot(NDISP1,NDISP2,NDISP2*1+t+2*NDISP2*(II-1));
        imagesc(xx,yy,reshape(p_tp_s{t},N1,N2),[0 max(pR(:))]);axis xy;title(sprintf('x_{%d}',t));axis off
        axis square
    end

end
%%

% MAKING TITLES FOR FIGURE 2 (IGNORE IMAGESC, WE SIMLY USE THE TITLE)
% figure
% subplot(2,2,1);
% imagesc(xx,yy,reshape(pR,N1,N2));title('q(x_{t}|x_{t-1})')
% set(gca,'FontSize',14)
% subplot(2,2,2);
% imagesc(xx,yy,reshape(pR,N1,N2));title('p_s(x_{t-1}|x_{t})')
% set(gca,'FontSize',14)
