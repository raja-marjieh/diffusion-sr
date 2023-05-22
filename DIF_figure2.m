
clear all;close all;clc; % "To begin at the beginning/ Under Milk Wood - Dylan Thomas"

%%
%JMP=0.04;  % low resolution
JMP=0.025; %  This defines the grid resolution (smaller number higher resolution)
NOISE_TYPE_GAUSSIAN=1; NOISE_TYPE_STRECH=2;NOISE_TYPE_TRIMODAL=3;NOISE_TYPE_FIXED=4;NOISE_TYPE_FIXED_UNIFORM=5; % noise types


% SET PARAMETERS
QQ0=[];  % shared parameters
QQ0.JMP=JMP;
QQ0.IS_PLOT=true;
QQ0.IS_SWISS_ROLL=true;


todoS=cell(1,1); % list of conditions to run
tcnt=0;

tcnt=tcnt+1;
QQ=QQ0;
QQ.NOISE_TYPE=NOISE_TYPE_GAUSSIAN;% Figures 2A-B
T=10;% number of steps
sigma0=0.03;sigma1=0.1; sigma_vec=linspace(sigma0,sigma1,T);
todoS{tcnt,1}.QQ=QQ;
todoS{tcnt,1}.sigma_vec=sigma_vec;
todoS{tcnt,1}.T=T;

tcnt=tcnt+1;
QQ=QQ0;
QQ.NOISE_TYPE=NOISE_TYPE_STRECH;% Streched noise: Figure 2C
T=10;% number of steps
sigma0=0.02;sigma1=0.07;sigma_vec=linspace(sigma0,sigma1,T);
todoS{tcnt,1}.QQ=QQ;
todoS{tcnt,1}.sigma_vec=sigma_vec;
todoS{tcnt,1}.T=T;


% RUN CONDITIONS:
for II=1:length(todoS)
    QQ=todoS{II}.QQ;
    sigma_vec=todoS{II}.sigma_vec;
    T=todoS{II}.T;

    figure(101+II);clf;
    set(gcf,'Units','normalized');
    set(gcf,'Position',[ 0,         0  ,  0.5363  ,  0.9143]);
    RES=DIF_simulated_once(sigma_vec,QQ); % RUN SIMULATION
    todoS{II}.RES=RES;
end
%%


NDISP1=4; % number of display subplots.
NDISP2=13;
figure(200);clf;
for II=1:length(todoS)
    QQ=todoS{II}.QQ; % get parameters for condition
    RES=todoS{II}.RES; % get results for conditions
    % UNPACK RESULTS:
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
    T=RES.T; % number of steps
    sigma_vec=RES.sigma_vec; % noise schedule
    stat=RES.stat; % further stats

    % PLOT STUFF:
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
