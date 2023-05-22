% Create figure 1 - also create figure S1 (part of the figure the involves
% the stationarity measure)

clear all;close all;clc; % "To begin at the beginning/ Under Milk Wood - Dylan Thomas"


%%
JMP=0.04; %  This defines the grid resolution (smaller number higher resolution)
%JMP=0.025; % high resolution - we use low resolution because of speed

NOISE_TYPE_GAUSSIAN=1; NOISE_TYPE_STRECH=2;NOISE_TYPE_TRIMODAL=3;NOISE_TYPE_FIXED=4;NOISE_TYPE_FIXED_UNIFORM=5; % noise types

TS=[2,3,6,12,24,48,100]; % possible time steps

QQ0=[];
QQ0.JMP=JMP;
QQ0.IS_PLOT=true;
QQ0.IS_SWISS_ROLL=true;
QQ0.NOISE_TYPE=NOISE_TYPE_GAUSSIAN;


todoS=cell(1,1);
tcnt=0;

for II=1:length(TS) % make list of all conditions based on time steps
    tcnt=tcnt+1;
    QQ=QQ0;
    T=TS(II);
    sigma0=0.01;sigma1=0.05; sigma_vec=linspace(sigma0,sigma1,T);
    mclr=0.9*mod([II*123213,(II+123)*12312,II*II+II+123+ II*1232131],256)/256+ [0.1,0.1,0.1];
    todoS{tcnt,1}.QQ=QQ;
    todoS{tcnt,1}.sigma_vec=sigma_vec;
    todoS{tcnt,1}.T=T;
    todoS{tcnt,1}.mclr=mclr;
end

% RUN SIMULATIONS:
for II=1:length(todoS)
    QQ=todoS{II}.QQ;
    sigma_vec=todoS{II}.sigma_vec;
    T=todoS{II}.T;
    fprintf('now in todo %d of %d\n',II,length(todoS));

    figure(101+II);clf;
    set(gcf,'Units','normalized');
    set(gcf,'Position',[ 0,         0  ,  0.5363  ,  0.9143]);
    RES=DIF_simulated_once(sigma_vec,QQ); % Run simulation.

    todoS{II}.RES=RES;
end
%%

% PLOT RESULTS
NDISP1=3; % number of display subplots
NDISP2=3;

figure(300);clf;

res_all=[]; % aggregated performance measures

for II=1:length(todoS)
    QQ=todoS{II}.QQ;
    sigma_vec=todoS{II}.sigma_vec;
    T=todoS{II}.T;

    RES=todoS{II}.RES;
    pR=RES.pR;
    pU=RES.pU;
    Z=RES.Z;
    xx=RES.xx;
    yy=RES.yy;
    N1=RES.N1;
    N2=RES.N2;
    M=RES.M;

    p_t_tp_s=RES.p_t_tp_s;%denoiser conditional at each stage
    q_tp_s=RES.q_tp_s;% forward process marginal that gets noisier
    q_tp_t_s=RES.q_tp_t_s;%noiser conditional at each stage
    p_tp_s=RES.p_tp_s; % generatd marginals across T Number of steps
    ps_t_t_s=RES.ps_t_t_s;% kernel of the sampling process
    pjoint_t_t_s=RES.pjoint_t_t_s; %joint distribution kernel of the sampling process
    T=RES.T;
    sigma_vec=RES.sigma_vec;
    stat=RES.stat; % further stats

    res_all=[res_all; T,max(stat.mdkl_dif_vec),stat.mdkl_score,stat.my_int_score];

    figure(300);

    subplot(NDISP1,NDISP2,1);

    plot(1:T,sigma_vec,'o-','MarkerFaceColor','y','MarkerSize',12, 'LineWidth',2); hold all;title('noise schedule');
    plot([1 T],[min(sigma_vec),min(sigma_vec)],'g--'); hold on;
    plot([1 T],[max(sigma_vec),max(sigma_vec)],'g--'); hold on;
    xlabel('Iteration')


    subplot(NDISP1,NDISP2,2);
    plot(1:T,stat.H_t,'o-','MarkerFaceColor','y','MarkerSize',12, 'LineWidth',2); hold on;title('cond entropy');
    Hmax=-sum(pU.*log2(pU+eps));
    xlabel('Iteration')

    hold on; plot(TS,ones(size(TS))*Hmax,'g--')
    subplot(NDISP1,NDISP2,3);
    plot(1:T,stat.I_t,'o-','MarkerFaceColor','y','MarkerSize',12, 'LineWidth',2); hold on;title('mutual info');
    Imax=-sum(pU.*log2(pU+eps));
    hold on; plot(TS,ones(size(TS))*Imax,'g--')
    xlabel('Iteration')

    subplot(NDISP1,NDISP2,4);
    J_t=stat.mjsd_dif_vec;
    plot(1:(T-1),J_t,'o-','MarkerFaceColor','y','MarkerSize',12, 'LineWidth',2); hold on;title('JSD Number of steps');
    xlabel('Iteration')


    subplot(NDISP1,NDISP2,5);
    D_t=stat.mdkl_dif_vec;
    plot(1:(T-1),D_t,'o-','MarkerFaceColor','y','MarkerSize',12, 'LineWidth',2); hold on;title('DKL Number of steps');
    xlabel('Iteration')

    subplot(NDISP1,NDISP2,6);
    DS_t=stat.mdkl_score;
    plot(T,DS_t,'o-','MarkerFaceColor','y','MarkerSize',12, 'LineWidth',2); hold on;title('DKL score');
    xlabel('Iteration')

    subplot(NDISP1,NDISP2,7);
    DS_t=stat.mdkl_score;
    D_t_max=max(stat.mdkl_dif_vec);subplot(NDISP1,NDISP2,7);
    DS_t=stat.mdkl_score;
    D_t_max=max(stat.mdkl_dif_vec);
    plot(D_t_max,DS_t,'o-','MarkerFaceColor','y','MarkerSize',8, 'LineWidth',2); hold on;title('Performance tradeoff');
    xlabel('Complexity');ylabel('Performance');
    plot(D_t_max,DS_t,'o-','MarkerFaceColor','y','MarkerSize',8, 'LineWidth',2); hold on;title('Performance tradeoff');
    xlabel('Complexity');ylabel('Performance');

    subplot(NDISP1,NDISP2,8);
    plot(T,stat.my_int_score,'o-','MarkerFaceColor','y','MarkerSize',8, 'LineWidth',2); hold on;title('Int score');
    xlabel('Steps');ylabel('Int score');

end

%%
figure(400);clf;
plot(res_all(:,2),res_all(:,3),'k-+','MarkerFaceColor','y','MarkerSize',8, 'LineWidth',2);hold on;
for ll=1:size(res_all,1)
    text(res_all(ll,2),res_all(ll,3),sprintf('    T=%d',res_all(ll,1)),'HorizontalAlignment','left');hold on;
end

figure(401);clf; % Figure 3

subplot(3,1,2);
plot(res_all(:,1),res_all(:,3),'k-+','MarkerFaceColor','y','MarkerSize',8, 'LineWidth',2);hold on; set(gca,'FontSize',14);
xlabel('Number of steps');ylabel('Distance to distribution (bits)');title('Performance (error)');

subplot(3,1,3);

plot(res_all(:,1),res_all(:,2),'k-+','MarkerFaceColor','y','MarkerSize',8, 'LineWidth',2);hold on;set(gca,'FontSize',14)
xlabel('Number of steps');ylabel('Complexity (bits)');title('Distance between Number of steps (complexity) ');


for II=1:length(todoS)
    QQ=todoS{II}.QQ; % get params
    mclr=todoS{II}.mclr; % color of condition
    RES=todoS{II}.RES; % get results
    
    % Un pack results:
    pR=RES.pR;
    pU=RES.pU;
    Z=RES.Z;
    xx=RES.xx;
    yy=RES.yy;
    N1=RES.N1;
    N2=RES.N2;
    M=RES.M;
    p_t_tp_s=RES.p_t_tp_s;%denoiser conditional at each stage
    q_tp_s=RES.q_tp_s;% forward process marginal that gets noisier
    q_tp_t_s=RES.q_tp_t_s;%noiser conditional at each stage
    p_tp_s=RES.p_tp_s; % generatd marginals across T Number of steps
    ps_t_t_s=RES.ps_t_t_s;% kernel of the sampling process
    pjoint_t_t_s=RES.pjoint_t_t_s; %joint distribution kernel of the sampling process
    T=RES.T;
    sigma_vec=RES.sigma_vec;
    stat=RES.stat; % further stats

    % PLOT:
    subplot(3,1,1)
    plot(1:T,sigma_vec,'o-','MarkerFaceColor',mclr,'MarkerSize',8, 'LineWidth',2,'Color',mclr); hold all;title('Noise schedule');

    xlabel('Steps');set(gca,'FontSize',14);ylabel('Noise (sigma)')
    set(gcf,'Units','Normalized')
    set(gcf,'Position',[0.2137    0.3080    0.5084    0.5857])

    subplot(3,1,2);
    plot(res_all(II,1),res_all(II,3),'o','MarkerFaceColor',mclr,'MarkerSize',10, 'LineWidth',2,'Color',mclr);hold on;
    subplot(3,1,3);title('Complexity')
    plot(res_all(II,1),res_all(II,2),'o','MarkerFaceColor',mclr,'MarkerSize',10, 'LineWidth',2,'Color',mclr);hold on;

end

% plot min and max:
subplot(3,1,1)
plot([1 T],[min(sigma_vec),min(sigma_vec)],'g--'); hold on;
plot([1 T],[max(sigma_vec),max(sigma_vec)],'g--'); hold on;

% make legend
mleg=cell(size(TS));
for ll=1:length(TS)
    mleg{ll}=sprintf('T=%d',TS(ll));
end
legend(mleg,'AutoUpdate','off')
%%
figure(402);clf; % Figure 3

subplot(4,1,2);
plot(res_all(:,1),res_all(:,3),'k-+','MarkerFaceColor','y','MarkerSize',8, 'LineWidth',2);hold on; set(gca,'FontSize',14);
xlabel('Number of steps');ylabel('Distance to distribution (bits)');title('Performance (error)');

subplot(4,1,3);

plot(res_all(:,1),res_all(:,2),'k-+','MarkerFaceColor','y','MarkerSize',8, 'LineWidth',2);hold on;set(gca,'FontSize',14)
xlabel('Number of steps');ylabel('Complexity (bits)');title('Distance between Number of steps (complexity) ');

subplot(4,1,4);

plot(res_all(:,1),res_all(:,4),'k-+','MarkerFaceColor','y','MarkerSize',8, 'LineWidth',2);hold on;set(gca,'FontSize',14)
xlabel('Number of steps');ylabel('Stationarity score');title('Stationarity');


for II=1:length(todoS)
    QQ=todoS{II}.QQ; % get params
    mclr=todoS{II}.mclr; % color of condition
    RES=todoS{II}.RES; % get results
    
    % Un pack results:
    pR=RES.pR;
    pU=RES.pU;
    Z=RES.Z;
    xx=RES.xx;
    yy=RES.yy;
    N1=RES.N1;
    N2=RES.N2;
    M=RES.M;
    p_t_tp_s=RES.p_t_tp_s;%denoiser conditional at each stage
    q_tp_s=RES.q_tp_s;% forward process marginal that gets noisier
    q_tp_t_s=RES.q_tp_t_s;%noiser conditional at each stage
    p_tp_s=RES.p_tp_s; % generatd marginals across T Number of steps
    ps_t_t_s=RES.ps_t_t_s;% kernel of the sampling process
    pjoint_t_t_s=RES.pjoint_t_t_s; %joint distribution kernel of the sampling process
    T=RES.T;
    sigma_vec=RES.sigma_vec;
    stat=RES.stat; % further stats

    % PLOT:
    subplot(4,1,1)
    plot(1:T,sigma_vec,'o-','MarkerFaceColor',mclr,'MarkerSize',8, 'LineWidth',2,'Color',mclr); hold all;title('Noise schedule');

    xlabel('Steps');set(gca,'FontSize',14);ylabel('Noise (sigma)')
    set(gcf,'Units','Normalized')
    set(gcf,'Position',[0.2137    0.3080    0.5084    0.5857])

    subplot(4,1,2);
    plot(res_all(II,1),res_all(II,3),'o','MarkerFaceColor',mclr,'MarkerSize',10, 'LineWidth',2,'Color',mclr);hold on;
    subplot(4,1,3);title('Complexity')
    plot(res_all(II,1),res_all(II,2),'o','MarkerFaceColor',mclr,'MarkerSize',10, 'LineWidth',2,'Color',mclr);hold on;
     subplot(4,1,4);title('Stationarity')
    plot(res_all(II,1),res_all(II,4),'o','MarkerFaceColor',mclr,'MarkerSize',10, 'LineWidth',2,'Color',mclr);hold on;

end

% plot min and max:
subplot(4,1,1)
plot([1 T],[min(sigma_vec),min(sigma_vec)],'g--'); hold on;
plot([1 T],[max(sigma_vec),max(sigma_vec)],'g--'); hold on;

% make legend
mleg=cell(size(TS));
for ll=1:length(TS)
    mleg{ll}=sprintf('T=%d',TS(ll));
end
legend(mleg,'AutoUpdate','off')
