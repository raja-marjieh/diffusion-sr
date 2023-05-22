function RES=DIF_simulated_once_v03(sigma_vec,QQ)


%%% Read parameters:
T=length(sigma_vec);

if isfield(QQ,'JMP')
    JMP=QQ.JMP;
else
    JMP=0.04;
end

if isfield(QQ,'IS_PLOT')
    IS_PLOT=QQ.IS_PLOT;
else
    IS_PLOT=True;
end

if isfield(QQ,'IS_SWISS_ROLL')
    IS_SWISS_ROLL=QQ.IS_SWISS_ROLL;
else
    IS_SWISS_ROLL=True;
end

NOISE_TYPE_GAUSSIAN=1; NOISE_TYPE_STRECH=2;NOISE_TYPE_TRIMODAL=3;NOISE_TYPE_FIXED=4;NOISE_TYPE_FIXED_UNIFORM=5;

if isfield(QQ,'NOISE_TYPE')
    NOISE_TYPE=QQ.NOISE_TYPE;
else
    NOISE_TYPE=NOISE_TYPE_GAUSSIAN;
end

if isfield(QQ,'PROIDE_pR')
    PROIDE_pR=QQ.PROIDE_pR;
else
    PROIDE_pR=[];
end

if isfield(QQ,'IS_DENOISE_ONLY')
    IS_DENOISE_ONLY=QQ.IS_DENOISE_ONLY;
else
    IS_DENOISE_ONLY=false;
end
%%

%%% Set parameters:
xx=-0.5:JMP:0.5; %% x grid
yy=-0.5:JMP:0.5; %% y griz
[X,Y] = meshgrid(xx,yy); % grids as a matrices
Z=[X(:), Y(:)];  % vector of all coordinate locations ordered as (very) long  columns (x|,y|) for the x and y coordinates
N1=length(yy); %% length of the y-grid
N2=length(xx);%% length of the x-grid
M=length(Z); %% length of the z-grid


%%% compute elegant number of displays subplots (display purpose only):
T_MORE=5;
NDISP1=floor(sqrt(T+T_MORE));
NDISP2=floor(sqrt(T+T_MORE));
if NDISP1*NDISP2<(T+T_MORE)
    NDISP1=NDISP1+1;
end
if NDISP1*NDISP2<(T+T_MORE)
    NDISP2=NDISP2+1;
end

%%

pU=ones(size(Z,1),1);pU=pU/sum(pU(:)); % uniform distribution


comb=0.1; % The proportion of uniform noise in the prior (this makes learning a bit more challanging)


% A distribution with 3 modes we sometime use:
W=8*(0.04).^2; % relative width
x1=[0.1 0.2]; % mode (1)
x2=[-0.2 -0.1];% mode (2)
x3=[0.2 -0.2];% mode (3)
pR3=mvnpdf(Z, x1, eye(2)*W) + mvnpdf(Z, x2, eye(2)*W)+  mvnpdf(Z, x3, [2 1; 1 2]*W); % the distribution(arbitrary)
pR3=pR3/sum(pR3(:)); % normalization

if IS_SWISS_ROLL % compute Swiss roll
    roll_width=0.04;
    t=3*pi/2 * (1 + 2*(0:0.01:1));
    ts=[t.*sin(t);-t.*cos(t)]';
    ts(:,1)=0.8*(ts(:,1)-min(ts(:,1)))/(max(ts(:,1))-min(ts(:,1)))-0.4;
    ts(:,2)=0.8*(ts(:,2)-min(ts(:,2)))/(max(ts(:,2))-min(ts(:,2)))-0.4;

    % make the roll with some width
    min_dist=999*ones(size(X));
    for ll=1:length(t)
        pnt=ts(ll,:);
        min_dist=min(min_dist,sqrt( ((X-pnt(1)).^2) +  ((Y-pnt(2)).^2)) );

    end
    pR=min_dist<roll_width;
    pR=pR/sum(pR(:)); % normalization to a distribution

    pR=pR(:);
    pR=(1-comb)*pR+(comb)*pU;   % Now add a bit of uniform distribution (this makes learning a bit more challanging)

else
    pR=pR3; % start with 3 modes distribution
    pR=pR/sum(pR(:)); % normalize distribution (just to make sure)
    pR=(1-comb)*pR+(comb)*pU;   % Now add a bit of uniform distribution (this makes learning a bit more challanging)

end

if ~isempty(PROIDE_pR) % use predetermined distribution.
    pR=PROIDE_pR;
end

pR=pR/sum(pR(:)); % normalize distribution (just to make sure)

% pR is the intial distribution
% pF is the final distribution
switch NOISE_TYPE
    case NOISE_TYPE_FIXED
        pF=pR3;
    otherwise
        pF=pU;
end
% %%%%%%%%%%%%%%%%%%
% %%% PLOT initial distribution (debugging)
%
%     figure(1);clf;
%     subplot(NDISP1,NDISP2,1);
%     imagesc(xx,yy,reshape(pR,N1,N2),[0 max(pR(:))]);axis xy;title('Intial distribution');
% %%%%
% %%%%%%%%%%%%%%%%%%

%%
p_t_tp_s=cell(T,1); %denoiser conditional at each stage
q_tp_s=cell(T,1); % forward process marginal that gets noisier
q_tp_t_s=cell(T,1); %noiser conditional at each stage

fprintf('forward pass\t ...')
for t=1:(T)
    fprintf('.');
    if t==1
        q_t=pR;

    else
        q_t=q_tp_s{t-1};
    end

    sigma_iter=sigma_vec(t);

    SIGMA_iter=eye(2)*(sigma_iter^2);
    SIGMA_iter_inv=inv(SIGMA_iter);
    pNoise=nan(M,M);
    for I=1:M
        z=Z(I,:);
        %Q=mvnpdf(Z,z,SIGMA_iter); % we can't use it becuase of the
        %circular boundary conditions

        % circular distance
        min_Zmu=999*ones(size(Z));
        for xshift=-1:1
            for yshift=-1:1
                Zmu=(Z-repmat(z+[xshift,yshift],M,1));
                pos=sum(Zmu.*Zmu,2)<sum(min_Zmu.*min_Zmu,2);
                min_Zmu(pos,:)=Zmu(pos,:);
            end
        end
        Zmu=min_Zmu;

        switch NOISE_TYPE
            case NOISE_TYPE_GAUSSIAN
                Q=exp( -0.5*sum((Zmu*SIGMA_iter_inv).*Zmu,2) ); % Gaussian noise on circular boundaries

            case NOISE_TYPE_STRECH % bimodal noise
                SHIFT=[0.05,0.05];
                Zmu1=Zmu+repmat( SHIFT,M,1);
                Zmu2=Zmu+repmat(-SHIFT,M,1);

                Q1=exp( -0.5*sum((Zmu1*SIGMA_iter_inv).*Zmu1,2) );
                Q2=exp( -0.5*sum((Zmu2*SIGMA_iter_inv).*Zmu2,2) );
                Q1=Q1/sum(Q1);
                Q2=Q2/sum(Q2);
                Q=Q1+Q2;
            case NOISE_TYPE_TRIMODAL % trimodal noise
                SHIFT=[-0.05,-0.05];
                Zmu1=Zmu+repmat( SHIFT,M,1);
                Zmu2=Zmu+repmat(-SHIFT,M,1);

                Q1=exp( -0.5*sum((Zmu1*SIGMA_iter_inv).*Zmu1,2) );
                Q2=exp( -0.5*sum((Zmu2*SIGMA_iter_inv).*Zmu2,2) );
                Q3=exp( -0.5*sum((Zmu*SIGMA_iter_inv).*Zmu,2) );
                Q1=Q1/sum(Q1);
                Q2=Q2/sum(Q2);
                Q3=Q3/sum(Q3);
                Q=0.5*Q1+0.5*Q2+Q3;
            case {NOISE_TYPE_FIXED,NOISE_TYPE_FIXED_UNIFORM} % fade-in noise
                SIGMA_delta=eye(2)*(0.01^2);
                SIGMA_delta_inv=inv(SIGMA_delta);
                Q1=exp( -0.5*sum((Zmu*SIGMA_delta_inv).*Zmu,2) );
                Q2=pF;
                Q=(1-sigma_iter)*Q1+sigma_iter*Q2;

            otherwise
                assert(1==0);
        end

        Q=Q./sum(Q(:)); % normalize in any case;
        pNoise(:,I)=Q;
    end
    q_tp_t=pNoise; % iteration noising kernel

    %%% perform Bayesian inference with noise as likelihood
    pDen=nan(M,M);
    for I=1:M
        pDen(:,I)=q_tp_t(I,:).*(q_t');
        pDen(:,I)=pDen(:,I)/sum(pDen(:,I));
    end
    p_t_tp=pDen;

    % compute forward iteration:
    q_tp=q_tp_t*q_t;

    % save data:
    p_t_tp_s{t}=p_t_tp;
    q_tp_s{t}=q_tp;
    q_tp_t_s{t}=q_tp_t;

end
fprintf('\n')

p_tp_s=cell(T,1); % generatd marginals across T iterations
ps_t_t_s=cell(T,1); % kernel of the sampling process
pjoint_t_t_s=cell(T,1); %joint distribution kernel of the sampling process
fprintf('backward pass\t ...')
for t=T:-1:1
    fprintf('.');
    if t==T
        p_t=pF;

    else
        %p_t=q_tp_s{t+1};
        p_t=p_tp_s{t+1};
        
    end

    if IS_DENOISE_ONLY
        my_kernel=p_t_tp_s{t};
        p_tp=p_t_tp_s{t}*p_t; % this use my kernel but for efficncy we compute it in another way.
        p_tp_s{t}=p_tp;
        ps_t_t_s{t}=my_kernel;

    else
        my_kernel=p_t_tp_s{t}*q_tp_t_s{t};
        p_tp=p_t_tp_s{t}*(q_tp_t_s{t}*p_t); % this use my kernel but for efficncy we compute it in another way.
        %pjoint_t_t_s{t}=q_tp_t_s{t}.*repmat(p_t',M,1);
        p_tp_s{t}=p_tp;
        ps_t_t_s{t}=my_kernel;
    end
end
fprintf('\n');

%%
fprintf('computing stats...\n')

%%
% Compute stats:
% q are the generated noise marginals
%
t=T;
vec=p_tp_s{t}./(p_tp_s{t-1}+eps);

for t=(T-1):-1:1
    vec2=sum(ps_t_t_s{t}.*repmat(vec,1,M));
    vec=vec2';
end
vec_int_dif=vec;
my_int_score=norm(vec_int_dif-ones(size(vec_int_dif)))/sqrt(length(vec_int_dif));

%%
%%% make performance scores.

mdkl_dif_vec=nan(T-1,1); % differenc between consecutive iterations in terms of DKL
mjsd_dif_vec=nan(T-1,1); % differenc between consecutive iterations in terms of JSD

for t=(T-1):-1:1
    mdkl_dif_vec(t)=DKL2(q_tp_s{t},q_tp_s{t+1});
    mjsd_dif_vec(t)=JSD2(q_tp_s{t},q_tp_s{t+1});

end

mdkl_score=DKL2(p_tp_s{1},pR);
mjsd_score=JSD2(p_tp_s{1},pR);

%%
% %%% cond  entropy and mutual info
if ~IS_DENOISE_ONLY
H_t=nan(T,1);
I_t=nan(T,1);
for t=T:-1:1
    if t==T
        p_t=pF;

    else
        %p_t=q_tp_s{t+1}; 
        p_t=p_tp_s{t+1};
    end

    p_tp=p_t_tp_s{t}*(q_tp_t_s{t}*p_t); % this use my kernel but for efficncy we compute it in another way.
    pjoint_t_t_s{t}=q_tp_t_s{t}.*repmat(p_t',M,1);
    %p_tp_s{t}=p_tp;
    H_t(t)=-sum((q_tp_t_s{t}.*log2(q_tp_t_s{t}+eps)),1)*p_t;
    I_t(t)=sum(sum(   (pjoint_t_t_s{t}.*log2 ( eps+ pjoint_t_t_s{t}./  (eps+repmat(p_t',M,1).*repmat(p_tp,1,M)) )   )     ));

end
end
%%

fprintf('saving results...\n')
RES.pR=pR;
RES.pU=pU;
RES.pF=pF;
RES.Z=Z;
RES.xx=xx;
RES.yy=yy;
RES.X=X;
RES.Y=Y;
RES.N1=N1;
RES.N2=N2;
RES.M=M;

RES.p_t_tp_s=p_t_tp_s;%denoiser conditional at each stage
RES.q_tp_s=q_tp_s;% forward process marginal that gets noisier
RES.q_tp_t_s=q_tp_t_s;%noiser conditional at each stage
RES.p_tp_s=p_tp_s; % generatd marginals across T iterations
RES.ps_t_t_s=ps_t_t_s;% kernel of the sampling process
RES.pjoint_t_t_s=pjoint_t_t_s; %joint distribution kernel of the sampling process
RES.T=T;
RES.sigma_vec=sigma_vec;
RES.stat.my_int_score=my_int_score;
RES.stat.vec_int_dif=vec_int_dif;% intergral dif value as vector
RES.stat.mdkl_score=mdkl_score; % DKL-generated vs true distribution
RES.stat.mjsd_score=mjsd_score; % JSD-generated vs true distribution
RES.stat.mdkl_dif_vec=mdkl_dif_vec; % DKL-difference between consecutive iterations
RES.stat.mjsd_dif_vec=mjsd_dif_vec; % JSD-difference between consecutive iterations
if ~IS_DENOISE_ONLY
    RES.stat.H_t=H_t;
    RES.stat.I_t=I_t;
end


%%
%%%%%%%%%%%%%%%%%%
%%% PLOT results:

if IS_PLOT
    fprintf('plotting results...\n')
    %%%% noising process
    subplot(2*NDISP1,NDISP2,1);
    imagesc(xx,yy,reshape(pR,N1,N2),[0 max(pR(:))]);axis xy;axis off;title('Intial Distribution Noising');

    for t=1:T
        subplot(2*NDISP1,NDISP2,2+t);
        imagesc(xx,yy,reshape(q_tp_s{t},N1,N2),[0 max(pR(:))]);axis xy;title(sprintf('noising iter=%d',t));
        axis off
    end

    subplot(2*NDISP1,NDISP2,NDISP1*NDISP2);
    imagesc(xx,yy,reshape(q_tp_s{end},N1,N2),[0 max(pR(:))]);axis xy;title('noising last iteration');
    axis off

    for t=1:T
        subplot(2*NDISP1,NDISP2,NDISP1*NDISP2+2+t);
        mdkl=DKL2(p_tp_s{t},pR);
        mjsd=JSD2(p_tp_s{t},pR);

        imagesc(xx,yy,reshape(p_tp_s{t},N1,N2),[0 max(pR(:))]);axis xy;title(sprintf('denoising iter=%d',t));axis off
    end

    subplot(2*NDISP1,NDISP2,NDISP1*NDISP2+NDISP1*NDISP2);
    imagesc(xx,yy,reshape(pU,N1,N2),[0 max(pR(:))]);axis xy;title('denoising first iteration');axis off

    subplot(2*NDISP1,NDISP2,NDISP1*NDISP2+1);
    imagesc(xx,yy,reshape(p_tp_s{1},N1,N2),[0 max(pR(:))]);axis xy;title(sprintf('gen-distrib. DKL/JSD= %3.5g / %3.5 ',mdkl_score,mjsd_score));axis off


    subplot(2*NDISP1,NDISP2,NDISP1*NDISP2+2);
    imagesc(xx,yy,reshape(pR,N1,N2),[0 max(pR(:))]);axis xy;title('comparison initial distribution');axis off

    ttt=1:T;
    subplot(2*NDISP1,NDISP2,NDISP1*NDISP2+NDISP1*NDISP2-1);plot(sigma_vec,'+-');hold on;
    sigma0=min(sigma_vec);
    sigma1=max(sigma_vec);
    plot(ttt,ones(size(ttt))*sigma0,'g--','LineWidth',3);hold on;
    plot(ttt,ones(size(ttt))*sigma1,'g--','LineWidth',3);hold on;

    subplot(2*NDISP1,NDISP2,NDISP1*NDISP2+NDISP1*NDISP2-2)

    plot(vec_int_dif,'-','LineWidth',1.5);hold on;
    plot([1 length(vec_int_dif)],[1 1 ],'g--','LineWidth',2);


    title(sprintf('int-score: %g',my_int_score));

    subplot(2*NDISP1,NDISP2,NDISP1*NDISP2+NDISP1*NDISP2-3)
    imagesc(xx,yy,reshape(vec_int_dif,N1,N2));axis xy;axis off;
    title("Deviation of integral from 1")
    colorbar
end


