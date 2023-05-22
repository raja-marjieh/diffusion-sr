%%% Example of simulation run

JMP=0.025; %  This defines the grid resolution (smaller number higher resolution)

% possible noise type to select from:
NOISE_TYPE_GAUSSIAN=1; NOISE_TYPE_STRECH=2;NOISE_TYPE_TRIMODAL=3;NOISE_TYPE_FIXED=4;NOISE_TYPE_FIXED_UNIFORM=5;


T=10; % number of steps
sigma0=0.03;sigma1=0.1; sigma_vec=linspace(sigma0,sigma1,T);

% set parameters:
QQ0=[];
QQ.JMP=JMP; % resolution
QQ.IS_PLOT=true; % plot intrmidiate results
QQ.IS_SWISS_ROLL=true; % Swiss roll initial distribution
QQ.NOISE_TYPE=NOISE_TYPE_GAUSSIAN; % chose noise type (here Gaussian)

RES=DIF_simulated_once(sigma_vec,QQ); % run sumulation and plot it

return % we don't want to do anything automatically but you can try to run each of the cells below
%%

%CREATE  FIGURES OF THE PAPER:
%%
% figure 2 panels A-C
DIF_figure2;
%%
% figure 2 panels D-E
DIF_figure2b;
%%
% figure 3
DIF_figure3;
%The stationarity measure is reported also in figure S1
%%
% figure 4 results panel summary
DIF_figure4_panel;

%%

%%%% Supplumental figures:
%%
DIF_figure_S2 % Comparing with sampling with noise injection and without noise injection
%%
DIF_figure_S3 % Bad noise display
%%
DIF_figure_S4 % Bad noise stats:  running all conditions and plotting also the stationarity score
%%
DIF_figure_S5 % Bimodal noise