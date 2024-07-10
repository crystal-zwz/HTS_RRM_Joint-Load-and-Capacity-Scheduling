% 二维卫星场景2*3
% - BW-POW: Flexible POW and BW allocation. Fixed user assignment to the beams
% - MAP:    Fixed resources per beam.       Flexible beam user assignment.
% - BW-MAP: Flexible BW allocation.         Flexible beam user assignment.
% - CR
% - CR_BW_alternate

%%
close all
clear all
clc
orig_state = warning;
warning off
%% Monte-Carlo parameters
Nsims_v= 1;   
formatOut='dd_mm_yyyy_HH_MM';
text_file=['System_sim_batch' datestr(datetime('now'),formatOut)]; % Data labelling

% Parent folder for storing the numerical results
parentpath=['System simulations '  datestr(datetime('now'),formatOut)];
mkdir(parentpath)

%% Selection of resource allocation strategies
en.SR=1;
en.SR_BW_alternate=1;
en.POW=0; % Boolean variable. Enables ( with value equal to 1) the optimization of flexible power allocation to cope with the traffic demand.  
en.BW=0; % Boolean variable. Enables ( with value equal to 1) the optimization of flexible bandwidth allocation to cope with the traffic demand.  
en.MAP=0; % Boolean variable. Enables ( with value equal to 1) the optimization of flexible beam-user mapping( with fixed resources) to cope with the traffic demand.  
en.BW_MAP=0; % Boolean variable. Enables( with value equal to 1) the joint optimization of bandwidth and beam-user mapping to cope with the traffic demand.  
en.BW_POW=0; % Boolean variable. Enables ( with value equal to 1) the joint optimization of bandwidth and power to cope with the traffic demand.  
 
%%    Two Dimesion scenario
K=6;% Number beams. IT MUST BE EVEN)
%%% Beam modeling- Bessel
R=50; % km,Beam Radius
beam_cross_roll_off=3;  % dB, roll off value at beam-crossover

% Place beams in the one-dimesional scenario
Roll_off_bound=-beam_cross_roll_off;
% b-th波束孔径距n-th用户位置的距离
% 对贝塞尔模型中波束半径为R的波束，在排布局中波束的位置使相邻波束之间的距离为2R
d=1:0.001:2*R;
u=2.07123*d/R;
% 从b-th波束到中心波束内n-th用户的信道增益
a= besselj(1,u)./(2*u)+ 36*besselj(3,u)./(u.^3);
G_dB_ref=10*log10( abs(besselj(1,u)./(2*u)+ 36*besselj(3,u)./(u.^3)).^2  );
% 查找第一个非零值 初始
ind_roll_off=find(G_dB_ref<=Roll_off_bound,1,'first');
d_rollOff=d(ind_roll_off);

% 可以允许的信噪比下降范围
ind_roll_off_f =find(G_dB_ref<=(8.7-15),1,'first');

% 半径
R_0=d(ind_roll_off);
R_f=d(ind_roll_off_f);

Distance_beam=2*d_rollOff;

% 波束中心坐标
Center_Beams=zeros(K,2);
Center_Beams(1:K/2,1)=0;
Center_Beams(1:K/2,2)=Distance_beam*(0:K/2-1);
Center_Beams(K/2+1:K,1)=Distance_beam*sqrt(3)/2;
Center_Beams(K/2+1:K,2)=Distance_beam*(0:K/2-1)+Distance_beam/2;

%% Traffic distributions definition
 
% 同质流量
alpha_v{1}=ones(1,K); % Scenario 1
label_alpha{1}='HT'; % Scenario 1 label
mkdir(parentpath,label_alpha{1}) % Create folder to store the data

% 热点
alpha_v{2}=[5 5 30 5 5 5]; % Scenario 2
label_alpha{2}='HS'; % Scenario 2 label
mkdir(parentpath,label_alpha{2})  % Create folder to store the data

% 宽热点
alpha_v{3}=[4 4 1 1 1 1]; % Scenario 3
label_alpha{3}='WHS';  % Scenario 3 label
mkdir(parentpath,label_alpha{3})  % Create folder to store the data

Nscen=1;
% Nscen=length(alpha_v);

%% System parameters

M=4; % Number of carrier per colour ( Two-colour scheme in the one-dimensional scenario)
Delta_W=  1/(2*M); % Portion of the carrier bandwidth

total_band=500e6; % Total Bandwidth
% Traffic requested per user in Mbps, 
% normalized to the avalaible bandwidth : 25 Mbps/ 500 MHz
Req_user_ref= 25e6/(total_band); 
% log2(1+ SNR_eff_beam) ,log(1+10^1.35)/log(2) with SNR_eff_beam the average SNR for 
% a uniform distribution of user within a beam ( approx 13.5 dB )
SE_av=4.5271; 
% Average spectral effiency with a two-Colour scheme 6.8e3/500/6
SE_beam_uniform=SE_av/2;
max_user_per_beam = round(SE_beam_uniform/Req_user_ref);

% Number of user to achieve the  capacity per beam
Nuser_beam=  SE_beam_uniform./Req_user_ref; 
Nuser_tot=round(K*Nuser_beam); 
% Number of total users to achieve the system capacity. 
% We consider the satellite capacity, 
% interpreted as the total average throughput that the satellite can deliver 
% when resources (power and bandwidth) are uniformly allocated across beams.

%% Satellite Parameters

% Maximum saturation power
Max_P_sat= 200*K/6;
% Max Saturation downlink power per Sat in Watts ( beam)
Max_Pb_Sat = 2*Max_P_sat/K;
% Reference value of Saturation downlink power per Sat in Watts
Ref_Pb_Sat = Max_P_sat/K; % 200 W for 6 beams
% Max Power per HPA
Max_Pamp=Ref_Pb_Sat*2*2;


% Maximum antenna Gain Satellite (dBi)
Max_G_Sat= 52;
% Satellite repeater output loss (dB)
Tx_sat_repeater_output_loss_dB = 2;
% Satellite repeater antenna loss (dB)
Tx_sat_repeater_antenna_loss_dB = 0.05;
% Downlink Polarization Loss (dB)
DL_polarization_loss_dB = 0.2;
% Maximum satellite antenna gain (dBi)
G_max_Sat= 52;

% Terminal Antenna Diameter (m)
terminal_antenna_diameter = 0.6;
% Terminal Antenna Efficiency (<1)
terminal_antenna_efficiency = 0.65;
% Losses due to terminal depointing (dB)
terminal_depointing_losses = 0.5;
% Terminal component losses (dB)
terminal_component_loss_dB = 0;
% Terminal component Noise Temperature (K)
terminal_component_T = 310;
% Noise Figure of the terminal LNB (dB)
terminal_LNA_NF_dB = 2;


symbol_rate= total_band;

% Average free space losse (dB)
FSL= 210;
% Average atmospheric losses (dB)
L_atm=0.4296;
% Ground Noise Temperature (K)
ground_Tground = 45;
% Average clear sky noise temperature at the receiver(K)
T_sky_av=28.4082;
% Average cloud noise temperature at the receiver(K)
T_cloud_av=0.6712;

%%% Constants
% Terminal temperature (K)
terminal_default_T = 290;
% Speed of light (m/s)
speed_of_light = 299792458;
% Boltzmann Constant (J/K)
boltzmann_constant = 1.3806503e-23;

%%% Input parameters used for the generation of the provided files, do not change
% Downlink Frequency Hz
DL_freq = 20e9;

%% Link Budget

% Compute the antenna temperature
T_ta_rx =T_sky_av+T_cloud_av;
% Compute the total noise temperature for each point on the grid
T_tot_rx    = T_ta_rx+ground_Tground+(10^(terminal_component_loss_dB/10)-1)*terminal_component_T+(10^(terminal_LNA_NF_dB/10)-1)*terminal_default_T/(10^(-terminal_component_loss_dB/10));

% Compute the maximum antenna gain for the user
terminal_antenna_Gmax = 10*log10(terminal_antenna_efficiency.*(pi*terminal_antenna_diameter*DL_freq/speed_of_light).^2);


% Compute the G/kT (in dB) for each grid point
GkT  = terminal_antenna_Gmax- FSL - DL_polarization_loss_dB - terminal_depointing_losses-10*log10(T_tot_rx)-10*log10(boltzmann_constant);
% Compute a refernce C/N factor (in dB) without considering the antenna pattern (assumic isotropic antenna) and without the transmit power
L_CN0_dB = -Tx_sat_repeater_output_loss_dB - Tx_sat_repeater_antenna_loss_dB + GkT - L_atm;

% Compute a refernce SNR factor taking into account the symbol rate
symbol_rate_dB = 10*log10(symbol_rate);
L_SNR_dB = L_CN0_dB - symbol_rate_dB;

%% Solver variables
Max_time_opt_second=10; % s, Maxium simulation time for the second-step process. One value for each simulated scenario
optionsGA.PoP= 2000*2; % Population for the genetic algorithm
optionsGA.Elite_num=ceil(0.1*optionsGA.PoP);% Number of elite members for the genetic algorithm
optionsGA.pmut=0.1; % Mutation probability for the genetic algorithm
optionsGA.Max_Gen=5000; % Maximum number of generations for the genetic algorithm


%% Batch simulations


for ind_batch=1:length(Nsims_v)
    Nsims=Nsims_v(ind_batch); % Obtain batch size
    for ind_scen=Nscen  % Select Scenario
    % for ind_scen=1:Nscen  % Select Scenario
        % Variables that stores the simulated data
        Nuser_beam_c=cell(1,Nsims);
        Req_user_c=cell(1,Nsims);
        R_off_POW_c=cell(1,Nsims);
        R_off_BW_c=cell(1,Nsims);
        R_off_CR_c=cell(1,Nsims);
        R_off_MAP_c=cell(1,Nsims);
        R_off_BW_CR_c=cell(1,Nsims);
        R_off_BW_CR_c_1=cell(1,Nsims);
        R_off_BW_MAP_c=cell(1,Nsims);
        R_off_BW_POW_c=cell(1,Nsims);
        Assig_MAP_c=cell(1,Nsims);
        Assig_CR_c=cell(1,Nsims);
        Assig_BW_CR_c=cell(1,Nsims);
        Assig_BW_CR_c_1=cell(1,Nsims);
        Assig_BW_MAP_c=cell(1,Nsims);
        users_locations_c=cell(1,Nsims);
        snr_car_c=cell(1,Nsims);
        user_beams_c=cell(1,Nsims);
        Bandwidth_Allo_BW=cell(1,Nsims);
        Bandwidth_Allo_BW_CR=cell(1,Nsims);
        Bandwidth_Allo_BW_MAP=cell(1,Nsims);
        Power_Allo_POW=cell(1,Nsims);
        Power_Allo_BW_POW=cell(1,Nsims);
        Bandwidth_Allo_BW_POW=cell(1,Nsims);
        
        res_IterNum_Av = 0;
        res_NQU_Av=0;
        res_NU_Av=0;
        res_OffRate_Av=0;
        res_MinUserRate_Av=0;
        
      sprintf('Start of Batch simulations %i of Scenario %i',ind_batch,ind_scen) % Display to keep track of the simulations     
      tic
        for indsims=1:Nsims % Monte-Carlo simulations for the selected scenario and batch size.
            
            sprintf('Monte-Carlo simulation %i of %i',indsims,Nsims) % Display to keep track of the simulations
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Number of user per beam    %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aux_alpha_v= alpha_v{ind_scen}; % Values to model the traffic. Dirichlet distribution
            
            % Obtain random number following a Dirichlet distribution
            n=1; % Number of generated scenarios
            % 从 gamma 分布中生成一个随机数数组
              r = gamrnd(repmat(aux_alpha_v,n,1),1,n,K); % Generate n vectors of longitude K following a gamma distribution Gamma(alpha_i,1)
            Dir_rand = r ./ repmat(sum(r,2),1,K); % Normalization to obtain random numbers following a Dirichlet distribution.
            
            % Obtain a integer number of users per beam
            Nuser_beam=round( Nuser_tot*Dir_rand);
            % 假设每个活跃用户要求25 Mbps的瞬时下行速率
            % 为了尝试用完整体容量T，在不同的波束上容纳Nuser_tot=272个同时活跃的用户
            if sum(Nuser_beam)>Nuser_tot
                in=1;
                while(in)
                    index=randi([1 K]);
                    if  Nuser_beam(index)>0
                        Nuser_beam(index)=Nuser_beam(index)-1;
                    end
                    if sum(Nuser_beam)==Nuser_tot
                        in=0;
                    end
                end
            elseif sum(Nuser_beam)<Nuser_tot
                in=1;
                while(in)
                    index=randi([1 K]);
                    Nuser_beam(index)=Nuser_beam(index)+1;
                    if sum(Nuser_beam)==Nuser_tot
                        in=0;
                    end
                end
            end
            N=sum(Nuser_beam);
            % Auxiliar variable 辅助变量
            cum_Nu_sim_index=cumsum(Nuser_beam);
            
            %%  User location generation   %
            % User are randomly placed within the beam radius

            user_beams=zeros(1,N);   % Vector that indicates the domminant beam for each user
            users_b_index=cell(K,1); % Cell with the user indexes for each beam
            users_b_loc=cell(K,1);
            for i=1:K
                x0_beam=Center_Beams(i,1);
                y0_beam=Center_Beams(i,2);
                
                if Nuser_beam(i)~=0
                    % Generate random user locations within a beam
                    t = 2*pi*rand(Nuser_beam(i),1);
                    r = R*sqrt(rand(Nuser_beam(i),1)); % Beam with radius R from Bessel modeling
                    % Obtain user position
                    x = x0_beam + r.*cos(t);
                    y = y0_beam + r.*sin(t);
                    % Generate auxiliar variables
                    switch i
                        case 1
                            users_b_index{i}=1:cum_Nu_sim_index(i);
                        otherwise
                            users_b_index{i}=cum_Nu_sim_index(i-1)+1:cum_Nu_sim_index(i);
                    end
                    user_beams(users_b_index{i})=i;
                    users_b_loc{i}=[x y];
                end
            end
            users_locations=cell2mat(users_b_loc);

            %%  Obtain channel values   %

            %%% Compute distances between user locations and beam centers for Bessel modeling
            distance_All=zeros(N,K);
            for i=1:K
                distance_All(:,i)= vecnorm(users_locations-Center_Beams(i,:),2,2);
            end
            
            %%% Aplay Bessel modeling
            u=2.07123*distance_All/R;
            indZ=find(~distance_All);

            g= besselj(1,u)./(2*u)+ 36*besselj(3,u)./(u.^3);
            G_dB=10*log10( abs(besselj(1,u)./(2*u)+ 36*besselj(3,u)./(u.^3)).^2  );
            G_dB(indZ)=zeros(1,length(indZ));
            
            %%% Magnitude Channel matrix, normalized respect to the Noise with bandwidth W
            Gamma_dB_wo_P=L_SNR_dB+(G_dB+G_max_Sat); %  in dB
            Gamma_wo_P=10.^(Gamma_dB_wo_P/10); % in natural
            
            % Carrier SNR for uniform resource allocation
            snr_car= Gamma_wo_P*Ref_Pb_Sat/0.5; % Carrier SNR
            snr_car_unfilter=snr_car; 
            SNR_car=10*log10(snr_car); % 272*6
            Th_car= 8.7000; % SNR threshold to obtain a C/I=24 dB or higher.
            ind_snr=find(SNR_car<Th_car);
            snr_car(ind_snr)=0;   % Filter SNR values
            
            %%  Auxiliar varibles for the resource managment  %

            % Requested traffic per user
            Req_user=Req_user_ref*ones(N,1);
            Req_user_c{ indsims}=Req_user+0;
            
            Total_Req=sum(Req_user); % Total requested traffic
            
            % Traffic requested per beam
            Req_b=zeros(1,K);

            for i=1:K
                users_index=users_b_index{i};
                if ~isempty(users_index)
                    Req_b(i)=sum(Req_user(users_index));
                end

            end

            % hold off;
            %% Structure with the simulated scenario
            scen_data.K=K;
            scen_data.M=M;
            scen_data.N=N;
            scen_data.Delta_W=Delta_W;
            scen_data.Nuser_beam=Nuser_beam;
            scen_data.user_beams=user_beams;
            scen_data.users_b_index=users_b_index;
            scen_data.Max_P_sat=Max_P_sat;
            scen_data.Ref_Pb_Sat=Ref_Pb_Sat;
            scen_data.Max_Pamp=Max_Pamp;
            scen_data.aux_Cn=(0.5/Ref_Pb_Sat)*snr_car_unfilter;
            scen_data.Req_beam=Req_b;
            scen_data.Req_user=Req_user;
            scen_data.Gamma_wo_P=Gamma_wo_P;
            scen_data.snr_car= snr_car;
            scen_data.Max_time_opt_second=Max_time_opt_second;
            scen_data.SE_av=SE_av;

            % scen_data.R_S=R_S;
            scen_data.R_0=R_0;
            scscen_data.Center_Beams=Center_Beams;
            scen_data.R_f=R_f;
            scen_data.distance_All=distance_All;
            scen_data.max_user_per_beam=max_user_per_beam;
            scen_data.users_locations=users_locations;
            scen_data.Center_Beams=Center_Beams;
            scen_data.L_SNR_dB=L_SNR_dB;
            
            % Create a "ResourceAssignmentOneDimesion" object  with the simulated data
            Sim_object= ResourceAssignmentTwoDimesion_r7_paper1(scen_data);
            
            %%              Resource assignment                %
   
            %% Flexible Bandwidth. Two-step optimization process.
            if en.BW
                sprintf('Solution: Flexible Bandwidth')
                [ R_off_BW,M_BW]=FlexibleBandwidth(Sim_object);
                % Store user rates
                R_off_BW_c{ indsims}=R_off_BW+0;
                % Store bandwidth allocation
                Bandwidth_Allo_BW{indsims}=M_BW';
                
                
                % Obtain the Normalized Quadratic unment demand (NQU)
                NQU_BW=sum(  ( Req_user- R_off_BW).^2  )*(total_band^2)/( N*(total_band*Req_user_ref)^2 );
                % Obtain the Normalized Unment demand (NU)
                NU_BW=sum(  ( Req_user- R_off_BW)  )/(SE_beam_uniform*K);
                % Obtain the total offered rate, Mbps
                Offer_rate_BW=sum( R_off_BW)*total_band/1e6;
                % Obtain the minimum user rate offered rate, Mbps
                Min_User_rate_BW=min( R_off_BW)*total_band/1e6;
                
            else
                NQU_BW=NaN;
                NU_BW=NaN;
                Offer_rate_BW=NaN;
                Min_User_rate_BW=NaN;
            end
            
            %% Flexible Power. Two-step optimization process
            if en.POW
                 sprintf('Solution: Flexible Power')
                [R_off_POW,P_POW]=FlexiblePower(Sim_object);
                % Store user rates
                R_off_POW_c{ indsims}=R_off_POW+0;
                % Store power allocation
                Power_Allo_POW{indsims}=P_POW;
                
                % Obtain the Normalized Quadratic unment demand (NQU)
                NQU_POW=sum(  ( Req_user- R_off_POW).^2  )*(total_band^2)/( N*(total_band*Req_user_ref)^2 );
                % Obtain the Normalized Unment demand (NU)
                NU_POW=sum(  ( Req_user- R_off_POW)  )/(SE_beam_uniform*K);
                % Obtain the total offered rate, Mbps
                Offer_rate_POW=sum( R_off_POW)*total_band/1e6;
                % Obtain the minimum user rate offered rate, Mbps
                Min_User_rate_POW=min( R_off_POW)*total_band/1e6;
            else
                NQU_POW=NaN;
                NU_POW=NaN;
                Offer_rate_POW=NaN;
                Min_User_rate_POW=NaN;
            end

            %% SR
            if en.SR
                sprintf('Solution: Flexible RangeCover')
                [R_off_CR,Assig_CR,iterationNum_CR] = FlexibleCoverRange(Sim_object);
                % Store user rates
                R_off_CR_c{indsims}=R_off_CR+0;
                % Store beam-user mapping
                Assig_CR_c{indsims}=Assig_CR';

                % Obtain the Normalized Quadratic unment demand (NQU)
                NQU_CR=sum(  ( Req_user- R_off_CR).^2  )*(total_band^2)/( N*(total_band*Req_user_ref)^2 );
                % Obtain the Normalized Unment demand (NU)
                NU_CR=sum(  ( Req_user- R_off_CR)  )/(SE_beam_uniform*K);
                % Obtain the total offered rate, Mbps
                Offer_rate_CR=sum( R_off_CR)*total_band/1e6;
                % Obtain the minimum user rate offered rate, Mbps
                Min_User_rate_CR=min( R_off_CR)*total_band/1e6;
            else
                NQU_CR=NaN;
                NU_CR=NaN;
                Offer_rate_CR=NaN;
                Min_User_rate_CR=NaN;
                iterationNum_CR = NaN;
            end

            %% Flexible Beam-user mapping with fixed resources. Two-step optimization process.
            if en.MAP
                sprintf('Solution: Flexible beam-user mapping with fixed resources')
                [R_off_MAP,Assig_MAP]=FixResFlexibleMapping(Sim_object);
                % Store user rates
                R_off_MAP_c{ indsims}=R_off_MAP+0;
                % Store beam-user mapping
                Assig_MAP_c{indsims}=Assig_MAP;
                
                % Obtain the Normalized Quadratic unment demand (NQU)
                NQU_MAP=sum(  ( Req_user- R_off_MAP).^2  )*(total_band^2)/( N*(total_band*Req_user_ref)^2 );
                % Obtain the Normalized Unment demand (NU)
                NU_MAP=sum(  ( Req_user- R_off_MAP)  )/(SE_beam_uniform*K);
                % Obtain the total offered rate, Mbps
                Offer_rate_MAP=sum( R_off_MAP)*total_band/1e6;
                % Obtain the minimum user rate offered rate, Mbps
                Min_User_rate_MAP=min( R_off_MAP)*total_band/1e6;
            else
                NQU_MAP=NaN;
                NU_MAP=NaN;
                Offer_rate_MAP=NaN;
                Min_User_rate_MAP=NaN;
            end
            
            %% SR_BW_alternate
            if en.SR_BW_alternate
                sprintf('Solution: Flexible bandwidth and RangeCover')
                [R_off_BW_CR,Assig_BW_CR,M_BW_CR,iterationNum_CRBW] = FlexBandwidthFlexibleCoverRange(Sim_object);
                % Store user rates
                R_off_BW_CR_c{indsims}=R_off_BW_CR+0;
                % Store beam-user mapping
                Assig_BW_CR_c{indsims}=Assig_BW_CR';
                % Store bandwidth allocation
                Bandwidth_Allo_BW_CR{indsims}=M_BW_CR';

                % Obtain the Normalized Quadratic unment demand (NQU)
                NQU_BW_CR=sum(  ( Req_user- R_off_BW_CR).^2  )*(total_band^2)/( N*(total_band*Req_user_ref)^2 );
                % Obtain the Normalized Unment demand (NU)
                NU_BW_CR=sum(  ( Req_user- R_off_BW_CR)  )/(SE_beam_uniform*K);
                % Obtain the total offered rate, Mbps
                Offer_rate_BW_CR=sum( R_off_BW_CR)*total_band/1e6;
                % Obtain the minimum user rate offered rate, Mbps
                Min_User_rate_BW_CR=min( R_off_BW_CR)*total_band/1e6;
            else
                NQU_BW_CR=NaN;
                NU_BW_CR=NaN;
                Offer_rate_BW_CR=NaN;
                Min_User_rate_BW_CR=NaN;
                iterationNum_CRBW = NaN;
            end
        
            %% Flexible Bandwidth and Beam-user mapping. Two-step optimization process.
            if en.BW_MAP
                sprintf('Solution: Flexible bandwidth and beam-user mapping') 
                [R_off_BW_MAP,Assig_BW_MAP,M_BW_MAP]=FlexBandwidthFlexMapping(Sim_object);
                % Store user rates
                R_off_BW_MAP_c{ indsims}=R_off_BW_MAP+0;
                % Store beam-user mapping
                Assig_BW_MAP_c{indsims}=Assig_BW_MAP;
                % Store bandwidth allocation
                Bandwidth_Allo_BW_MAP{indsims}=M_BW_MAP';
                
                % Obtain the Normalized Quadratic unment demand (NQU)
                NQU_BW_MAP=sum(  ( Req_user- R_off_BW_MAP).^2  )*(total_band^2)/( N*(total_band*Req_user_ref)^2 );
                % Obtain the Normalized Unment demand (NU)
                NU_BW_MAP=sum(  ( Req_user- R_off_BW_MAP)  )/(SE_beam_uniform*K);
                % Obtain the total offered rate, Mbps
                Offer_rate_BW_MAP=sum( R_off_BW_MAP)*total_band/1e6;
                % Obtain the minimum user rate offered rate, Mbps
                Min_User_rate_BW_MAP=min( R_off_BW_MAP)*total_band/1e6;
            else
                NQU_BW_MAP=NaN;
                NU_BW_MAP=NaN;
                Offer_rate_BW_MAP=NaN;
                Min_User_rate_BW_MAP=NaN;
            end
            
            % Flexible Bandwidth and Power. Genetic Algorithm
            if en.BW_POW
                sprintf('Solution: Flexible bandwidth and power') 
                
                % Default values if other strategies with fixed mapping are not simualted.
                %可通过控制flag选择要仿真的策略
                if ~en.POW
                    P_POW=Ref_Pb_Sat*ones(1,K);
                end
                
                if ~en.BW
                    M_BW=M*ones(K,1);
                end
                
                [R_off_BW_POW,M_BW_POW,P_BW_POW]=FlexBandwidthPower(Sim_object,P_POW,M_BW,optionsGA);
                % Store user rates
                R_off_BW_POW_c{indsims} =R_off_BW_POW+0;
                % Store bandwidth allocation
                Bandwidth_Allo_BW_POW{indsims}=M_BW_POW;
                % Store power allocation
                Power_Allo_BW_POW{indsims}=P_BW_POW;
                
                % Obtain the Normalized Quadratic unment demand (NQU)
                NQU_BW_POW=sum(  ( Req_user- R_off_BW_POW).^2  )*(total_band^2)/( N*(total_band*Req_user_ref)^2 );
                % Obtain the Normalized Unment demand (NU)
                NU_BW_POW=sum(  ( Req_user- R_off_BW_POW)  )/(SE_beam_uniform*K);
                % Obtain the total offered rate, Mbps
                Offer_rate_BW_POW=sum( R_off_BW_POW)*total_band/1e6;
                % Obtain the minimum user rate offered rate, Mbps
                Min_User_rate_BW_POW=min( R_off_BW_POW)*total_band/1e6;
   
                
            else
                NQU_BW_POW=NaN;
                NU_BW_POW=NaN;
                Offer_rate_BW_POW=NaN;
                Min_User_rate_BW_POW=NaN;
            end
            
            
            % Store basic data to replicate the scenario
            user_beams_c{ indsims}=user_beams;
            Nuser_beam_c{ indsims}=Nuser_beam;
            users_locations_c{ indsims}=users_locations;
            
            
            % Display Normalized Quadratic unment demand (NQU)
          
            sprintf('Simulation Results:')
            res_IterNum = [  iterationNum_CR iterationNum_CRBW];
            res_NQU = [ NQU_POW NQU_BW NQU_BW_POW NQU_MAP NQU_BW_MAP NQU_CR NQU_BW_CR]';
            res_NU = [ NU_POW NU_BW NU_BW_POW NU_MAP NU_BW_MAP NU_CR NU_BW_CR]';
            res_OffRate = [ Offer_rate_POW Offer_rate_BW Offer_rate_BW_POW Offer_rate_MAP Offer_rate_BW_MAP Offer_rate_CR  Offer_rate_BW_CR]';
            res_MinUserRate = [  Min_User_rate_POW Min_User_rate_BW Min_User_rate_BW_POW Min_User_rate_MAP Min_User_rate_BW_MAP Min_User_rate_CR Min_User_rate_BW_CR]';
            
            T = table(res_NQU,res_NU,res_OffRate,res_MinUserRate,'VariableNames', ...
                {'NQU','NU','Offered Rate [Mps]','Min. User Rate [Mbps]'}, ...
                'RowName',{'POW','BW','BW-POW','MAP','BW-MAP','CR','CR-BW'});
            disp(T)
            
            res_IterNum_Av = res_IterNum_Av + res_IterNum;
            res_NQU_Av=res_NQU_Av+res_NQU;
            res_NU_Av=res_NU_Av+res_NU;
            res_OffRate_Av=res_OffRate_Av+res_OffRate;
            res_MinUserRate_Av=res_MinUserRate_Av+res_MinUserRate;
        toc
        end % loop Nsims
      
         sprintf('End of Batch simulations %i of Scenario %i',ind_batch,ind_scen) % Display to keep track of the simulations
      
         
         sprintf('Batch simulations Results:')
            T = table(res_NQU_Av/Nsims,res_NU_Av/Nsims,res_OffRate_Av/Nsims,res_MinUserRate_Av/Nsims, ...
                'VariableNames',{'NQU','NU','Offered Rate [Mps]','Min. User Rate [Mbps]'}, ...
                'RowName',{'POW','BW','BW-POW','MAP','BW-MAP','CR','CR-BW'});
            disp(T)
            IterNum_Av = res_IterNum_Av/Nsims;
            disp(IterNum_Av)
            
            
        sprintf('Storing Batch simulation results')
        save([ parentpath  '\' label_alpha{ind_scen} '\'  text_file  '_Batch_' num2str(ind_batch) '_alpha_' label_alpha{ind_scen} '.mat'],...
            'Nsims_v','M','K','Delta_W','Nuser_beam_c','Req_user_c','R_off_POW_c','R_off_BW_c','R_off_CR_c','R_off_MAP_c',...
            'R_off_BW_CR_c','R_off_BW_CR_c_1','R_off_BW_MAP_c','Bandwidth_Allo_BW','Bandwidth_Allo_BW_CR','Bandwidth_Allo_BW_MAP','Power_Allo_POW','Power_Allo_BW_POW',...
            'Bandwidth_Allo_BW_POW','user_beams_c','users_locations_c','Assig_BW_CR_c','Assig_BW_CR_c_1','Assig_MAP_c','Assig_CR_c','Assig_BW_MAP_c','R_off_BW_POW_c','IterNum_Av')
 
        
    end % loop Nscen
    
    
end % loop batch


warning(orig_state)
