% The script plot the average performance in the console window and displays different cumulativedistribution functions of the performance metrics.
% If multiple PC were employed to perform the simulations, all the data have to be stored in a folder with the same name as the label given to the scenario.
%
% This Matlab script was developed to generate simulation results to:
%
% Tomas Ramirez, Carlos Mosquera, Nader Alagha,
% "Flexible User Mapping for Radio Resource Assignment in Advanced Satellite Payloads"
% arXiv, 2021
%
% License: This code is licensed under the GPLv3 license. If you in any way
% use this code for research that results in publications, please cite our
% work as described above.
%%
clear all
close all
clc
%% Scenario selection

disp('Please, select a scenario subfolder with the stored numerical results');
selpath = uigetdir;

Scenario_label=selpath(find('\'==selpath',1,'last')+1:end);
 
% Replace "_" for a white space 
Scenario_label_text= Scenario_label;
Scenario_label_text( Scenario_label_text=='_')=' ';


%% Auxiliar variables to store the simulated data

Nuser_beam_c_aux=[];
Req_user_c_aux=[];
user_beams_c_aux=[];

R_off_POW_c_aux=[];
R_off_BW_c_aux=[];
R_off_CR_c_aux=[];
R_off_MAP_c_aux=[];
R_off_BW_CR_c_aux=[];
R_off_BW_MAP_c_aux=[];
R_off_BW_POW_c_aux=[];

Assig_BW_MAP_c_aux=[];

M_pre=[];
K_pre=[];
%% Process simulated results

filelist=dir([ selpath '\*' Scenario_label '.mat']);

if length(filelist)==0
    
    error('The selected folder does not contain any data')
end

for ind_file=1:length(filelist)
    
    load([ filelist(ind_file).folder '\' filelist(ind_file).name])
    
    if isempty(M_pre) % If it is empty, this is the first processed file
        M_pre=M;
        K_pre=K;
        
        en_POW=~any(cellfun( @(x) isempty(x),R_off_POW_c));
        en_BW=~any(cellfun( @(x) isempty(x),R_off_BW_c));
        en_CR=~any(cellfun( @(x) isempty(x),R_off_CR_c));
        en_MAP=~any(cellfun( @(x) isempty(x),R_off_MAP_c));
        en_BW_CR=~any(cellfun( @(x) isempty(x),R_off_BW_CR_c));
        en_BW_MAP=~any(cellfun( @(x) isempty(x),R_off_BW_MAP_c));
        en_BW_POW=~any(cellfun( @(x) isempty(x),R_off_BW_POW_c));
    else
        
        if (M_pre~=M)|| (K_pre~=K)
            error('Different simulation parameters');
        end
        
    end
    
    
    
    Nuser_beam_c_aux=[Nuser_beam_c_aux Nuser_beam_c];
    Req_user_c_aux=[Req_user_c_aux  Req_user_c];
    user_beams_c_aux=[user_beams_c_aux user_beams_c];
    
    if en_POW
        R_off_POW_c_aux=[R_off_POW_c_aux R_off_POW_c];
    end
    
    if en_BW
        R_off_BW_c_aux=[R_off_BW_c_aux R_off_BW_c];
    end
    
    if en_CR
        R_off_CR_c_aux=[R_off_CR_c_aux R_off_CR_c];
    end

    if en_MAP
        R_off_MAP_c_aux=[R_off_MAP_c_aux R_off_MAP_c];
    end
    
    if en_BW_CR
        R_off_BW_CR_c_aux=[R_off_BW_CR_c_aux R_off_BW_CR_c];
    end

    if en_BW_MAP
        R_off_BW_MAP_c_aux=[R_off_BW_MAP_c_aux R_off_BW_MAP_c];
    end
    
    if en_BW_POW
        R_off_BW_POW_c_aux=[R_off_BW_POW_c_aux R_off_BW_POW_c];
    end
    
    
    
end

Nuser_beam_c=Nuser_beam_c_aux;
Req_user_c=Req_user_c_aux ;
user_beams_c=user_beams_c_aux ;

R_off_POW_c=R_off_POW_c_aux ;
R_off_BW_c=R_off_BW_c_aux ;

R_off_CR_c=R_off_CR_c_aux ;
R_off_MAP_c=R_off_MAP_c_aux ;
R_off_BW_CR_c=R_off_BW_CR_c_aux ;
R_off_BW_MAP_c=R_off_BW_MAP_c_aux ;

R_off_BW_POW_c=R_off_BW_POW_c_aux;


clear Nuser_beam_c_aux  Req_user_c_aux R_off_POW_c_aux ...
    R_off_BW_c_aux R_off_CR_c_aux R_off_MAP_c_aux R_off_BW_CR_c_aux R_off_BW_MAP_c_aux R_off_BW_POW_c_aux ...
    user_beams_c_aux 
 

%% Auxiliar variables to plot the results
total_band=500e6;
SE_av=4.5271; % log2(1+ SNR_eff_beam) , with SNR_eff_beam  the average SNR for a uniform distribution of user within a beam (  approx 13.5 dB )

tech_label={'POW', 'BW', 'SR','MAP','BW-SR','BW-MAP','BW-POW'};
marker_v=['v','o','s','^','+','d','*','x','x','v','p','.'];
Ntech=length(tech_label);
c6 = categorical(tech_label);
c6 = reordercats(c6,tech_label);


sim_sol= find([en_POW en_BW en_CR en_MAP en_BW_CR en_BW_MAP en_BW_POW]);


NQU_sims=cell(1,Ntech);
NQU_av=ones(1,5)*NaN;

NU_sims=cell(1,Ntech);
NU_av=ones(1,5)*NaN;

Offered_sims=cell(1,Ntech);
Offered_av=ones(1,5)*NaN;

Min_user_sims=cell(1,Ntech);
Min_user_av=ones(1,5)*NaN;


colors(1,:)=[111 111 111]/256;
colors(2,:)=[91 181 172]/256;
colors(3,:)=[84 123 180]/256;
colors(4,:)=[98 156 53]/256;
colors(5,:)=[221 124 79]/256;
colors(6,:)=[192 50 26]/256;
colors(7,:)=[108 97 175]/256;


% Variables for pdf/cdf computation
F=cell(1,Ntech);
f=cell(1,Ntech);
x=cell(1,Ntech);
factor_amp=25;
stat_band=0.15;
Npoints=5e3;

N_total=272;

R_av=25/500; % Traffic requested per user, normalized to the avalaible bandwidth : 25 Mbps/ 500 MHz

% Performance functions
NQU_fun=@(x,y) ((x-y).^2)*(total_band)^2/( N_total*(total_band*R_av)^2 );
NU_fun=@(x,y) ((x-y))/((K/2)*SE_av);
Offered_fun=@(x) sum(x)*total_band/1e6;   % Mbps
Min_user_fun=@(x) min(x)*total_band/1e6; %Mbps


%% Beam Traffic demand std

std_Req_beam=cellfun( @(x) std( x*R_av*500), Nuser_beam_c);

[ff,xx] = ksdensity( std_Req_beam);
FF=cumsum(ff)/sum(ff);

figure,plot(xx,ff,'LineWidth',2)
xlabel(' Beam Traffic Demand std, Mbps')
ylabel('pdf')
grid on
grid minor
xlim([0 1e3*sqrt(5)*6.8/6])
title([ Scenario_label_text ': Standard deviation of the traffic distribution'])
%% Plot normalized quadtratic unmet (NQU) cdf

% Obtain techniques performance
if en_POW
    NQU_sims{1}=cellfun(@(x,y) sum(NQU_fun(x,y)),Req_user_c,R_off_POW_c);
    NQU_av(1)=mean(NQU_sims{1});
end

if en_BW
    NQU_sims{2}=cellfun(@(x,y) sum(NQU_fun(x,y)),Req_user_c,R_off_BW_c);
    NQU_av(2)=mean(NQU_sims{2});
end

if en_CR
    NQU_sims{3}=cellfun(@(x,y) sum(NQU_fun(x,y)),Req_user_c,R_off_CR_c);
    NQU_av(3)=mean(NQU_sims{3});
    
end

if en_MAP
    NQU_sims{4}=cellfun(@(x,y) sum(NQU_fun(x,y)),Req_user_c,R_off_MAP_c);
    NQU_av(4)=mean(NQU_sims{4});
    
end

if en_BW_CR
    NQU_sims{5}=cellfun(@(x,y) sum(NQU_fun(x,y)),Req_user_c,R_off_BW_CR_c);
    NQU_av(5)=mean(NQU_sims{5});
end

if en_BW_MAP
    NQU_sims{6}=cellfun(@(x,y) sum(NQU_fun(x,y)),Req_user_c,R_off_BW_MAP_c);
    NQU_av(6)=mean(NQU_sims{6});
end

if en_BW_POW
    NQU_sims{7}=cellfun(@(x,y) sum(NQU_fun(x,y)),Req_user_c,R_off_BW_POW_c);
    NQU_av(7)=mean(NQU_sims{7});
end


% Plot cdf
figure,hold on
for ind_tech=sim_sol
    
    [f,x] = ksdensity(factor_amp*NQU_sims{ind_tech},'NumPoints',Npoints,'Bandwidth',stat_band);
    F =cumsum(f)/sum(f);
    
    plot(x/factor_amp,F ,'LineWidth',2,'Color',colors(ind_tech,:),'DisplayName',tech_label{ind_tech},...
        'Marker',marker_v(ind_tech),'MarkerIndices', 1+20:200:length(x) );
    
    legend('Location','SouthEast');
    
    grid on
    grid minor
    ylabel('cdf')
    xlabel('NQU')
end

title([ Scenario_label_text ': NQU cumulative distribution function'])
%% Plot normalzied Unmet (NU) cdf

% Obtain techniques performance
if en_POW
    NU_sims{1}=cellfun(@(x,y) sum(NU_fun(x,y)),Req_user_c,R_off_POW_c);
    NU_av(1)=mean(NU_sims{1});
end

if en_BW
    NU_sims{2}=cellfun(@(x,y) sum(NU_fun(x,y)),Req_user_c,R_off_BW_c);
    NU_av(2)=mean(NU_sims{2});
end

if en_CR
    NU_sims{3}=cellfun(@(x,y) sum(NU_fun(x,y)),Req_user_c,R_off_CR_c);
    NU_av(3)=mean(NU_sims{3});
    
end

if en_MAP
    NU_sims{4}=cellfun(@(x,y) sum(NU_fun(x,y)),Req_user_c,R_off_MAP_c);
    NU_av(4)=mean(NU_sims{4});
    
end

if en_BW_CR
    NU_sims{5}=cellfun(@(x,y) sum(NU_fun(x,y)),Req_user_c,R_off_BW_CR_c);
    NU_av(5)=mean(NU_sims{5});
end

if en_BW_MAP
    NU_sims{6}=cellfun(@(x,y) sum(NU_fun(x,y)),Req_user_c,R_off_BW_MAP_c);
    NU_av(6)=mean(NU_sims{6});
end

if en_BW_POW
    NU_sims{7}=cellfun(@(x,y) sum(NU_fun(x,y)),Req_user_c,R_off_BW_POW_c);
    NU_av(7)=mean(NU_sims{7});
end


pts=0:1e-3:factor_amp;


% Plot cdf
figure,hold on
for ind_tech=sim_sol
    
    [f,x] = ksdensity(factor_amp*NU_sims{ind_tech},pts);
    F =cumsum(f)/sum(f);
    
    plot(x/factor_amp,F ,'LineWidth',2,'Color',colors(ind_tech,:),'DisplayName',tech_label{ind_tech},...
        'Marker',marker_v(ind_tech),'MarkerIndices', 1+20:600:length(x) );
    
    legend('Location','SouthEast');
    
    grid on
    grid minor
    ylabel('cdf')
    xlabel('NU')
end


title([ Scenario_label_text ': NU cumulative distribution function'])



%% Plot Total Offered rate cdf

% Obtain techniques performance
if en_POW
    Offered_sims{1}=cellfun(@(x) Offered_fun(x),R_off_POW_c);
    Offered_av(1)=mean(Offered_sims{1});
end

if en_BW
    Offered_sims{2}=cellfun(@(x) Offered_fun(x),R_off_BW_c);
    Offered_av(2)=mean(Offered_sims{2});
end

if en_CR
    Offered_sims{3}=cellfun(@(x) Offered_fun(x),R_off_CR_c);
    Offered_av(3)=mean(Offered_sims{3});
    
end

if en_MAP
    Offered_sims{4}=cellfun(@(x) Offered_fun(x),R_off_MAP_c);
    Offered_av(4)=mean(Offered_sims{4});
    
end

if en_BW_CR
    Offered_sims{5}=cellfun(@(x) Offered_fun(x),R_off_BW_CR_c);
    Offered_av(5)=mean(Offered_sims{5});
end

if en_BW_MAP
    Offered_sims{6}=cellfun(@(x) Offered_fun(x),R_off_BW_MAP_c);
    Offered_av(6)=mean(Offered_sims{6});
end

if en_BW_POW
    Offered_sims{7}=cellfun(@(x) Offered_fun(x),R_off_BW_POW_c);
    Offered_av(7)=mean(Offered_sims{7});
end


% Plot cdf
pts=2000:1:6800;
figure,hold on
plot(total_band/1e6*[sum(Req_user_c{1}) sum(Req_user_c{1})],[0 1],'k','LineWidth',2,'DisplayName','T, satellite capacity')
for ind_tech=sim_sol
    
    [f,x] = ksdensity(Offered_sims{ind_tech},pts);
    F =cumsum(f)/sum(f);
    
    plot(x,F ,'LineWidth',2,'Color',colors(ind_tech,:),'DisplayName',tech_label{ind_tech},...
        'Marker',marker_v(ind_tech),'MarkerIndices', 1+20:200:length(x) );
    
    legend('Location','northwest');
    
    grid on
    grid minor
    ylabel('cdf')
    xlabel('Total Offered Rate [Mbps]')
end


title([ Scenario_label_text ': Total offered rate cumulative distribution function'])

%% Min

% Obtain  performances
if en_POW
    Min_user_sims{1}=cellfun(@(x) Min_user_fun(x),R_off_POW_c);
    Min_user_av(1)=mean(Min_user_sims{1});
end

if en_BW
    Min_user_sims{2}=cellfun(@(x)Min_user_fun(x),R_off_BW_c);
    Min_user_av(2)=mean(Min_user_sims{2});
end

if en_CR
    Min_user_sims{3}=cellfun(@(x)Min_user_fun(x),R_off_CR_c);
    Min_user_av(3)=mean(Min_user_sims{3});
    
end

if en_MAP
    Min_user_sims{4}=cellfun(@(x)Min_user_fun(x),R_off_MAP_c);
    Min_user_av(4)=mean(Min_user_sims{4});
    
end

if en_BW_CR
    Min_user_sims{5}=cellfun(@(x)Min_user_fun(x),R_off_BW_CR_c);
    Min_user_av(5)=mean(Min_user_sims{5});
end

if en_BW_MAP
    Min_user_sims{6}=cellfun(@(x)Min_user_fun(x),R_off_BW_MAP_c);
    Min_user_av(6)=mean(Min_user_sims{6});
end

if en_BW_POW
    Min_user_sims{7}=cellfun(@(x)Min_user_fun(x),R_off_BW_POW_c);
    Min_user_av(7)=mean(Min_user_sims{7});
end


% Plot cdf
bw_kernel=0.01;
pts=0:bw_kernel/10:R_av*500;
figure,hold on
for ind_tech=sim_sol
    
    [f,x] = ksdensity(Min_user_sims{ind_tech},pts,'Bandwidth',1);
    F =cumsum(f)/sum(f);
    
    plot(x,F ,'LineWidth',2,'Color',colors(ind_tech,:),'DisplayName',tech_label{ind_tech},...
        'Marker',marker_v(ind_tech),'MarkerIndices', 1+20:900:length(x) );
    
    legend('Location','northwest');
    
    grid on
    grid minor
    ylabel('cdf')
    xlabel('Minimum user rate [Mbps]')
end


title([ Scenario_label_text ': Minimum user rate cumulative distribution function'])



%% User offered rate cdf


if en_POW
    aux_Rates{1}=total_band/1e6*cell2mat(R_off_POW_c(:));
end

if en_BW
    aux_Rates{2}=total_band/1e6*cell2mat(R_off_BW_c(:));
end

if en_CR
    aux_Rates{3}=total_band/1e6*cell2mat(R_off_CR_c(:));
end

if en_MAP
    aux_Rates{4}=total_band/1e6*cell2mat(R_off_MAP_c(:));
end

if en_BW_CR
    aux_Rates{5}=total_band/1e6*cell2mat(R_off_BW_CR_c(:));
end

if en_BW_MAP
    aux_Rates{6}=total_band/1e6*cell2mat(R_off_BW_MAP_c(:));
end

if en_BW_POW
    aux_Rates{7}=total_band/1e6*cell2mat(R_off_BW_POW_c(:));
end



% Plot cdf
bw_kernel=0.01;
pts=0:bw_kernel/10:R_av*500+5;
figure,hold on
plot(total_band/1e6*[R_av R_av],[0 1],'k','LineWidth',2,'DisplayName','Requested')
for ind_tech=sim_sol
    
    [f,x] = ksdensity(aux_Rates{ind_tech},pts,'Bandwidth',bw_kernel,'BoundaryCorrection','reflection');
    F=cumsum(f)/sum(f);
    
    plot(x,F,'LineWidth',2,'Color',colors(ind_tech,:),'DisplayName',tech_label{ind_tech},...
        'Marker',marker_v(ind_tech),'MarkerIndices', 1+20:1000:length(x) );
    
    legend('Location','SouthEast');
    
    
end


ylabel('cdf')
xlabel('Offered user rate [Mbps]')
legend('Location','NorthWest')
grid on
grid minor
xlim([0 500*R_av+5])


title([ Scenario_label_text ': Offered user rate cumulative distribution function'])


%% Average requested and offerd traffic per beam


Nsims=size(Req_user_c,2);
Offer_beam=cell(1,Ntech);


Req_b=zeros(Nsims,K);
 
 
for ind_sim=1:Nsims
 
    user_dombeam=user_beams_c{ind_sim};
    
    for i=1:K
        
        user_index=find(user_dombeam==i);
        
        if en_POW
            Offer_beam{1}(ind_sim,i)= total_band/1e6*sum ( R_off_POW_c{ind_sim}(user_index));
        end
        
        if en_BW
           Offer_beam{2}(ind_sim,i)= total_band/1e6*sum ( R_off_BW_c{ind_sim}(user_index));
        end
        
        if en_MAP
            Offer_beam{3}(ind_sim,i)= total_band/1e6*sum ( R_off_CR_c{ind_sim}(user_index));
        end

        if en_MAP
            Offer_beam{4}(ind_sim,i)= total_band/1e6*sum ( R_off_MAP_c{ind_sim}(user_index));
        end

        if en_BW_CR
            Offer_beam{5}(ind_sim,i)= total_band/1e6*sum ( R_off_BW_CR_c{ind_sim}(user_index));
        end
        
        if en_BW_MAP
            Offer_beam{6}(ind_sim,i)= total_band/1e6*sum ( R_off_BW_MAP_c{ind_sim}(user_index));
        end
        
        if en_BW_POW
             Offer_beam{7}(ind_sim,i)= total_band/1e6*sum ( R_off_BW_POW_c{ind_sim}(user_index));
        end
 
    end
    Req_b(ind_sim,:)= Nuser_beam_c{ind_sim}*25;
end

res=[];
for ind_tech=sim_sol   
    res= [ res ; mean(Offer_beam{ind_tech},1)];   
end
tech_label2={'Offered: POW', 'Offered: BW','Offered: SR', 'Offered: MAP','Offered: BW-SR','Offered: BW-MAP','Offered: BW-POW'};
tech_label3={'Requested'};
tech_label3(2:1+length(sim_sol))=tech_label2(sim_sol);

% %% just for 64 beams
% figure,hold on
% t = tiledlayout(8,1);
% for i = 1:8
%     nexttile
%     Ybar = [  mean(Req_b,1) ;res ]';
%     hbar=bar((i-1)*8+1:i*8,Ybar((i-1)*8+1:i*8,:));
%     set(gca,'yscale','log')
%     hbar(1).FaceColor=[ 0 0 0];
%     for ind_tech=1:length(sim_sol)
%         hbar(1+ind_tech).FaceColor= colors(sim_sol(ind_tech),:);
%     end
%     grid on
%     grid minor
% end
% xlabel(t,'Beam Index')
% ylabel(t,'Average Rates [Mbps]')
% t.TileSpacing = 'tight';
% t.Padding = 'tight';
% title(t,[ Scenario_label_text ': Offered and Requested rates per beam'])
% legend(tech_label3);

%% for 6 beams
figure,hbar=bar(1:6,[  mean(Req_b,1) ;res ]');
ylabel('Average Rates [Mbps]')
xlabel('Beam Index')

legend(tech_label3);

grid on
grid minor

hbar(1).FaceColor=[ 0 0 0];
for ind_tech=1:length(sim_sol)   
    hbar(1+ind_tech).FaceColor= colors(sim_sol(ind_tech),:);
end

title([ Scenario_label_text ': Offered and Requested rates per beam'])



%% Display metrics

disp([ '  Average performance : ' Scenario_label_text])
T = table(NQU_av',NU_av',Offered_av',Min_user_av', ...
    'VariableNames',{'NQU','NU','Offered Rate [Mps]', ...
    'Min. User Rate [Mbps]'},'RowName',{'POW','BW','SR','MAP','BW-CR','BW-MAP','BW-POW'});
disp(T)





