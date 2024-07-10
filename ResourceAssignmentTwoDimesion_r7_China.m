classdef ResourceAssignmentTwoDimesion_r7_China

    properties
        K % Number of beams
        M % Number of carriers per color
        N % Number of users
        Delta_W; % Carrier Bandwidth ( normalized to the total bandwidth)

        Nuser_beam %  Vector with the number of user per beam 每个波束的用户数量
        user_beams % Vector that indicates the dominant beam for each user 每个用户的主波束
        users_b_index % Cell with the user indexes for each beam 单元，包含每个波束的用户索引
        Gamma_wo_P  % Matrix with the channel values 信道参数
        snr_car % Carrier SNR for uniform carrier power allocation 载波信噪比，用于均匀载波功率分配。
        Max_P_sat % Overall power constraint 全部功率约束
        Ref_Pb_Sat % Value for uniform power allocation 均匀功率分布值
        Max_Pamp % Maximum amplifier power 最大放大器功率
        L_SNR_dB



        Req_beam  % Vector with the user requested traffic 用户请求流量
        Req_user % Vector with the beam requested traffic 波束请求流量

        Max_time_opt_second
        % Auxiliar variable for the second step process.
        % Sets the maximum time for the optimization.
        effective_snr_beam
        % Auxiliar vector for the beam resource assignment.
        % Effective SNR per beam
        tol
        % Auxiliar variable. Tolerance value 容差值

        SE_av
        % Auxilixar variable.  log2(1+ SNR_eff_beam) , with SNR_eff_beam
        % the average SNR for a uniform distribution of user within a beam (  approx 13.5 dB )
        % 波束内用户均匀分布的平均信噪比
        User_per_car
        % Auxiliar variable.
        % Estimation of the number of user which can be served by one carrier
        % 单载波估计可以服务的用户数量
        ub_B
        % Auxiliar variable.
        % Upper bound for the bandwidth allocation in the genetic algorithm
        % 遗传算法中带宽分配的上界
        aux_Cn
        % Auxiliar variable for the genetic algorithm. Matrix with SNR values.
        % 遗传算法中信噪比矩阵

        R_f
        R_0
        distance_All
        max_user_per_beam
        users_locations
        Center_Beams

    end


    methods
        function obj = ResourceAssignmentTwoDimesion_r7_China(arg)
            % Construct function. Load the scenario data from the argument into the class object
            obj.K=arg.K;
            obj.M= arg.M;
            obj.N=arg.N;
            obj.Delta_W=arg.Delta_W;

            obj.Nuser_beam=  arg.Nuser_beam;
            obj.user_beams= arg.user_beams;
            obj.users_b_index= arg.users_b_index;
            obj.Gamma_wo_P=  arg.Gamma_wo_P;
            obj.snr_car=  arg.snr_car;
            obj.Max_P_sat=  arg.Max_P_sat;
            obj.Ref_Pb_Sat= arg.Ref_Pb_Sat;
            obj.Max_Pamp= arg.Max_Pamp;
            obj.L_SNR_dB=arg.L_SNR_dB;

            obj.Req_user=arg.Req_user;
            obj.Req_beam = arg.Req_beam;
            obj.Max_time_opt_second=arg.Max_time_opt_second;
            obj.tol=1e-5;

            obj.effective_snr_beam= ComputEffectiveSNR(obj);

            obj.SE_av=arg.SE_av;
            obj.User_per_car= round(1./(obj.Req_user(1)/(obj.SE_av/(2*obj.M))));
            obj.ub_B =   min(ceil(obj.Nuser_beam'/obj.User_per_car),2*obj.M) ;
            obj.aux_Cn=arg.aux_Cn;

            obj.R_f= arg.R_f;
            obj.R_0= arg.R_0;
            obj.distance_All= arg.distance_All;
            obj.max_user_per_beam= arg.max_user_per_beam;
            obj.users_locations = arg.users_locations;
            obj.Center_Beams=arg.Center_Beams;

        end

        %% CR
        function [R_off_CR,Nuser_beam_CR,iterationNum_CR] = FlexibleCoverRange(obj)
            % 定义初始负载均衡度 loadBalancing_0
            for index_b = 1:obj.K
                to_average_user_load(index_b)= abs(45-obj.Nuser_beam(index_b))/45;
            end
            loadBalancing_0 = sum(to_average_user_load);
            % 初始信噪比
            snr_car_CR = obj.snr_car;
            flag_loadBalancing = 0;

            % 初始波束覆盖图
            figure('name','uniform')
            scatter(obj.users_locations(:,1),obj.users_locations(:,2),...
                6,'o','MarkerEdgeColor','none','MarkerFaceColor',"#DD7C4F");         
            hold on
            for ind_beam  =1:obj.K
                pos = [obj.Center_Beams(ind_beam,1)-obj.R_0, ...
                    obj.Center_Beams(ind_beam,2)-obj.R_0, ...
                    2*obj.R_0, 2*obj.R_0];
                r1 = rectangle('Position',pos,'Curvature',[1 1]);
                r1.LineStyle = '-';
                r1.EdgeColor = '#6F6F6F';
                r1.LineWidth = 1.2;
                hold on
            end
           
            title('Initial uniform beam pattern')
            % ax = gca;
            % ax.XLim = [0 710];
            % xlim([0 720])
            % ylim([-50 800])
            % xticks(0:50:750)
            % yticks(-50:50:800)
            xlim tight
            xlabel('x,km')
            ylabel('y,km')        
            axis equal

            % 当负载均衡度没有比上一次负载均衡度小时，跳出循环
            iterationNum_CR = 0;
            tic
            while flag_loadBalancing == 0
                iterationNum_CR = iterationNum_CR+1;
                %% FIRST STEP: Beam resource allocation
                cvx_solver default
                cvx_begin quiet
                variable w_beam_user(obj.K,obj.N)
                expression Roff_sim2(obj.K,obj.N)
                expression obj_fun
                for index_beam=1:obj.K
                    for index_user=1:obj.N
                        Roff_sim2( index_beam,index_user)= ...
                            w_beam_user(index_beam,index_user)*obj.Delta_W* ...
                            log(1 + snr_car_CR(index_user,index_beam) )/log(2);
                    end
                end
                obj_fun= ( obj.Req_user.' - sum(Roff_sim2,1)).^2;
                minimize( sum(obj_fun) );
                subject to
                0<=w_beam_user<=1
                sum(w_beam_user,1)<= 1 %  Constraint of one carrier per user
                sum(w_beam_user,2)<=obj.M % Constraint of fixed bandwidth ( M carriers) per beam
                cvx_end

                % Number of user per beam
                Assig_UM_FlexB=zeros(1,obj.N);
                for ind_n=1:obj.N
                    [~,indmax]=max( Roff_sim2(:,ind_n) );
                    Assig_UM_FlexB(ind_n)= indmax;
                end

                % Obtain association between users and beams.柔性映射
                Nuser_beam_CR=hist(Assig_UM_FlexB,1:obj.K);

                % 改变波束中心，画出波束覆盖图
                figure('name','CR')
                Center_Beams_CR = zeros(obj.K,2);
                
                for ind_beam=1:obj.K
                    if Nuser_beam_CR(ind_beam)> 0
                        users_index=find(Assig_UM_FlexB==ind_beam);
                        x_CR = obj.users_locations(users_index,1);
                        y_CR = obj.users_locations(users_index,2);
                        scatter(x_CR,y_CR,...
                            6,'o','MarkerEdgeColor','none','MarkerFaceColor',"#DD7C4F");
                        hold on
                        [center_CR, radius_CR] = min_cover_circle(x_CR, y_CR, Nuser_beam_CR(ind_beam));
                       
                        pos = [center_CR(1)-radius_CR, ...
                            center_CR(2)-radius_CR, ...
                            2*radius_CR, 2*radius_CR];
                        r1 = rectangle('Position',pos,'Curvature',[1 1]);
                        r1.LineStyle = '-';
                        r1.EdgeColor = '#6F6F6F';
                        r1.LineWidth = 1.2;
                        hold on
                        
                        Center_Beams_CR(ind_beam,:) = center_CR';
                    end
                end
                title(['CR:flexible beam pattern result of No.',num2str(iterationNum_CR),' iteration'])
                xlim tight
                xlabel('x,km')
                ylabel('y,km')
                axis equal

                s = 1;
                for ind_beam=1:obj.K
                    if Nuser_beam_CR(ind_beam)== 0                      
                        idBeamOverload = find(Nuser_beam_CR > 45);
                        [~,I] = sort(Nuser_beam_CR(idBeamOverload),'descend');
                        idBeamOverload = idBeamOverload(I);
                        Center_Beams_CR(ind_beam,:) = 0.5*(Center_Beams_CR(idBeamOverload(s))+Center_Beams_CR(idBeamOverload(s+1)));
                        s = s+1;
                    end
                end


                % Compute distances between  user locations and beam centers for Bessel modeling
                distance_All_CR=zeros(obj.N,obj.K);
                for i=1:obj.K
                    distance_All_CR(:,i)= vecnorm(obj.users_locations-Center_Beams_CR(i,:),2,2);
                end
                R = 50;
                u=2.07123*distance_All_CR/R;
                indZ=find(~distance_All_CR);
                % 在贝塞尔函数辐射方向图模型下，从b-th波束到中心波束内n-th用户的信道增益可以表示为
                g= besselj(1,u)./(2*u)+ 36*besselj(3,u)./(u.^3);
                G_dB=10*log10( abs(besselj(1,u)./(2*u)+ 36*besselj(3,u)./(u.^3)).^2  );
                G_dB(indZ)=zeros(1,length(indZ));
                %%% Magnitude Channel matrix, normalized respect to the Noise with bandwidth W
                G_max_Sat= 52;
                Gamma_dB_wo_P=obj.L_SNR_dB+(G_dB+G_max_Sat); %  in dB
                Gamma_wo_P_CR=10.^(Gamma_dB_wo_P/10); % in natural
                % Carrier SNR for uniform resource allocation
                snr_car_CR= Gamma_wo_P_CR*obj.Ref_Pb_Sat/0.5; % Carrier SNR
                SNR_car_CR=10*log10(snr_car_CR); % 272*6
                Th_car= 8.7000; % SNR threshold to obtain a C/I=24 dB or higher.
                ind_snr=find(SNR_car_CR<Th_car);
                snr_car_CR(ind_snr)=0;

                for index_b = 1:obj.K
                    to_average_user_load(index_b)= abs(45-Nuser_beam_CR(index_b))/45;
                end
                loadBalancing = sum(to_average_user_load);
                Assig_UM_FlexB_last = Assig_UM_FlexB;
                Nuser_beam_CR_last = Nuser_beam_CR;
                if loadBalancing<loadBalancing_0
                    flag_loadBalancing = 0;
                    
                else
                    flag_loadBalancing = 1;
                    Assig_UM_FlexB = Assig_UM_FlexB_last;
                    Nuser_beam_CR = Nuser_beam_CR_last;
                end
                loadBalancing_0 = loadBalancing;

            end
toc
            %% Second STEP: User-carrier assignment within the beam

            R_off_CR=zeros(obj.N,1);  % Variable with the offerd traffic to each user

            for ind_beam=1:obj.K % Index of the selected beam
                N_user_beam=Nuser_beam_CR(ind_beam); % Number of user in the selected beam
                M_b=obj.M; % Fixed carrier allocation per beam
                users_index=find(Assig_UM_FlexB==ind_beam); % Auxiliar variable
                Req_user_aux= obj.Req_user(users_index).'; % % Requested user traffic in the selected beam
                if N_user_beam>0
                    % Mixed-integer problem is solved with CVX and Mosek
                    R_off_users= SecondStep_CR(obj,N_user_beam,M_b,ind_beam,users_index,Req_user_aux, obj.Ref_Pb_Sat, 0.5,Gamma_wo_P_CR);
                    R_off_CR(users_index)=R_off_users.';
                end
            end

        end


        %% BW
        function [R_off,M_beam] = FlexibleBandwidth(obj)
            % OUTPUT:
            %       R_off: Vector with the offered user rate with the simulated strategy 提供给波束服务的用户的速率
            %       M_beam: Vector with the number of allocated carriers per beam 载波数

            %% FIRST STEP: Beam resource allocation
            % 第一步:波束资源分配

            % Define convex optimization problem with CVX
            cvx_solver default
            cvx_begin quiet
            variable w_beam(obj.K) % Number beams
            expression R_off_beam(obj.K)
            expression obj_fun

            % C3：避免在相邻波束中使用相同部分的频谱以减少同频干扰
            % 二维与一维不同
                C_OverB_s=zeros(obj.K-1,obj.K);
                v_OverB_s= 2*obj.M*ones(obj.K-1,1);
                for aux_C=1:obj.K-1
                    C_OverB_s(aux_C,:)= [  zeros(1, aux_C-1)  1  1 zeros(1,obj.K-aux_C-1)];
                end
                for aux_C=9:8:obj.K-1
                    C_OverB_s(aux_C,:) = zeros(1,obj.K);
                end

            % Offered traffic per beam
            for index=1:obj.K
                if obj.effective_snr_beam(index)~=0
                    % 提供的速率Roff
                    R_off_beam(index)=w_beam(index)* log( 1+ obj.effective_snr_beam(index)/obj.K   )/log(2);
                else
                    R_off_beam(index)=0;
                end
            end

            % Objective function, Quadratic unment demand (QU)
            % 二次未满足的目标函数
            obj_fun=  (obj.Req_beam.'-R_off_beam).^2; % QU
            % Minimization problem
            minimize( sum(obj_fun) );

            % Apply optimization constraints
            subject to
            w_beam>=0
            C_OverB_s*w_beam<=v_OverB_s

            cvx_end

            % Obtain discrete carrier allocation for each beam from the continuos bandwidth from the first step( w_beam).
            % 从第一步的连续带宽中得到每个波束的离散载波分配
            M_beam=BandwidthDisc(obj,w_beam,1);

            %%%%%%%%%%%% Second STEP: User-carrier assignment within the beam
            % 第二步:波束内用户载波分配

            R_off=zeros(size(obj.Req_user)); % Variable with the offerd traffic to each user 提供根据提供给每个用户的流量而变化

            for ind_beam=1:obj.K   % Beam selection
                N_user_beam=obj.Nuser_beam(ind_beam); % Number of users  in the selected beam
                M_b=M_beam(ind_beam); % Number of carriers in the selected beam
                users_index=obj.users_b_index{ind_beam}; % Auxiliar variable
                Req_user_aux= obj.Req_user(users_index).'; % Requested user traffic in the selected beam
                if N_user_beam>0
                    % Mixed-integer problem is solved with CVX and Mosek
                    R_off_users= SecondStep(obj,N_user_beam,M_b,ind_beam,users_index,Req_user_aux, obj.Ref_Pb_Sat, 0.5);
                    R_off(users_index)=R_off_users.';
                end
            end


        end


        %% POW
        function [R_off,P_flex] = FlexiblePower(obj)
            %       R_off: Vector with the offered user rate with the simulated strategy
            %       P_flex: Vector with the allocated power per beam.

            %%%%%%%%%%%% FIRST STEP: Beam resource allocation

            % Optimization of power per HPA
            Nvar=obj.K/2; % One amplifier per pair of beams.
            p0=zeros(1,Nvar);
            ub=obj.Max_Pamp/obj.Max_P_sat*ones(1,Nvar); % Upper bound for the power allocation
            lb=zeros(1,Nvar); % Lower bound for the power allocation

            % Settings for the optimization
            alg='sqp';
            options_FlexP_Beam = optimoptions('fmincon','OptimalityTolerance',1e-10,'StepTolerance',1e-10,'ConstraintTolerance',1e-10,'MaxFunctionEvaluations',4e5,'Display','off','Diagnostics','off','Algorithm',alg);
            A=ones(1,Nvar);
            b=1;

            % Optimization ( with fmincon)
            [P_flex_aux,~,~,~] =  fmincon( @(p) Flex_Power_Beam(p,obj.Req_beam,obj.K,obj.effective_snr_beam),p0,A,b,[],[],lb,ub,[],options_FlexP_Beam);

            % Transalte power per HPA to power per beam
            P_flex=obj.Max_P_sat*repelem(P_flex_aux,2)/2; % Power per beam ( Same power for each beam sharing the HPA)


            %%%%%%%%%%%% Second STEP: User-carrier assignment within the beam
            R_off=zeros(size(obj.Req_user)); % Variable with the offerd traffic to each user

            for index_beam=1:obj.K % Index of the selected beam
                N_user_beam=obj.Nuser_beam(index_beam); % Number of users  in the selected beam
                users_index=obj.users_b_index{index_beam};  % Auxiliar variable
                M_b=obj.M; % Fixed number of carriers per beam
                Req_user_aux=obj.Req_user(users_index).'; % Requested user traffic in the selected beam
                if N_user_beam>0
                    % Mixed-integer problem is solved with CVX and Mosek
                    R_off_users= SecondStep(obj,N_user_beam,M_b,index_beam,users_index,Req_user_aux, P_flex(index_beam), 0.5);
                    R_off(users_index)=R_off_users.';
                end
            end


        end


        %% MAP
        function [R_off,Nuser_beam_UM_FixedRes]=FixResFlexibleMapping(obj)
            %       R_off: Vector with the offered user rate with the simulated strategy
            %       Nuser_beam_UM_FixedRes: Vector that contains the association between users and beams.

            % Relaxed problem to obtain the beam-user assignment
            cvx_solver default
            cvx_begin quiet
            variable w_beam_user(obj.K,obj.N)
            expression Roff_sim2(obj.K,obj.N)
            expression obj_fun
            for index_beam=1:obj.K
                for index_user=1:obj.N
                    Roff_sim2( index_beam,index_user)= ...
                        w_beam_user(index_beam,index_user)*obj.Delta_W* ...
                        log(1 + obj.snr_car(index_user,index_beam) )/log(2);
                end
            end
            obj_fun= ( obj.Req_user.' - sum(Roff_sim2,1)).^2;
            minimize( sum(obj_fun) );
            subject to
            0<=w_beam_user<=1
            sum(w_beam_user,1)<= 1 %  Constraint of one carrier per user
            sum(w_beam_user,2)<=obj.M % Constraint of fixed bandwidth ( M carriers) per beam
            cvx_end

            % Number of user per beam
            Assig_UM_FlexB=zeros(1,obj.N);
            for ind_n=1:obj.N
                [~,indmax]=max( Roff_sim2(:,ind_n) );
                Assig_UM_FlexB(ind_n)= indmax;
            end

            % Obtain association between users and beams.柔性映射
            Nuser_beam_UM_FixedRes=hist(Assig_UM_FlexB,1:obj.K);

            %% Second STEP: User-carrier assignment within the beam

            R_off=zeros(obj.N,1);  % Variable with the offerd traffic to each user

            for ind_beam=1:obj.K % Index of the selected beam
                N_user_beam=Nuser_beam_UM_FixedRes(ind_beam); % Number of user in the selected beam
                M_b=obj.M; % Fixed carrier allocation per beam
                users_index=find(Assig_UM_FlexB==ind_beam); % Auxiliar variable
                Req_user_aux= obj.Req_user(users_index).'; % % Requested user traffic in the selected beam
                if N_user_beam>0
                    % Mixed-integer problem is solved with CVX and Mosek
                    R_off_users= SecondStep(obj,N_user_beam,M_b,ind_beam,users_index,Req_user_aux, obj.Ref_Pb_Sat, 0.5);
                    R_off(users_index)=R_off_users.';
                end
            end


        end


        %% CR-BW
        function  [R_off,Nuser_beam_UM_FlexB,M_beam,iterationNum_CRBW]=FlexBandwidthFlexibleCoverRange(obj)
            % 定义初始负载均衡度 loadBalancing_0
            for index_b = 1:obj.K
                to_average_user_load(index_b)= abs(45-obj.Nuser_beam(index_b))/45;
            end
            loadBalancing_0 = sum(to_average_user_load);
            % 初始信噪比
            snr_car_CR = obj.snr_car;
            flag_loadBalancing = 0;
            % 当负载均衡度没有比上一次负载均衡度小时，跳出循环
            iterationNum_CRBW = 0;
            while flag_loadBalancing == 0
                iterationNum_CRBW = iterationNum_CRBW + 1;
                % Relaxed problem to obtain the beam-user assignment and bandwidth allocation
                C_OverB_sim=zeros(obj.K-1,obj.K);
                v_OverB_sim= 2*obj.M*ones(obj.K-1,1);
                for aux_C=1:obj.K-1
                    C_OverB_sim(aux_C,:)= [  zeros(1, aux_C-1)  1  1 zeros(1,obj.K-aux_C-1)];
                end
                for aux_C=9:8:obj.K-1
                    C_OverB_sim(aux_C,:) = zeros(1,obj.K);
                end
                cvx_solver default
                cvx_begin quiet
                variable w_sim(obj.K,obj.N)
                expression Roff_sim(obj.K,obj.N)
                expression obj_fun
                for index_beam=1:obj.K
                    for index_user=1:obj.N
                        Roff_sim( index_beam,index_user)= ...
                            w_sim(index_beam,index_user)*obj.Delta_W* ...
                            log(1 +  snr_car_CR(index_user,index_beam)   )/log(2);
                    end
                end
                obj_fun= ( obj.Req_user.' - sum(Roff_sim,1)).^2;
                minimize( sum(obj_fun) );
                subject to
                0<=w_sim<=1
                sum(w_sim,1)<= 1 %  Constraint of one carrier per user
                % 错误使用  * Inner matrix dimensions must agree.
                C_OverB_sim*sum(w_sim,2)<=v_OverB_sim % Constraint to avoid the reuse of bandwith in adjacent beams
                cvx_end

                % Extract beam user mapping from the results of the relaxed problem. If no rates are assigned to a user, the user is assigned to its dominant beam
                Assig_UM_FixedRes=zeros(1,obj.N);
                for ind_n=1:obj.N
                    [vmax,indmax]=max( Roff_sim(:,ind_n) );
                    if vmax<obj.tol
                        Assig_UM_FixedRes(ind_n)=obj.user_beams(ind_n);
                    else
                        Assig_UM_FixedRes(ind_n)= indmax;
                    end
                end

                % Obtain beam-user mapping
                Nuser_beam_UM_FlexB=hist(Assig_UM_FixedRes,1:obj.K);

                % Obtain discrete carrier allocation for each beam from the continuous bandwidth from the first step( w_beam).
                w_beam=sum(w_sim,2)+0;
                M_beam=BandwidthDisc_CR(obj,w_beam,2,Nuser_beam_UM_FlexB);

                % 改变波束中心，画出波束覆盖图
                figure('name','BW-CR-alternate')
                Center_Beams_CR = zeros(obj.K,2);
                for ind_beam=1:obj.K
                    if Nuser_beam_UM_FlexB(ind_beam)> 0
                        users_index=find(Assig_UM_FixedRes==ind_beam);
                        x_CR = obj.users_locations(users_index,1);
                        y_CR = obj.users_locations(users_index,2);
                        scatter(x_CR,y_CR,...
                            6,'o','MarkerEdgeColor','none','MarkerFaceColor',"#DD7C4F");
                        hold on

                        [center_CR, radius_CR] = min_cover_circle(x_CR, y_CR, Nuser_beam_UM_FlexB(ind_beam));

                        pos = [center_CR(1)-radius_CR, ...
                            center_CR(2)-radius_CR, ...
                            2*radius_CR, 2*radius_CR];
                        r1 = rectangle('Position',pos,'Curvature',[1 1]);
                        r1.LineStyle = '-';
                        r1.EdgeColor = '#6F6F6F';
                        r1.LineWidth = 1.2;
                        hold on

                        Center_Beams_CR(ind_beam,:) = center_CR';
                    end
                end
                title(['BW-CR:flexible beam pattern result of No.',num2str(iterationNum_CRBW),' iteration'])
                xlim tight
                xlabel('x,km')
                ylabel('y,km')
                axis equal

                s = 1;
                for ind_beam=1:obj.K
                    if Nuser_beam_UM_FlexB(ind_beam)== 0
                        idBeamOverload = find(Nuser_beam_UM_FlexB > 45);
                        [~,I] = sort(Nuser_beam_UM_FlexB(idBeamOverload),'descend');
                        idBeamOverload = idBeamOverload(I);
                        Center_Beams_CR(ind_beam,:) = 0.5*(Center_Beams_CR(idBeamOverload(s))+Center_Beams_CR(idBeamOverload(s+1)));
                        s = s+1;
                    end
                end

                % Compute distances between  user locations and beam centers for Bessel modeling
                distance_All_CR=zeros(obj.N,obj.K);
                for i=1:obj.K
                    distance_All_CR(:,i)= vecnorm(obj.users_locations-Center_Beams_CR(i,:),2,2);
                end
                R = 50;
                u=2.07123*distance_All_CR/R;
                indZ=find(~distance_All_CR);
                % 在贝塞尔函数辐射方向图模型下，从b-th波束到中心波束内n-th用户的信道增益可以表示为
                g= besselj(1,u)./(2*u)+ 36*besselj(3,u)./(u.^3);
                G_dB=10*log10( abs(besselj(1,u)./(2*u)+ 36*besselj(3,u)./(u.^3)).^2  );
                G_dB(indZ)=zeros(1,length(indZ));
                %%% Magnitude Channel matrix, normalized respect to the Noise with bandwidth W
                G_max_Sat= 52;
                Gamma_dB_wo_P=obj.L_SNR_dB+(G_dB+G_max_Sat); %  in dB
                Gamma_wo_P_CR=10.^(Gamma_dB_wo_P/10); % in natural
                % Carrier SNR for uniform resource allocation
                snr_car_CR= Gamma_wo_P_CR*obj.Ref_Pb_Sat/0.5; % Carrier SNR
                SNR_car_CR=10*log10(snr_car_CR); % 272*6
                Th_car= 8.7000; % SNR threshold to obtain a C/I=24 dB or higher.
                ind_snr=find(SNR_car_CR<Th_car);
                snr_car_CR(ind_snr)=0;

                for index_b = 1:obj.K
                    to_average_user_load(index_b)= abs(45-Nuser_beam_UM_FlexB(index_b))/45;
                end
                loadBalancing = sum(to_average_user_load);
                Assig_UM_FixedRes_last = Assig_UM_FixedRes;
                Nuser_beam_UM_FlexB_last = Nuser_beam_UM_FlexB;
                if loadBalancing<loadBalancing_0
                    flag_loadBalancing = 0;
                else
                    flag_loadBalancing = 1;
                    Assig_UM_FixedRes = Assig_UM_FixedRes_last;
                    Nuser_beam_UM_FlexB = Nuser_beam_UM_FlexB_last;
                end
                loadBalancing_0 = loadBalancing;

            end

            %%%%%%%%%%%% Second STEP: User-carrier assignment within the beam
            R_off=zeros(obj.N,1);  % Variable with the offerd traffic to each user

            for ind_beam=1:obj.K % Index of the selected beam
                N_user_beam=Nuser_beam_UM_FlexB(ind_beam); % Number of users  in the selected beam
                M_b=M_beam(ind_beam); % Number of carriers  in the selected beam
                users_index=find(Assig_UM_FixedRes==ind_beam); % Auxiliar variable
                Req_user_aux= obj.Req_user(users_index).'; % Requested uset traffic in the selected beam
                if N_user_beam>0
                    % Mixed-integer problem is solved with CVX and Mosek
                    R_off_users= SecondStep_CR(obj,N_user_beam,M_b,ind_beam,users_index,Req_user_aux, obj.Ref_Pb_Sat, 0.5,Gamma_wo_P_CR);
                    R_off(users_index)=R_off_users.';
                end
            end
        end
        

        %% BW-MAP
        function  [R_off,Nuser_beam_UM_FlexB,M_beam]=FlexBandwidthFlexMapping(obj)
            %       R_off: Vector with the offered user rate with the simulated strategy
            %       Nuser_beam_UM_FlexB: Vector that contains the association between users and beams
            %       M_beam: Vector with the number of allocated carriers per beam

            % Relaxed problem to obtain the beam-user assignment and bandwidth allocation
            C_OverB_sim=zeros(obj.K-1,obj.K);
            v_OverB_sim= 2*obj.M*ones(obj.K-1,1);
            for aux_C=1:obj.K-1
                C_OverB_sim(aux_C,:)= [  zeros(1, aux_C-1)  1  1 zeros(1,obj.K-aux_C-1)];
            end
            for aux_C=9:8:obj.K-1
                C_OverB_sim(aux_C,:) = zeros(1,obj.K);
            end

            cvx_solver default
            cvx_begin quiet
            variable w_sim(obj.K,obj.N)
            expression Roff_sim(obj.K,obj.N)
            expression obj_fun
            for index_beam=1:obj.K
                for index_user=1:obj.N
                    Roff_sim( index_beam,index_user)=w_sim(index_beam,index_user)*obj.Delta_W* log(1 +  obj.snr_car(index_user,index_beam)   )/log(2);
                end
            end
            obj_fun= ( obj.Req_user.' - sum(Roff_sim,1)).^2;
            minimize( sum(obj_fun) );
            subject to
            0<=w_sim<=1
            sum(w_sim,1)<= 1 %  Constraint of one carrier per user
            C_OverB_sim*sum(w_sim,2)<=v_OverB_sim % Constraint to avoid the reuse of bandwith in adjacent beams
            cvx_end


            % Extract beam user mapping from the results of the relaxed problem. If no rates are assigned to a user, the user is assigned to its dominant beam
            Assig_UM_FixedRes=zeros(1,obj.N);
            for ind_n=1:obj.N

                [vmax,indmax]=max( Roff_sim(:,ind_n) );
                if vmax<obj.tol
                    Assig_UM_FixedRes(ind_n)=obj.user_beams(ind_n);
                else
                    Assig_UM_FixedRes(ind_n)= indmax;
                end
            end

            % Obtain discrete carrier allocation for each beam from the continuous bandwidth from the first step( w_beam).
            w_beam=sum(w_sim,2)+0;
            M_beam=BandwidthDisc(obj,w_beam,2);

            % Obtain beam-user mapping
            Nuser_beam_UM_FlexB=hist(Assig_UM_FixedRes,1:obj.K);


            %%%%%%%%%%%% Second STEP: User-carrier assignment within the beam
            R_off=zeros(obj.N,1);  % Variable with the offerd traffic to each user

            for ind_beam=1:obj.K % Index of the selected beam
                N_user_beam=Nuser_beam_UM_FlexB(ind_beam); % Number of users  in the selected beam
                M_b=M_beam(ind_beam); % Number of carriers  in the selected beam
                users_index=find(Assig_UM_FixedRes==ind_beam); % Auxiliar variable
                Req_user_aux= obj.Req_user(users_index).'; % Requested uset traffic in the selected beam
                if N_user_beam>0
                    % Mixed-integer problem is solved with CVX and Mosek
                    R_off_users= SecondStep(obj,N_user_beam,M_b,ind_beam,users_index,Req_user_aux, obj.Ref_Pb_Sat, 0.5);
                    R_off(users_index)=R_off_users.';
                end
            end
        end


        %% BW-POW
        function  [R_off,M_GA,P_GA]=FlexBandwidthPower(obj,P_TwoStep,M_beam,optionsGA)

            % INPUT:
            %       obj: The class object that contains the scenario information
            % OUTPUT:
            %       R_off: Vector with the offered user rate with the simulated strategy
            %       M_GA: Vector with the number of allocated carriers per beam
            %       M_beam: Vector with the allocated power per beam.

            %%%%%%%%%%%% FIRST STEP: Beam resource allocation
            % Genetic algorithm settings
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            aux_P=1:obj.K/2;  % Auxiliar variable
            aux_B=obj.K/2+1:obj.K+obj.K/2; % Auxiliar variable
            lb_P = zeros(1,obj.K/2); % Lower bound for the power allocation
            ub_P =obj.Max_P_sat/(obj.K*obj.M)*2*ones(1,obj.K/2);  % Upper bound for the power allocation
             lb_B = zeros(1,obj.K);  % Lower bound for the bandwidth allocation
            lb= [lb_P lb_B ];  % Lower bound for the genetic algorithm
            ub= [ub_P obj.ub_B]; % Upper bound for the genetic algorithm
            aux_pflex_car= P_TwoStep(1:2:obj.K)/obj.M; % Auxiliar variable. Power allocation for the fixed bandwidth and flexible power assignment.
            un_case=  [ obj.Max_P_sat/(obj.K*obj.M)*ones(1,obj.K/2) obj.M*ones(1,obj.K) ; % Initial individuals for the genetic algorithm. Solutions from other resource assigments
                aux_pflex_car obj.M*ones(1,obj.K);                                        % are employed as individuals.
                obj.Max_P_sat/(obj.K*obj.M)*ones(1,obj.K/2)   M_beam' ];
            options = optimoptions('ga','Display','none','InitialPopulationMatrix',un_case,'MaxGenerations',optionsGA.Max_Gen,'UseParallel', false,'CreationFcn','gacreationuniform','EliteCount',optionsGA.Elite_num,'MutationFcn', {@mutationadaptfeasible,optionsGA.pmut},'PopulationSize',optionsGA.PoP,'SelectionFcn', {@selectiontournament,5});


            % Genetic algorithm
            [x_ga,~,~,~,~,~] = ga(@(x) ga_function(x,obj),length(lb),A,b,Aeq,beq,lb,ub,[],4:9,options);


            %  Process the results from the GA from MATLAB ( Cheking of the values is needed due to the GA implementation in MATLAB)
            % Check Bandwidth
            [B] = Check_B(x_ga(aux_B), obj);
            % Check Power
            P=repelem(x_ga(aux_P),2).*B*2*obj.M;
            P_amp=zeros(1,obj.K/2);
            for index=1:obj.K/2
                P_amp(index)=sum(P( 2*(index-1)+1:2*index));

                if    P_amp(index)> obj.Max_Pamp
                    k=  obj.Max_Pamp/P_amp(index);
                    P( 2*(index-1)+1:2*index)= k*P( 2*(index-1)+1:2*index);
                end
            end
            if sum(P_amp)>obj.Max_P_sat
                k=  obj.Max_P_sat/sum(P_amp);
                P= k*P;
            end

            % Final resource assignment
            P_GA=P;
            M_GA=B*2*obj.M;


            %%%%%%%%%%%% Second STEP: User-carrier assignment within the beam
            R_off=zeros(obj.N,1);  % Variable with the offerd traffic to each user
            for index_beam=1:obj.K
                N_user_beam=obj.Nuser_beam(index_beam); % Number of users  in the selected beam
                users_index=obj.users_b_index{index_beam}; % Auxiliar variable
                Req_user_aux=obj.Req_user(users_index).'; % Requested user traffic in the selected beam
                M_b=M_GA(index_beam); % Number of carriers in the selected beam
                %%%% Opt problem with explicit carriers
                if N_user_beam>0
                    % Mixed-integer problem is solved with CVX and Mosek
                    R_off_users= SecondStep(obj,N_user_beam,M_b,index_beam,users_index,Req_user_aux, P_GA(index_beam), M_b*obj.Delta_W);
                    R_off(users_index)=R_off_users.';
                end
            end


        end


    end


    methods (Access=protected )

        function effective_snr=ComputEffectiveSNR (obj)
            % Compute effective snr per beam
            % OUTPUT:
            %       effective_snr: Vector with effective snr per beam

            effective_snr=zeros(1,obj.K);
            for index=1:obj.K
                users_index=obj.users_b_index{index};
                if ~isempty(users_index)
                    aux=obj.Gamma_wo_P( users_index,index)*obj.Max_P_sat/0.5  ;
                    effective_snr(index)=  geo_mean(aux) ;
                end
            end

        end
        
        
        %% CR
        function effective_snr=ComputEffectiveSNR_CR (obj,Gamma_wo_P_CR,users_b_index_CR)
            % Compute effective snr per beam
            % OUTPUT:
            %       effective_snr: Vector with effective snr per beam

            effective_snr=zeros(1,obj.K);
            for index=1:obj.K
                users_index=users_b_index_CR{index};
                if ~isempty(users_index)
                    aux=Gamma_wo_P_CR( users_index,index)*obj.Max_P_sat/0.5  ;
                    effective_snr(index)=  geo_mean(aux) ;
                end
            end

        end


        %%
        function M_beam= BandwidthDisc(obj,W_cont,mode)
            % Obtain discrete carrier allocation for each beam from the obtained results( W_cont).
            % 从得到的结果(W_cont)中得到每个波束的离散载波分配。
            % If two beams reuse the same carrier after the discretization, the conflict is resolve by allocating the carrier to the beam with higher traffic demand.
            % 如果两个波束在离散化后重用同一载波，则通过将载波分配给业务量需求较大的波束来解决冲突。
            % INPUT:
            %       obj: The class object that contains the scenario information
            %       W_cont: Continous bandwidth allocation
            %       mode: Value 1 when discretizing the bandwidh from the solution with flexible bandwidth allocation.
            %             Value 2 when discretizing the bandwidh from the solution with flexible bandwidth and beam-user mapping.
            % OUTPUT:
            %       M_beam: Discrete bandwidth. Vector with the number of carriers per beam
            M_beam = zeros(obj.K,1);
            % Process the input
            switch mode
                case 1 % Flexible Bandwidth
                    M_beam_aux= W_cont*obj.M/0.5;
                    M_beam_aux( obj.Nuser_beam==0)=0;
                case 2 % Flexible Bandwidth. Flexibe Mapping
                    M_beam_aux=W_cont;
            end

            for col = 1 : 8
                idcol = (col-1)*8+1;
                % First step of discretization
                M_beamPerCol=floor(M_beam_aux(idcol:idcol+7));
                % Perform a floor operation as first step for the discretization.
                Number_carrier_conflict =  obj.M*8-sum(M_beamPerCol);
                % Check if more carriers can be allocated or the first step violates the allowed numbers of carriers.
                [~,indS]=sort(M_beam_aux(idcol:idcol+7)-M_beamPerCol,'descend');
                % Order the beams regarding the unsastiefied bandwidth


                % Second step of the discretization.
                % Carriers area allocated following the previous sorting order.
                % The carrier allocation takes into account the bandwidth of the adjacent beam.
                aux_cont=0; % Auxiliar variable.
                Nuser_beamPerCol = obj.Nuser_beam(idcol:idcol+7);
                for index=1:8
                    if Nuser_beamPerCol(indS(index))>0
                        if aux_cont==Number_carrier_conflict
                            break
                        end
                        switch  indS(index)
                            case 1
                                if  (M_beamPerCol(indS(index)+1)+ M_beamPerCol(indS(index))+1)<=2*obj.M
                                    M_beamPerCol(indS(index))= M_beamPerCol(indS(index))+1;
                                    aux_cont=aux_cont+1;
                                end
                            case 8
                                if  (M_beamPerCol(indS(index)-1)+ M_beamPerCol(indS(index))+1)<=2*obj.M
                                    M_beamPerCol(indS(index))= M_beamPerCol(indS(index))+1;
                                    aux_cont=aux_cont+1;
                                end
                            otherwise
                                if  (M_beamPerCol(indS(index)-1)+ M_beamPerCol(indS(index))+1 )<=2*obj.M && (M_beamPerCol(indS(index)+1)+ M_beamPerCol(indS(index))+1)<=2*obj.M
                                    M_beamPerCol(indS(index))= M_beamPerCol(indS(index))+1;
                                    aux_cont=aux_cont+1;
                                end
                        end
                    end
                end
                M_beam(idcol:idcol+7) = M_beamPerCol;
            end
        end


        %% CR
        function M_beam= BandwidthDisc_CR(obj,W_cont,mode,Nuser_beam_UM_FlexB)
            % Obtain discrete carrier allocation for each beam from the obtained results( W_cont).
            % 从得到的结果(W_cont)中得到每个波束的离散载波分配。
            % If two beams reuse the same carrier after the discretization, the conflict is resolve by allocating the carrier to the beam with higher traffic demand.
            % 如果两个波束在离散化后重用同一载波，则通过将载波分配给业务量需求较大的波束来解决冲突。
            % INPUT:
            %       obj: The class object that contains the scenario information
            %       W_cont: Continous bandwidth allocation
            %       mode: Value 1 when discretizing the bandwidh from the solution with flexible bandwidth allocation.
            %             Value 2 when discretizing the bandwidh from the solution with flexible bandwidth and beam-user mapping.
            % OUTPUT:
            %       M_beam: Discrete bandwidth. Vector with the number of carriers per beam
            
            M_beam = zeros(obj.K,1);
            % Process the input
            switch mode
                case 1 % Flexible Bandwidth
                    M_beam_aux= W_cont*obj.M/0.5;
                    M_beam_aux( Nuser_beam_UM_FlexB==0)=0;
                case 2 % Flexible Bandwidth. Flexibe Mapping
                    M_beam_aux=W_cont;
            end

            for col = 1 : 8
                idcol = (col-1)*8+1;
                % First step of discretization
                M_beamPerCol=floor(M_beam_aux(idcol:idcol+7));
                % Perform a floor operation as first step for the discretization.
                Number_carrier_conflict =  obj.M*8-sum(M_beamPerCol);

                % Check if more carriers can be allocated or the first step violates the allowed numbers of carriers.
                [~,indS]=sort(M_beam_aux(idcol:idcol+7)-M_beamPerCol,'descend');
                % Order the beams regarding the unsastiefied bandwidth


                % Second step of the discretization.
                % Carriers area allocated following the previous sorting order.
                % The carrier allocation takes into account the bandwidth of the adjacent beam.
                aux_cont=0; % Auxiliar variable.
                Nuser_beamPerCol = Nuser_beam_UM_FlexB(idcol:idcol+7);
                for index=1:8
                    if Nuser_beamPerCol(indS(index))>0
                        if aux_cont==Number_carrier_conflict
                            break
                        end
                        switch  indS(index)
                            case 1
                                if  (M_beamPerCol(indS(index)+1)+ M_beamPerCol(indS(index))+1)<=2*obj.M
                                    M_beamPerCol(indS(index))= M_beamPerCol(indS(index))+1;
                                    aux_cont=aux_cont+1;
                                end
                            case 8
                                if  (M_beamPerCol(indS(index)-1)+ M_beamPerCol(indS(index))+1)<=2*obj.M
                                    M_beamPerCol(indS(index))= M_beamPerCol(indS(index))+1;
                                    aux_cont=aux_cont+1;
                                end
                            otherwise
                                if  (M_beamPerCol(indS(index)-1)+ M_beamPerCol(indS(index))+1 )<=2*obj.M && (M_beamPerCol(indS(index)+1)+ M_beamPerCol(indS(index))+1)<=2*obj.M
                                    M_beamPerCol(indS(index))= M_beamPerCol(indS(index))+1;
                                    aux_cont=aux_cont+1;
                                end
                        end
                    end

                end
                M_beam(idcol:idcol+7) = M_beamPerCol;
            end

        end


        %% CR
        function R_off_users=SecondStep_CR(obj,N,M,ind_beam,users_index,Req_user_aux,P,B,Gamma_wo_P_CR)
            % User-Carrier assignment within the beams

            %       N: Number of users served by the selected beam
            %       M: Number of allocated carriers to the selected beam
            %       users_index: Auxiliar variable. It associates the users to the beam
            %       Req_user_aux: Vector with the user demanded rates
            %       P: Vector with the power per beam
            %       B: Vector with the bandwidth per beam ( Normalized to the total avaliable bandwidth)

            %       R_off_users: Offered user rate to the user within the selected beam

            % Mixed Binary Quadratic Program (MBQP) problem, solved with CVX and MOSEK
            cvx_solver Mosek
            cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', obj.Max_time_opt_second,'MSK_DPAR_LOWER_OBJ_CUT',-1.0e30);
            cvx_begin quiet
            variable w(M,N)
            variable u(M,N) binary
            expression Roff_cvx(M,N)
            expression obj_fun
            for index_car=1:M
                for index_user=1:N
                    Roff_cvx(index_car,index_user)=w(index_car,index_user)*obj.Delta_W ...
                        * log(1 + Gamma_wo_P_CR(users_index(index_user),ind_beam)*P/B)/log(2);
                end
            end
            obj_fun= ( Req_user_aux  - sum(Roff_cvx,1)).^2; % QU
            minimize( sum(obj_fun) );
            subject to
            u>=w
            sum(u,1)<= 1 % Constraint of one carrier per user
            sum(w,2)<=1 % Constraint of overall carrier time allocation
            0<=w<=1
            cvx_end
            % Store the obtained values
            R_off_users=sum(Roff_cvx,1);

        end


        %%
        function R_off_users=SecondStep(obj,N,M,ind_beam,users_index,Req_user_aux,P,B)
            % User-Carrier assignment within the beams

            %       N: Number of users served by the selected beam
            %       M: Number of allocated carriers to the selected beam
            %       users_index: Auxiliar variable. It associates the users to the beam
            %       Req_user_aux: Vector with the user demanded rates
            %       P: Vector with the power per beam
            %       B: Vector with the bandwidth per beam ( Normalized to the total avaliable bandwidth)
            %       R_off_users: Offered user rate to the user within the selected beam

            % Mixed Binary Quadratic Program (MBQP) problem, solved with CVX and MOSEK
            % cvx_solver Mosek
            cvx_solver Mosek
            cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', obj.Max_time_opt_second,'MSK_DPAR_LOWER_OBJ_CUT',-1.0e30);
            cvx_begin quiet
            variable w(M,N)
            variable u(M,N) binary
            expression Roff_cvx(M,N)
            expression obj_fun
            for index_car=1:M
                for index_user=1:N
                    Roff_cvx(index_car,index_user)=w(index_car,index_user)*obj.Delta_W ...
                        * log(1 + obj.Gamma_wo_P(users_index(index_user),ind_beam)*P/B)/log(2);
                end
            end
            obj_fun= ( Req_user_aux  - sum(Roff_cvx,1)).^2; % QU
            minimize( sum(obj_fun) );
            subject to
            u>=w
            sum(u,1)<= 1 % Constraint of one carrier per user
            sum(w,2)<=1 % Constraint of overall carrier time allocation
            0<=w<=1
            cvx_end
            % Store the obtained values
            R_off_users=sum(Roff_cvx,1);

        end


    end

end


%% Auxiliar functions

function [fval] = Flex_Power_Beam(x,Req_b,K,effective_snr_beam)
% Objective function for the flexible power allocation, for the optimization process with fmincon.
% INPUT:
%       x: Solution input ( Power allocation in this case)
%       Req_b: Requested traffic per beam
%       K: Number of beams
%       effective_snr_beam: Vector with the effective snr per beam
% OUTPUT:
%       fval: Evaluation of the objective function ( Quadratic Unment demand)


% Process the input
x=repelem(x,2)/2;
h=zeros(1,K);

% Measure the quadratic unmnet per beam
Roff=zeros(1,K);
for index=1:K
    if effective_snr_beam(index)~=0
        Roff(index)= 0.5*log2(  1+ effective_snr_beam(index)*x(index)   );
    else
        Roff(index)=0;
    end

    h(index)= ( Req_b(index)-Roff(index))^2;
end

% Aggregate the quadratic unmnet per beam
fval= sum(h);



end


function [f] = ga_function(x,obj )
% Objective function for genetic algorithm.
% INPUT:
%       x: Solution input ( Power allocation in this case)
%       obj: Struct with the simulated scenario information
% OUTPUT:
%       f: Evaluation of the objective function ( Quadratic Unment demand)



tol=1e-8;
nvars=length(x);
% Process the input
aux_P=1:obj.K/2;
aux_B=obj.K/2+1:nvars;

% Cheking of the values is needed due to the GA implementation
% Verify Bandwidth
[B] = Check_B(x(aux_B),obj);
% Carrier power to beam power
P=repelem(x(aux_P),2).*B*2*obj.M;
P_amp=zeros(1,obj.K/2);
for index=1:obj.K/2
    P_amp(index)=sum(P( 2*(index-1)+1:2*index));
    if    P_amp(index)> obj.Max_Pamp
        k=  obj.Max_Pamp/P_amp(index);
        P( 2*(index-1)+1:2*index)= k*P( 2*(index-1)+1:2*index);
    end
end
if sum(P_amp)>obj.Max_P_sat
    k=  obj.Max_P_sat/sum(P_amp);
    P= k*P;
end

% Obtain offered traffic per beam
R_off_beam=zeros(1,obj.K);
for index=1:obj.K

    users_index=find(  obj.user_beams==index);
    if ~isempty(users_index)
        if B(index)~=0
            aux=1+  obj.aux_Cn( users_index,index)*P(index)/(B(index));
        else
            aux=1;
        end
        snr_eff2_b=  geo_mean(aux) ;
    else
        snr_eff2_b=1;
    end
    if B(index)~=0
        R_off_beam(index)=B(index)* log(  snr_eff2_b   )/log(2);
    else
        R_off_beam(index)=0;
    end

end

% Measure the Quadratic Unmet demand
f= sum ( (R_off_beam-obj.Req_beam).^2);

if f<tol
    f=0;
end


end


function [B] = Check_B(B_in,obj)
% Verification and correction of the allocated bandwidth. This is described int Appendix B from as "Flexible User Mapping
% for Radio Resource Assignment in Advanced Satellite Payloads" arXiv, 2021. Is based on the process described in " A Genetic
% Algorithm for Joint Power and Bandwidth Allocation in Multibeam Satellite Systems", IEE Aeroespace Conference, 2019
% INPUT:
%       B_in: Allocated bandwidth
%       obj: Class Object
% OUTPUT:
%       B:  Allocated bandwidth after verification
B = zeros(1,64);
for col = 1 : 8
    % Process the input
    idcol = (col-1)*8+1;
    B_inpercol = B_in(idcol:idcol+7);
    Bpercol=B_inpercol/(2*obj.M);

    % Selects randomly a processing order
    coin_flip=rand(1);

    if coin_flip>0.5
        order_lim=1:length(Bpercol)-1;
        step=1;
    else
        order_lim=length(Bpercol):-1:2;
        step=-1;
    end


    % Ensure no bandwidth overlaping in neighbour beams
    for ind=order_lim

        if (Bpercol(ind)+Bpercol(ind+step))>1

            Bpercol(ind)=1-Bpercol(ind+step);
        end

    end


    % Ensure that the all avaliable bandwidth is employed
    Req_beampercol = obj.Req_beam(idcol:idcol+7);
    [~,order_req_Beam]=sort(Req_beampercol,'descend');

    for ind_center_b=order_req_Beam

        ind_left=ind_center_b-1;
        ind_right=ind_center_b+1;

        B_center=Bpercol(ind_center_b);

        if ind_left>0
            B_left=Bpercol(ind_left);
        else
            B_left=0;
        end

        if ind_right<=length(Bpercol)
            B_right=Bpercol(ind_right);
        else
            B_right=0;
        end


        B_un= 1- B_center -max(B_right,B_left);

        if  (B_center+B_un)> obj.ub_B(ind_center_b)/(2*obj.M)

            Bpercol(ind_center_b)= obj.ub_B(ind_center_b)/(2*obj.M);

        else
            Bpercol(ind_center_b)= B_center+B_un;
        end


        
    end
    B(idcol:idcol+7) = Bpercol;
end


end
