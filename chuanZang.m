% 用川藏交界处553600平方公里区域人口密度分布图转成相应密度2897个用户的散点图
% load("HuaBei_pop.mat");
% load("ChuanZang_pop.mat");
load("ChuanZang_pop1.mat");
figure("Name","Population density map of part of Sichuan-Tibet border area，692*800")
imshow(I1)
I11 = I1(sort(randperm(1003, 800)),:,:);
I12 = I11(:,sort(randperm(756, 690)),:);
figure("Name","Population density map")
imshow(I12)
I1 = I12;

I2 = rgb2gray(I1);
I3 = mat2gray(I2);
% figure(2)
% imshow(I3)
index = find(I3);
colorNot0 = I3(index);
% colorNot0_sort = sort(colorNot0);
% colorHist=histogram(colorNot0,7);
% 显示向量中数据的频数表。对于每个唯一值，显示该值的实例数和百分比
table = tabulate(colorNot0);
% 不同颜色对应的人口数量
% legend = [130000 5000 2500 500 100 50 25 5];
legend = [40000 2601 601 201 51 26 6 1];
% 真实总人口数
num_real_alluser = legend*table(:,2);
% 需要的散点数
num_spot = 2897;
% 真实总人口数/需要的散点数
rate_to_spot = num_real_alluser/num_spot;
% 每种颜色所占的真实人口数
for i = 1:8
    numPerLG = table(:,2);
    num_realuser_eachlabel(i) = legend(i)*numPerLG(i);
end
% 每种颜色所占的散点数
num_spot_eachlabel = round(num_realuser_eachlabel/rate_to_spot);
% num_spot_eachlabel(1) = num_spot_eachlabel(1)-1;
figure("Name","Generate Users Based on Population Density Map")
spotLegend1 = find(I3==table(1,1));
spotLegend2 = find(I3==table(2,1));
spotLegend3 = find(I3==table(3,1));
spotLegend4 = find(I3==table(4,1));
spotLegend5 = find(I3==table(5,1));
spotLegend6 = find(I3==table(6,1));
spotLegend7 = find(I3==table(7,1));
spotLegend8 = find(I3==table(8,1));

ispotLegend1 = spotLegend1(randi(numel(spotLegend1),1,num_spot_eachlabel(1)));
ispotLegend2 = spotLegend2(randi(numel(spotLegend2),1,num_spot_eachlabel(2)));
ispotLegend3 = spotLegend3(randi(numel(spotLegend3),1,num_spot_eachlabel(3)));
ispotLegend4 = spotLegend4(randi(numel(spotLegend4),1,num_spot_eachlabel(4)));
ispotLegend5 = spotLegend5(randi(numel(spotLegend5),1,num_spot_eachlabel(5)));
ispotLegend6 = spotLegend6(randi(numel(spotLegend6),1,num_spot_eachlabel(6)));
ispotLegend7 = spotLegend7(randi(numel(spotLegend7),1,num_spot_eachlabel(7)));
ispotLegend8 = spotLegend8(randi(numel(spotLegend8),1,num_spot_eachlabel(8)));

ispotLegend = [ispotLegend1;ispotLegend2;ispotLegend3;ispotLegend4;ispotLegend5; ...
    ispotLegend6;ispotLegend7;ispotLegend8];

y =800- mod(ispotLegend,800);
x = (ispotLegend-y)/800;
S = scatter(x,y','.');
xlabel('km')
ylabel('km')
xlim([0 700])
ylim([0 800])
axis equal

users_locations = [x,y];

%%  调整并计算每个波束的用户数量
% 计算8*8方格网中每个方格用户的数量，并匹配至波束
R = 50;
Nuser_beam = zeros(64,1);
for idx = 1:8
    xlb = 2*R*690/800  * (idx-1);
    xhb = xlb + 2*R*690/800;
    if idx < 8
        colspots = find(x>=xlb & x<xhb);
    else
        colspots = find(x>=xlb & x<=xhb);
    end
    for idy = 1:8
        ylb = 2*R * (idy-1);
        yhb = ylb + 2*R;
        if idy < 8
            rowspots = find(y>=ylb & y<yhb);
        else
            rowspots = find(y>=ylb & y<=yhb);
        end
        idspot = rowspots(ismember(rowspots,colspots));
        Nuser_beam((idx-1)*8 + idy) = length(idspot);
    end
end
Nuser_beam(1) = Nuser_beam(1)+1;
sum_Nuser_beam = sum(Nuser_beam);
%% 64波束
% 波束中心坐标
K = 64;
R_0 = 49.18;
Distance_beam = 2*R_0;
Center_Beams=zeros(K,2);

Center_Beams(1:sqrt(K),1)=50;
Center_Beams(1:sqrt(K),2)=Distance_beam*(0:sqrt(K)-1);

Center_Beams(sqrt(K)+1:2*sqrt(K),1)=50 + Distance_beam*sqrt(3)/2;
Center_Beams(sqrt(K)+1:2*sqrt(K),2)=Distance_beam*(0:sqrt(K)-1)+Distance_beam/2;

Center_Beams(2*sqrt(K)+1:3*sqrt(K),1)=50 +2*Distance_beam*sqrt(3)/2;
Center_Beams(2*sqrt(K)+1:3*sqrt(K),2)= Distance_beam*(0:sqrt(K)-1);

Center_Beams(3*sqrt(K)+1:4*sqrt(K),1)=50 +3*Distance_beam*sqrt(3)/2;
Center_Beams(3*sqrt(K)+1:4*sqrt(K),2)=Distance_beam*(0:sqrt(K)-1)+Distance_beam/2;

Center_Beams(4*sqrt(K)+1:5*sqrt(K),1)=50 +4*Distance_beam*sqrt(3)/2;
Center_Beams(4*sqrt(K)+1:5*sqrt(K),2)= Distance_beam*(0:sqrt(K)-1);

Center_Beams(5*sqrt(K)+1:6*sqrt(K),1)=50 +5*Distance_beam*sqrt(3)/2;
Center_Beams(5*sqrt(K)+1:6*sqrt(K),2)= Distance_beam*(0:sqrt(K)-1)+Distance_beam/2;

Center_Beams(6*sqrt(K)+1:7*sqrt(K),1)=50 +6*Distance_beam*sqrt(3)/2;
Center_Beams(6*sqrt(K)+1:7*sqrt(K),2)= Distance_beam*(0:sqrt(K)-1);

Center_Beams(7*sqrt(K)+1:8*sqrt(K),1)=50 +7*Distance_beam*sqrt(3)/2;
Center_Beams(7*sqrt(K)+1:8*sqrt(K),2)= Distance_beam*(0:sqrt(K)-1)+Distance_beam/2;
figure('name','uniform')

for ind_beam  =1:K
    % color = ["#0072BD","#A2142F","#EDB120","#7E2F8E","#77AC30","#4DBEEE"];

    pos = [Center_Beams(ind_beam,1)-R_0, ...
        Center_Beams(ind_beam,2)-R_0, ...
        2*R_0, 2*R_0];
    r1 = rectangle('Position',pos,'Curvature',[1 1],'LineStyle',':');
    % r1.EdgeColor = color(ind_beam);
    r1.LineWidth = 1.5;
    hold on
end
scatter(users_locations(:,1),users_locations(:,2),...
    8,'o','MarkerEdgeColor','none','MarkerFaceColor',"#D95319");
hold on

