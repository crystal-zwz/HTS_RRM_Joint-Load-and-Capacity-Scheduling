% close all
%% 获取川藏交界处图片，面积，经纬度
% I0= imread("C:\Users\dsczb\Desktop\mainlandlandscan\landscan-hd-china-v1-colorized.tif");
[A,R] = readgeoraster("C:\Users\dsczb\Desktop\landscan-global-2022-colorized.tif");
% A= A(:,:,1:3);
% I2 = imresize(A,0.001);
% I2 = A(1:100:end,1:100:end,1:3);
% figure(1)
% imshow(I2);

figure("Name","based on 1 point 1km")
%6.2/7.9*43200,1.2/3.9*21600
%成都平原(6601:7400,33701,34590)
I0 = A(6601:7400,33701:34390,1:3); 
imshow(I0);

% 33901:34590, 43200 Latitude;6640:7440 21600 Longitude
% 估计经纬度
CZ_Longitude(1) = 180-360*(43200-33701)/43200;
CZ_Longitude(2) = 180-360*(43200-34390)/43200;
CZ_Latitude(2) = 90-180*6601/21600;
CZ_Latitude(1) = 90-180*7400/21600;

% 估计面积
a = (CZ_Longitude(2)-CZ_Longitude(1))*111.1;
b = (CZ_Latitude(2)-CZ_Latitude(1))*111.1*cos(deg2rad(0.5*(CZ_Latitude(2)+CZ_Latitude(1))));
S = a*b;
S_compute = 690*800;
% diffS = S/S_compute;
R_a = 800/a;
R_b = 690/b;

% 修正点数
numspota = round(800*R_a); % 1003
numspotb = round(690*R_b); % 756
I1 = A(6601:7603,33701:34456,1:3); 
figure("Name","amendment")
imshow(I1);

% 再次估计经纬度
CZ_Longitude1(1) = 180-360*(43200-33701)/43200;
CZ_Longitude1(2) = 180-360*(43200-34456)/43200;
CZ_Latitude1(2) = 90-180*6601/21600;
CZ_Latitude1(1) = 90-180*7603/21600;
% 再次估计面积
a1 = (CZ_Longitude1(2)-CZ_Longitude1(1))*111.1;
b1 = (CZ_Latitude1(2)-CZ_Latitude1(1))*111.1*cos(deg2rad(0.5*(CZ_Latitude(2)+CZ_Latitude(1))));
S1 = a1*b1;