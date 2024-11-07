clc
clear
close all
%% 取得照片
img = imread('photo_059.jpg'); 
figure()
imshow(img)
title('原圖');

%% 拍攝當下姿態(實際運用應自動讀姿態，不能手動輸入)

roll = 16.65;
pitch = 9.05;
% yaw = 297.45;
yaw = 330;
roll = deg2rad(roll);
pitch = deg2rad(pitch);
yaw = deg2rad(yaw);
% totalAngle = sqrt(pitch^2 + roll^2);

%% 旋轉矩陣
% R = [cos(totalAngle) -sin(totalAngle); sin(totalAngle) cos(totalAngle)];
Rx = [1 0 0; 
      0 cos(roll) -sin(roll); 
      0 sin(roll) cos(roll)];
Ry = [cos(pitch) 0 sin(pitch); 
      0 1 0; 
      -sin(pitch) 0 cos(pitch)];
Rz = [cos(yaw) -sin(yaw) 0;
      sin(yaw) cos(yaw) 0;
      0 0 1;];
R = Rz*Ry*Rx;
%% 計算圖像的四個角點
[h, w, ~] = size(img);
corners = [0 0 1; w 0 1; w h 1; 0 h 1]';

%% 將角點進行旋轉變換
new_corners = R * corners;

%% 計算新的圖像邊界?
min_x = min(new_corners(1, :));
max_x = max(new_corners(1, :));
min_y = min(new_corners(2, :));
max_y = max(new_corners(2, :));
min_z = min(new_corners(3, :));

%% 計算新的圖像尺寸
new_width = ceil(max_x - min_x);   %% ceil 是往正無限大方向取最小整數
new_height = ceil(max_y - min_y);

%% 計算變換矩陣
tform = fitgeotrans(corners(1:2, :)', new_corners(1:2, :)', 'projective');

%% 應用透視變換
outputView = imref2d([new_height new_width], [min_x max_x], [min_y max_y]);
corrected_img = imwarp(img, tform, 'OutputView', outputView);

%% 校正後的圖像
figure()
subplot(2,1,1)
imshow(img)
title('原圖')
subplot(2,1,2)
imshow(corrected_img);
title('校正後的圖像');

%% 將校正後的圖像存起來
imwrite(corrected_img, 'corrected_img.jpg')
