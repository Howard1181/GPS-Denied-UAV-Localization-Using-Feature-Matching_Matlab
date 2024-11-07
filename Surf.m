%%  Surf 演算法
%   用 Surf 做特徵提取 並匹配圖像

% 1.  載入測試圖像
% 2.  轉灰階後分別檢測SURF特徵點
% 3.  分別提取SURF描述子，即特徵向量
% 4.  匹配兩個圖像的特徵向量
% 5.  利用匹配結果計算兩者之間的transform關係
% 6.  透過MSAC將錯誤的特徵點剔除
% 7.  將原本的小圖拼接至大圖中
% 8.  根據小圖在大圖上的位置與變換關係tform，在大圖上框出小圖位置
% 9.  計算小圖的中心點位置(pixel)=>算重心
% 10. 建立經緯度索引標籤，再透過小圖位置(pixel)對應真實地圖的經緯度，而且可能還要寫個內插算沒有經緯度座標的點
% ******************************** 目前進度到第10項 *************************************
% 11. 用磁力計(這並不準，但沒辦法，因為我們沒有GPS)、高度計和角度(飛機姿態角+鏡頭安裝角)推出水平距離，
%     再和前面算出的位置做補償可得飛機實際位置
% PS. 想一下該如何實作連續定位(即時性如何)，需要離線衛星地圖和空拍資料才能測試...

clc; clear; close all;

%% 讀取圖片
rgbImage = imread('G_FullMap_FJU5.png'); % Google 地圖(高度 X 寬度 X 3)，應自動讀取相片，開發階段先設1秒讀取一次
% rgbImage = histeq(rgbImage);
% img1 = imgaussfilt(img1, 2);

% rgbImage = cat(3, img1(:,:,1), img1(:,:,2), img1(:,:,3));


img2 = imread('逢甲空拍/094.jpg'); % 無人機實拍圖
% img2 = imrotate(img2, -180);
% img2 = histeq(img2);
% img2 = imgaussfilt(img2, 2);

% 調整img解析度
sfactor1 = 1;
sfactor2 = 0.1;
rgbImage = imresize(rgbImage, sfactor1);
% rgbImage = rgbImage + 10;    % 提高亮度
img2 = imresize(img2, sfactor2);
% img2 = img2 - 30;    % 降低亮度
% imwrite(img2, 'Resize_pic.jpg')
% figure;
% imshow(img2);
% title('Resized Image');
% img3 = imread('twdtm_asterV2_30m.tif'); % 仕豪給的台灣地圖(測試)

% 顯示結果
figure('Name','讀取彩色圖片')
subplot(2,1,1)
imshow(rgbImage) % 第3維的size = 3，是因為rgb???
title("離線地圖")
subplot(2,1,2)
imshow(img2)
title("無人機實拍畫面")  

%% 圖片姿態校正
roll = -2.69;   % 應自動讀取EXCEL表
pitch = 9.32;
yaw = 310;

% roll = 16.65;
% pitch = 9.05;
% yaw = 297.45;

roll = deg2rad(roll);
pitch = deg2rad(pitch);
yaw = deg2rad(yaw);
% totalAngle = sqrt(pitch^2 + roll^2);

% 旋轉矩陣
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
R = Rx*Ry*Rz;

% 計算圖像的四個角點
[h, w, ~] = size(img2);
corners = [0 0 1; w 0 1; w h 1; 0 h 1]';

% 將角點進行旋轉變換
new_corners = R * corners;

% 計算新的圖像邊界(?
min_x = min(new_corners(1, :));
max_x = max(new_corners(1, :));
min_y = min(new_corners(2, :));
max_y = max(new_corners(2, :));
min_z = min(new_corners(3, :));

% 計算新的圖像尺寸
new_width = ceil(max_x - min_x);   %% ceil 是往正無限大方向取最小整數
new_height = ceil(max_y - min_y);

% 計算變換矩陣
tform = fitgeotrans(corners(1:2, :)', new_corners(1:2, :)', 'projective');

% 應用透視變換
outputView = imref2d([new_height new_width], [min_x max_x], [min_y max_y]);
AttCorrected_img2 = imwarp(img2, tform, 'OutputView', outputView);

% 校正後的圖像
% figure()
% subplot(2,1,1)
% imshow(img2)
% title('原圖')
% subplot(2,1,2)
% imshow(AttCorrected_img2);
% title('姿態校正後的圖像');

%% 圖片切割
% New_img1 = img1; % 不想覆蓋原本的矩陣，所以另外存一個矩陣
% % 把矩陣大小補零，這樣在做切割時才不會有問題
% numrow = size(img1, 1); % 取img1的row數
% numcol = size(img1, 2);
% numrow = ceil(numrow / 100) * 100; % 將numrow往上取到100的倍數
% numcol = ceil(numcol / 100) * 100;
% 
% New_img1(numrow,numcol,3)=0;  % Gb1.JPG 這個要
% figure;
% imshow(New_img1)
% SkipStep = 100; % 每個切割後圖片塊的大小
% M_img1 = 100; % 圖片塊的長
% N_img1 = 100; % 圖片塊的寬
% n = 0; % 圖片塊的編號
% Double_img1 = im2double(New_img1); % 把圖片轉成雙精值
% [H,W,t] = size(Double_img1); % 得到矩陣的大小
% xStepNum = floor((W-N_img1)/SkipStep+1); % 朝負無窮方向取整，寬度方向block移動的次數
% yStepNum = floor((H-M_img1)/SkipStep+1); % 朝負無窮方向取整，長度方向block移動的次數
% for i = 1:xStepNum
%     for j = 1:yStepNum
%         n = n + 1;
%         P_img1 = Double_img1((j-1)*SkipStep+1:(j-1)*SkipStep+M_img1, (i-1)*SkipStep+1:(i-1)*SkipStep+N_img1, :); % 分割圖像
%         Cut_a = strcat('Cut_GFig\',num2str(n),'.jpg'); % 儲存的圖片位置及每幅圖片塊的命名
%         imwrite(double(P_img1),Cut_a);        
%     end
% end
% 
% %% 圖片拼接(ALL)：目前是將全部的地圖載入再進行拼接
% Puzzle_img1 = zeros(numrow,numcol,3); % 補零後的原圖片大小為 ()
% m = 0;
% for i = 1:xStepNum
%     for j = 1:yStepNum
%         m = m + 1;
%         % 拼接圖像
%         Puzzle_a = strcat('Cut_GFig\',num2str(m),'.jpg'); % 儲存的圖片位置及每幅圖片塊的命名
%         Puzzle_b = im2double(imread(Puzzle_a));     
%         Puzzle_img1((j-1)*SkipStep+1:(j-1)*SkipStep+M_img1, (i-1)*SkipStep+1:(i-1)*SkipStep+N_img1, :) = Puzzle_b;
%     end
% end
% w = strcat('Puzzle.jpg'); % 將合成圖片儲存至想要的位置並命名
% figure('Name','拼接完後的圖片[總圖]')
% imshow(double(Puzzle_img1));
% imwrite(double(Puzzle_img1),w);
% 
% %% 圖片拼接(3*3)：根據載入的圖號拼出周為8塊圖
% % % 因為現在只是圖片，所以只能直接給編號當位置，之後應該是根據經緯度去抓我們要的那張圖片，這部分之後再寫
% % promp1 = '請輸入圖片編號(1 ~ ';
% % promp2 = num2str(xyStepNum);
% % promp3 = ' 其中一個):';
% % promp = append(promp1,promp2,promp3);
% % % 選擇要的地圖編號
% % p = input(promp);
% 
% % 選擇要的地圖編號
% xyStepNum = xStepNum*yStepNum;
% p = input(append('請輸入圖片編號(1 ~ ',num2str(xyStepNum),' 其中一個):'));
% 
% % 再將那張圖片及他周遭的8張一併讀取出來，考量可能會抓到邊界或是計算到的圖片編號不在範圍內的情況，若發生這種現象則用零取代
% if (p-yStepNum-1)<1 || (p-yStepNum-1)>xyStepNum
%     MapNum1 =zeros([M_img1,N_img1,3]);
% else
%     MapNum1 = im2double(imread(['Cut_GFig\',num2str(p-yStepNum-1),'.jpg']));
% end
% 
% if (p-yStepNum)<1 || (p-yStepNum)>xyStepNum
%     MapNum2 =zeros([M_img1,N_img1,3]);
% else
%     MapNum2 = im2double(imread(['Cut_GFig\',num2str(p-yStepNum),'.jpg']));
% end
% 
% if (p-yStepNum+1)<1 || (p-yStepNum+1)>xyStepNum
%     MapNum3 =zeros([M_img1,N_img1,3]);
% else
%     MapNum3 = im2double(imread(['Cut_GFig\',num2str(p-yStepNum+1),'.jpg']));
% end
% 
% if (p-1)<1 || (p-1)>xyStepNum
%     MapNum4 =zeros([M_img1,N_img1,3]);
% else
%     MapNum4 = im2double(imread(['Cut_GFig\',num2str(p-1),'.jpg']));
% end
% 
% if (p)<1 || (p)>xyStepNum
%     MapNum5 =zeros([M_img1,N_img1,3]);
% else
%     MapNum5 = im2double(imread(['Cut_GFig\',num2str(p),'.jpg']));
% end
% 
% if (p+1)<1 || (p+1)>xyStepNum
%     MapNum6 =zeros([M_img1,N_img1,3]);
% else
%     MapNum6 = im2double(imread(['Cut_GFig\',num2str(p+1),'.jpg']));
% end
% 
% if (p+yStepNum-1)<1 || (p+yStepNum-1)>xyStepNum
%     MapNum7 =zeros([M_img1,N_img1,3]);
% else    
%     MapNum7 = im2double(imread(['Cut_GFig\',num2str(p+yStepNum-1),'.jpg']));
% end
% 
% if (p+yStepNum)<1 || (p+yStepNum)>xyStepNum
%     MapNum8 =zeros([M_img1,N_img1,3]);
% else    
%     MapNum8 = im2double(imread(['Cut_GFig\',num2str(p+yStepNum),'.jpg']));
% end
% 
% if (p+yStepNum+1)<1 || (p+yStepNum+1)>xyStepNum
%     MapNum9 =zeros([M_img1,N_img1,3]);
% else
%     MapNum9 = im2double(imread(['Cut_GFig\',num2str(p+yStepNum+1),'.jpg']));
% end
% 
% % 避免上下邊界的誤判
% for  count = 1:1:xStepNum % 第一排和最後一排分別有11個數字要判斷  
%     % 若 p 的值位於第一排
%     if p == (count-1)*yStepNum+1
%         % 就把 MapNum1、MapNum4 和 MapNum7 用零矩陣替換掉
%         MapNum1 = zeros([M_img1,N_img1,3]);  
%         MapNum4 = zeros([M_img1,N_img1,3]);   
%         MapNum7 = zeros([M_img1,N_img1,3]);
%     end
%     % 若 p 的值位於最末排
%     if p == (count)*yStepNum
%         % 就把 MapNum3、MapNum6 和 MapNum9 用零矩陣替換掉
%         MapNum3 = zeros([M_img1,N_img1,3]);  
%         MapNum6 = zeros([M_img1,N_img1,3]);   
%         MapNum9 = zeros([M_img1,N_img1,3]);
%     end
% end
% 
% % 再把九張圖片照順序拼回去
% Puzzle_c =[MapNum1 MapNum4 MapNum7; MapNum2 MapNum5 MapNum8; MapNum3 MapNum6 MapNum9]; 
% 
% % 把拼好的圖片顯示出來
% figure('Name','拼接完後的圖片[3*3]');
% imshow(double(Puzzle_c));
% 
% % 並把圖片儲存起來
% PictureName = strcat('Puzzle99.jpg'); % 將合成圖片儲存至想要的位置並命名
% imwrite(double(Puzzle_c),PictureName);

%% 圖片轉灰階
% img1 = imread('Puzzle9.jpg');  % 將9宮格地圖重新取代img1，須根據無人機位置選擇九宮格中心點(無人機目前位置讀取之周遭九宮格圖) -> 待確認
Gray_img1 = rgb2gray(rgbImage);
% Gray_img1 = histeq(Gray_img1);

Gray_img2 = rgb2gray(img2);
% Gray_img2 = histeq(Gray_img2);
% Gray_img2 = imgaussfilt(Gray_img2, 2);

% Gray_img3 = rgb2gray(AttCorrected_img2);
imageSize = size(rgbImage);
[M,N] = size(Gray_img2);
% 顯示結果
figure('Name','轉成灰階圖片')
subplot(2,1,1)
imshow(Gray_img1) 
title("Google地圖")
subplot(2,1,2)
imshow(Gray_img2)
title("無人機實拍地圖")

%% 圖片二值化(SURF只需要灰度圖，所以這部分用不到)
% Bi_mg1 = im2bw(Gray_img1,0.7); 
% Bi_mg2 = im2bw(Gray_img2,0.7);

%% 找圖像的特徵點
tic
points1 = detectSURFFeatures(Gray_img1);
points2 = detectSURFFeatures(Gray_img2);  % 可以試著調 Threshold
% points3 = detectSURFFeatures(Gray_img3); % 經姿態校正的空拍照片

% 顯示結果
figure('Name','對左右圖像分別進行特徵點提取')
subplot(2,1,1)
imshow(Gray_img1); hold on;
%title("The feature points of Google map")
[Lsize, m] = size(points1.Location);
plot(points1(1:2:Lsize,:),'showOrientation',true);
subplot(2,1,2)
imshow(Gray_img2); hold on;
%title("The feature points Microsoft map")
[Rsize, n] = size(points2.Location);
plot(points2(1:2:Rsize,:),'showOrientation',true); 

% figure('Name', '對進行姿態校正後的圖進行特徵點提取')
% imshow(Gray_img3); hold on;
% [Asize, l] = size(points3.Location);
% plot(points3(1:2:Asize, :), 'showOrientation', true);

%% 影像校正：因為相機會有個角度(安裝角+飛機姿態)，所以影像會變形(類似梯形)，因此需要把影像拉回成正投影的狀態
% 經過測試發現SURF演算法有類似幾何校正的功能，所以目前評估先不用寫這部分
% Theta = 0:180; % 0:179;
% R = radon(Gray_img2,Theta);
% % 求出圖像中心點至邊界的距離
% L = round(sqrt((M/2)^2+(N/2)^2));
% [C,angle] = max(R(L,:)); % 角度應該我們自己給會比較好
% % angle 為圖像傾斜角度
% angle = angle - 1;
% % 將圖片做校正
% A = imrotate(img2,angle,'nearest');
% 
% % 顯示校正後的圖片
% figure;imshow(A);

%% 計算描述向量
[f1, vpts1] = extractFeatures(Gray_img1, points1);
[f2, vpts2] = extractFeatures(Gray_img2, points2);
% [f3, vpts3] = extractFeatures(Gray_img3, points3);

%% 進行匹配；提取特徵點位置，SURF 特徵向量有先做正規化處理
% 經姿態校正前的空拍圖
% indexPairs1 = matchFeatures(f1, f2) ;
indexPairs1 = matchFeatures(f1, f2, 'Prenormalized', true) ; 
% indexPairs1 = matchFeatures(f1, f2, 'Method', 'Approximate') ;
% indexPairs1 = matchFeatures(f1, f2, 'Unique', true) ;
matched_pts1 = vpts1(indexPairs1(:, 1));
matched_pts2 = vpts2(indexPairs1(:, 2));

% 經姿態校正後的空拍圖
% % indexPairs2 = matchFeatures(f1, f2,'MatchThreshold',2) ;
% indexPairs2 = matchFeatures(f1, f3, 'Prenormalized', true) ; 
% % indexPairs2 = matchFeatures(f1, f2, 'Method', 'Approximate') ;
% matched_pts3 = vpts1(indexPairs2(:, 1));
% matched_pts4 = vpts3(indexPairs2(:, 2));

%% 用MSAC演算法去除錯誤的匹配點
% 通過特徵點匹配還得到了第二幅圖的變換矩陣tform，第二幅圖要經過變換矩陣變成和第一幅圖的座標一致
[tform, inlierimg2Points, inlierimg1Points] = estimateGeometricTransform(matched_pts2, matched_pts1, 'similarity'); %射影變換，tform映射點對1內點到點對2內點
% 該函數使用隨機樣本一致性（RANSAC，Random Sample Consensus）演算法的變體MSAC演算法實現，去除誤匹配點
% 返回的幾何映射矩陣映射第一參數內點到第二參數內點
toc
% 顯示對匹配結果，可以看到原先還有一些異常值，但經過MSAC後有明顯改善
  figure('Name','對左右圖像分別進行特徵點提取')
  subplot(2,1,1)
  showMatchedFeatures(Gray_img1,Gray_img2,matched_pts1,matched_pts2,'montage');
  legend('matched points 1','matched points 2'); 
  title("特徵值匹配結果")
  subplot(2,1,2)
  showMatchedFeatures(Gray_img1,Gray_img2,inlierimg1Points,inlierimg2Points,'montage');
  legend('matched points 1','matched points 2');
  title("去除錯誤的匹配點後的特徵值匹配結果")
  
% 同樣的事將地圖與姿態校正過後圖再做一次
%   [tform_attC, inlierimg2Points, inlierimg1Points] = estimateGeometricTransform(matched_pts4, matched_pts3, 'similarity', 'MaxDistance', 4); %射影變換，tform映射點對1內點到點對2內點
%   figure('Name', '經姿態校正後對左右圖特徵點提取')
%   subplot(2,1,1)
%   showMatchedFeatures(Gray_img1,Gray_img3,matched_pts3,matched_pts4,'montage');
%   title("特徵匹配結果")
%   subplot(2,1,2)
%   showMatchedFeatures(Gray_img1,Gray_img3,inlierimg1Points,inlierimg2Points,'montage');
%   title("去除錯誤點後的特徵匹配結果")
    
%% 進行圖像合併，tform是變換矩陣，以第一幅圖像為基準座標，第二幅圖要進行變換與其對應
% Rfixed 為第一幅圖的世界二維座標
  Rfixed = imref2d(size(rgbImage));
  [registered2, Rregistered] = imwarp(img2, tform);
  boxPolygon = [1, 1;... % 左上
    size(img2, 2), 1; ... % 右上
    size(img2, 2), size(img2, 1); ... % 右下
    1, size(img2, 1); ... % 左下
    1, 1]; % 重複左下，才可得到一個閉區間的多邊形
% 將多邊形變換到目標圖片上，變換的結果表示了物體的位置
  newBoxPolygon = transformPointsForward(tform, boxPolygon);
% 顯示被檢測到的物體
  figure()
  imshowpair(rgbImage,Rfixed,registered2,Rregistered,'blend');
%%% 計算拼接圖片的中心位置(pixel)：之後根據這個去找對應的索引值
% 要算原本的框還是背景的框？但背景的框我要想一下怎麼算
  hold on;
% [xlim, ylim] = outputLimits(tform, [1 imageSize(1)], [1 imageSize(2)]);
% 第一種算重心的方式：算術平均
  Xmdn = (newBoxPolygon(1, 1)+newBoxPolygon(2, 1)+newBoxPolygon(3, 1)+newBoxPolygon(4, 1))/4;
  Ymdn = (newBoxPolygon(1, 2)+newBoxPolygon(2, 2)+newBoxPolygon(3, 2)+newBoxPolygon(4, 2))/4;
  line(newBoxPolygon(:, 1), newBoxPolygon(:, 2), 'Color', 'y');
  hold on;
  InterestPoint1 = scatter(Xmdn, Ymdn, '*b'); % 和matlab算的重心有些微誤差
% 第二種算重心的方式：用 matlab 函數
  x1 = newBoxPolygon(1:4, 1)';
  y1 = newBoxPolygon(1:4, 2)';
  polyin = polyshape(x1,y1); 
  [CenterX,CenterY] = centroid(polyin);% matlab自帶計算重心的函數
  InterestPoint2 = scatter(CenterX, CenterY, '*g'); 
  title('圖像差異');

%% 應用姿態校正過後的圖再做一次
%   [registered_attC, Rregistered_attC] = imwarp(AttCorrected_img2, tform_attC);
%   boxPolygon_attC = [1, 1;... % 左上
%     size(AttCorrected_img2, 2), 1; ... % 右上
%     size(AttCorrected_img2, 2), size(AttCorrected_img2, 1); ... % 右下
%     1, size(AttCorrected_img2, 1); ... % 左下
%     1, 1]; % 重複左下，才可得到一個閉區間的多邊形
%   newBoxPolygon_attC = transformPointsForward(tform_attC, boxPolygon_attC);
%   figure()
%   imshowpair(img1,Rfixed,registered_attC,Rregistered_attC,'blend');
% % 要算原本的框還是背景的框？但背景的框我要想一下怎麼算
%   hold on;
% % [xlim, ylim] = outputLimits(tform, [1 imageSize(1)], [1 imageSize(2)]);
% % 第一種算重心的方式：算術平均
%   Xmdn = (newBoxPolygon_attC(1, 1)+newBoxPolygon_attC(2, 1)+newBoxPolygon_attC(3, 1)+newBoxPolygon_attC(4, 1))/4;
%   Ymdn = (newBoxPolygon_attC(1, 2)+newBoxPolygon_attC(2, 2)+newBoxPolygon_attC(3, 2)+newBoxPolygon_attC(4, 2))/4;
%   line(newBoxPolygon_attC(:, 1), newBoxPolygon_attC(:, 2), 'Color', 'y');
%   hold on;
%   InterestPoint1 = scatter(Xmdn, Ymdn, '*b'); % 和matlab算的重心有些微誤差
% % 第二種算重心的方式：用 matlab 函數
%   x1 = newBoxPolygon_attC(1:4, 1)';
%   y1 = newBoxPolygon_attC(1:4, 2)';
%   polyin = polyshape(x1,y1); 
%   [CenterX,CenterY] = centroid(polyin);% matlab自帶計算重心的函數
%   InterestPoint2 = scatter(CenterX, CenterY, '*g'); 
%   title('圖像差異');
 
%%% 計算拼接圖片的中心位置(pixel)：之後根據這個去找對應的索引值
% 要算原本的框還是背景的框？但背景的框我要想一下怎麼算
%   hold on;
% % [xlim, ylim] = outputLimits(tform, [1 imageSize(1)], [1 imageSize(2)]);
% % 第一種算重心的方式：算術平均
%   Xmdn = (newBoxPolygon_attC(1, 1)+newBoxPolygon_attC(2, 1)+newBoxPolygon_attC(3, 1)+newBoxPolygon_attC(4, 1))/4;
%   Ymdn = (newBoxPolygon_attC(1, 2)+newBoxPolygon_attC(2, 2)+newBoxPolygon_attC(3, 2)+newBoxPolygon_attC(4, 2))/4;
%   line(newBoxPolygon_attC(:, 1), newBoxPolygon_attC(:, 2), 'Color', 'y');
%   hold on;
%   InterestPoint1 = scatter(Xmdn, Ymdn, '*b'); % 和matlab算的重心有些微誤差
% % 第二種算重心的方式：用 matlab 函數
%   x1 = newBoxPolygon_attC(1:4, 1)';
%   y1 = newBoxPolygon_attC(1:4, 2)';
%   polyin = polyshape(x1,y1); 
%   [CenterX,CenterY] = centroid(polyin);% matlab自帶計算重心的函數
%   InterestPoint2 = scatter(CenterX, CenterY, '*g'); 
%   title('圖像差異');
  
 %% 以照片中心點座標計算無人機實際經緯度
 % 離線地圖四個角的經緯度
%   Pos_Map_UL = [24.265632, 120.815443];
%   Pos_Map_LL = [24.263959, 120.815445];
%   Pos_Map_UR = [24.265644, 120.817518];
%   Pos_Map_LR = [24.264007, 120.817536];
%  % 計算 X/Y 軸每個像素是多少經緯度(經緯度先x10e6)
%   Lon_per_pixel = (Pos_Map_UR(2) * 10^6 - Pos_Map_UL(2) * 10^6) / size(img1, 2);
%   Lat_per_pixel = (Pos_Map_UL(1) * 10^6 - Pos_Map_LL(1) * 10^6) / size(img1, 1);
%  % 計算實拍照片經緯度
%   Estimated_Lat = Pos_Map_UL(1) * 10^6 - Ymdn*Lat_per_pixel;
%   Estimated_Lon = Pos_Map_UL(2) * 10^6 + Xmdn*Lon_per_pixel;
%   Estimated_Pos = [Estimated_Lat Estimated_Lon] / 10^6;
%  % 將估算經緯度與實際經緯度畫出來(有可能不準，畢竟地圖四個角是用手點的)
%   hold on
%   plot(Xmdn, Ymdn, 'ro', 'MarkerSize', 12)
%   
%   X_shot_point = (120.8169532 - Pos_Map_UL(2))*10^6 / Lon_per_pixel;
%   Y_shot_point = (Pos_Map_UL(1) - 24.2643635)*10^6 / Lat_per_pixel;
%   hold on
%   plot(X_shot_point, Y_shot_point, 'bo', 'MarkerSize', 12)
%   hold off
  
 % 實際 [24.2643635000000,120.816953200000]
 % 估算 [24.264338,120.81686]
 % 誤差約為10公尺
 
%% 以圖片中心點座標計算實際飛機位置(這裡資訊不完整，還無法模擬)
% % 將像素位置對應到經緯度(我現在還沒有index，所以先填原本的值)
% PictureLon = CenterX;
% PictureLat = CenterY;
% % 鏡頭安裝角度是20度，實際未知
% CameraAngle = 20;
% % 飛機俯仰角
% Xsens_Pitch = Data(:,12); % 俯仰角
% %  Madgwick_Pitch = Data(:,15); % 俯仰角
% CPAngle = 90-(CameraAngle + Xsens_Pitch);
% % 氣壓高度，m(NED)，完了沒有LLA坐標系的高度!!!
% Air_Altitude = Data(:,19); 
% LLA_Altitude = 無法取得lla的高度壓;% 高度要轉LLA坐標系阿不然會錯@@!!!
% 
% % 方法一：近似法
% 計算斜距：這裡用近似法的原因是因為不需要座標轉換，這樣可以降低計算時間，又加上距離很短，所以可以視為圓形計算
%  ARC = 6371.393*1000; % 平均半徑(m)
%  TrueLon = PictureLon + LLA_Altitude/(ARC*cos(PictureLon)*2*pi/360);
%  TrueLat = PictureLat +(LLA_Altitude/cos(CPAngle))/ (ARC *2*pi/360);
% 
% % 方法二：座標轉換
% % LLA 座標系 轉 NED 座標系 [xNorth,yEast,zDown] = geodetic2ned(lat,lon,h,lat0,lon0,h0,spheroid) 
% % 起飛位置：
%   LonS_GPS = Data(1,8); % GPS LLA 座標下的 經度 位置
%   LatS_GPS = Data(1,9); % GPS LLA 座標下的 緯度 位置
%   HS_GPS = Data(1,10); % GPS LLA 座標下的 高度 位置
% % 磁力計
%   MagX = Data(:,8); % X 軸
%   MagY = Data(:,9); % Y 軸 
% % LLA 座標系 轉 NED 座標系 [xNorth,yEast,zDown] = geodetic2ned(lat,lon,h,lat0,lon0,h0,spheroid)，現在還沒有資料 
%  [GPS_N_S,GPS_E_S,GPS_D_S] = geodetic2ned(PictureLon,PictureLat,沒高度好鬱悶,LonS_GPS,LatS_GPS,HS_GPS,wgs84Ellipsoid);
% %  True_N_S = (MagX/MagY) *True_E_S;
% %  (True_N_S - GPS_N_S)^2 + (True_E_S - GPS_E_S)^2 = (Air_Altitude/cos(CPAngle))^2;
% % 一元二次方程式(計算距離)
%  syms True_N_S True_E_S;
%  f = 5*True_N_S+(-2*True_N_S*GPS_N_S-4*True_N_S*GPS_E_S)-(Air_Altitude/cos(CPAngle))^2+GPS_N_S^2+GPS_E_S^2;
%  % 求解方程式
%  True_N_S =solve(f==0);
%  True_E_S = 2*True_N_S;
%  % NED 座標系 轉 LLA 座標系
%  [True_LAT_S, True_LON_S, True_H] = ned2geodetic(True_N_S,True_E_S,又沒高度,LonS_GPS,LatS_GPS,HS_GPS,wgs84Ellipsoid); 
 

%%  *************************************************** 以下是自己的小研究，專案用不到 ***************************************************
%% 確定兩幅圖實際大小及位置 
% tic
% %輸出座標範圍
% [xlim, ylim] = outputLimits(tform, [1 imageSize(1)], [1 imageSize(2)]);
% 
% % 找到輸出空間限制的最大最小值
% xMin = min([1; xlim(:)]);
% xMax = max([imageSize(2); xlim(:)]); 
% yMin = min([1; ylim(:)]);
% yMax = max([imageSize(1); ylim(:)]);
%  
% % 全景圖的寬高
% width  = round(xMax - xMin);
% height = round(yMax - yMin); 
%  
% %創建2D空間參考物件定義全景圖尺寸
% xLimits = [xMin xMax];
% yLimits = [yMin yMax];
% panoramaView = imref2d([height width ], xLimits, yLimits);
%  
% % 變換圖片到全景圖
% unwarpedImage = imwarp(img1,projective2d, 'OutputView', panoramaView);
% warpedImage = imwarp(img2, tform, 'OutputView', panoramaView);
% 
% %% 漸入漸出融合：消除兩張圖之間的邊緣，讓其完美融入
% % 所謂漸入漸出就是將兩幅圖重合的區域按照距離兩幅圖的距離按照一定的權重重新分配重合部分圖畫的三原色權重
% % 比如最中間的就是0.5 0.5的比例。下一步就是找到他們重疊區域，也就是相同掩模區
% % MAKE MASKS FOR BOTH IMAGES
% % warpedImage(isnan(warpedImage))=0;
%  
% % newImage = zeros(size(warpedImage));
% % newImage(1:size(I1,1), 1: size(I1,2)+100,:) = I1;
% newImage=unwarpedImage;
% newImage=double(newImage);
%  
% balck1=(warpedImage(:,:,1)==0 & warpedImage(:,:,2)==0 & warpedImage(:,:,3)==0);
% balck2=(newImage(:,:,1)==0 & newImage(:,:,2)==0 & newImage(:,:,3)==0);
% black=and(balck1,balck2);
% black=~black;
%  
% maskA = (warpedImage(:,:,1)>0 |warpedImage(:,:,2)>0 | warpedImage(:,:,3)>0);%變換圖像掩膜
% mask1 = (newImage(:,:,1)>0 | newImage(:,:,2)>0 | newImage(:,:,3)>0);%非變換圖像掩膜
% mask1 = and(maskA, mask1);%重疊區掩膜
% % figure,imshow(mask1)
% [row,col] = find(mask1==1);
% left = min(col);
% right = max(col);%獲得重疊區左右範圍
% up=min(row);
% down=max(row);
% mask = ones(size(mask1));
% % figure()
% % imshow(mask)
% %mask(:,left:right) = repmat(linspace(0,1,right-left+1),size(mask,1),1);%複製平鋪矩陣
% mask(up:down,:) = repmat(linspace(1,0,down-up+1)',1,size(mask,2));%複製平鋪矩陣
% % BLEND EACH CHANNEL
% warpedImage=double(warpedImage);
% % figure()
% % warpedImage=uint8(warpedImage);
% % imshow(warpedImage)
% % figure()
% % imshow(mask)
% warpedImage(:,:,1) = warpedImage(:,:,1).*mask;
% warpedImage(:,:,2) = warpedImage(:,:,2).*mask;
% warpedImage(:,:,3) = warpedImage(:,:,3).*mask;
%  
% % REVERSE THE ALPHA VALUE
% %mask(:,left:right) = repmat(linspace(1,0,right-left+1),size(mask,1),1);
% mask(up:down,:) = repmat(linspace(0,1,down-up+1)',1,size(mask,2));%複製平鋪矩陣
% newImage(:,:,1) = newImage(:,:,1).*mask;
% newImage(:,:,2) = newImage(:,:,2).*mask;
% newImage(:,:,3) = newImage(:,:,3).*mask;
%  
% newImage(:,:,1) = warpedImage(:,:,1) + newImage(:,:,1);
% newImage(:,:,2) = warpedImage(:,:,2) + newImage(:,:,2);
% newImage(:,:,3) = warpedImage(:,:,3) + newImage(:,:,3);
%  
% % newImage(:,:,1) = newImage(:,:,1).*black;
% % newImage(:,:,2) = newImage(:,:,2).*black;
% % newImage(:,:,3) = newImage(:,:,3).*black;
% newImage=uint8(newImage);
% figure()
% imshow(newImage);
% title('漸入漸出融合');
% toc


