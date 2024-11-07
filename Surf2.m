%% 以不同特徵提取方法嘗試
clear; clc; close all;
tic

%% 讀取圖片
  img1 = imread('離線地圖/8K_FJU.png'); % 這裡在實際嵌入式運用中應讀取的是裁切完的照片
  img2FilePath = '逢甲空拍/123.jpg';
  img2 = imread(img2FilePath);
  % 使用 fileparts 提取文件夾路徑、文件名與副檔名
  [fileFolder, imageFileName, ext] = fileparts(img2FilePath);
  imageFileName = [imageFileName , ext];
% 提取TIF照片經緯度資訊(這裡在實際運用中應事前處理，先儲存在嵌入式系統內存中)
  GEOinfo = geotiffinfo('離線地圖/999.tif');
  TL_coor = [GEOinfo.BoundingBox(1), GEOinfo.BoundingBox(4)]; % 地圖四個角經緯度，經度先
  TR_coor = [GEOinfo.BoundingBox(2), GEOinfo.BoundingBox(4)];
  BL_coor = [GEOinfo.BoundingBox(1), GEOinfo.BoundingBox(3)];
  BR_coor = [GEOinfo.BoundingBox(2), GEOinfo.BoundingBox(3)];
  
% 定義地圖資訊
  Map_Width = 3677; 
  Map_Height = 4142;
  CropSize = 350; 
  numrow = ceil(Map_Width / CropSize) * CropSize; % 將numrow往上取到100的倍數
  numcol = ceil(Map_Height / CropSize) * CropSize;

% 圖像調整
%   h = fspecial('unsharp'); 
%   img1 = imgaussfilt(img1, 2);  
%   img1 = histeq(img1);
%   img1 = img1 + 30;  
%   img1 = cat(3, img1(:,:,1), img1(:,:,2), img1(:,:,3));

% 圖像調整
  % img2 = imgaussfilt(img2);
%   img2 = img2 - 30;
%   img2 = histeq(img2);
%   img2 = imhistmatch(img2, img1);

%%  
% 地圖圖片裁切(若有換底圖的話)
%  MapCrop()
  
% 提取空拍圖資訊(需自動尋找Excel表)
  excelFile = 'Log.csv';
  data = readtable(excelFile); % 使用 readtable 讀取 Excel
  fileNames = data{:, 1}; % 提取第一列的文件名
  idx = find(strcmp(fileNames, imageFileName)); % 找到文件名的位置
  if ~isempty(idx)
    AerialP_Lon = data{idx, 10};  % 用()提取table,用{}提取實際值
    AerialP_Lat = data{idx, 9};
    Flight_Yaw = data{idx, 12};
  else
    disp('未找到對應文件名');
  end
  
% 調整img解析度
  % sfactor1 = 1;
  sfactor2 = 0.2;  % 想想如何根據兩張照片性質自動調整scale?
  % img1 = imresize(img1, sfactor1);
  img2 = imresize(img2, sfactor2);
% 調整圖片方向
  img2 = imrotate(img2, -Flight_Yaw); % 根據excel yaw值旋轉圖片，imrotate(img, angle)->逆時針轉angle度

% 顯示結果
  % figure('Name','讀取彩色圖片')
  % subplot(2,1,1)
  % imshow(img1)
  % title("離線地圖")
  % subplot(2,1,2)
  % imshow(img2)
  % title("無人機實拍畫面") 

%% 根據無人機失去GPS前的選取週遭離線地圖
% 根據座標進行週遭地圖拼接
  % 計算一個像素是多少經緯度
  Lon_per_pix = (round(TR_coor(1),8) - round(TL_coor(1),8)) / Map_Width;  % 這些應要都是已知值
  Lat_per_pix = (round(TL_coor(2),8) - round(BL_coor(2),8)) / Map_Height;
  
  % 計算丟失GPS前經緯度在2離線地圖哪個Pixel，第一次迭代之後需更新成演算法計算出來的照片的位置
  AerialP_pixX = abs((round(AerialP_Lon,8) - round(BL_coor(1),8)) / Lon_per_pix);
  AerialP_pixY = abs((round(AerialP_Lat,8) - round(BL_coor(2),8)) / Lat_per_pix);
  
  % 計算pixel對應的圖片
  px = abs(AerialP_pixX - 1);   
  py = abs(Map_Height - AerialP_pixY - 1);
  block_X = ceil(px / CropSize);
  block_Y = ceil(py / CropSize);
%   block_X = floor(px / CropSize);
%   block_Y = floor(py / CropSize);
  imshow(img1)
  picture_num = (numcol/CropSize)*(block_X-1) + block_Y; % 座標點會落在picture_num這張小圖片上，圖片編號由上而下、由左至右編號
  hold on
  plot(AerialP_pixX-1, Map_Height-AerialP_pixY-1, 'ro', 'MarkerSize', 15)
  hold off
  
  % 從資料夾中讀取picture_num對應的小圖
  folderPath = 'Cropped_Map';
  imgName = sprintf('%d.png', picture_num);
  imagePath = fullfile(folderPath, imgName);
  if exist(imagePath, 'file')
    % 讀取圖片
    selectedImg = imread(imagePath);
  
    % 顯示圖片
%     figure()
%     imshow(selectedImg);
  else
    disp('指定的圖片檔案不存在');
  end
 
%% 拼接小圖週遭地圖
xStepNum = floor((numrow-CropSize)/CropSize+1); % 朝負無窮方向取整，寬度方向block移動的次數
yStepNum = floor((numcol-CropSize)/CropSize+1); % 朝負無窮方向取整，長度方向block移動的次數
xyStepNum = xStepNum*yStepNum;
% 將小圖片及他周遭的8張一併讀取出來，考量可能會抓到邊界或是計算到的圖片編號不在範圍內的情況，若發生這種現象則用零取代
if (picture_num-yStepNum-1)<1 || (picture_num-yStepNum-1)>xyStepNum
    MapNum1 =zeros([CropSize,CropSize,3]);
else
    MapNum1 = im2double(imread(['Cropped_Map\',num2str(picture_num-yStepNum-1),'.png']));
end

if (picture_num-yStepNum)<1 || (picture_num-yStepNum)>xyStepNum
    MapNum2 =zeros([CropSize,CropSize,3]);
else
    MapNum2 = im2double(imread(['Cropped_Map\',num2str(picture_num-yStepNum),'.png']));
end

if (picture_num-yStepNum+1)<1 || (picture_num-yStepNum+1)>xyStepNum
    MapNum3 =zeros([CropSize,CropSize,3]);
else
    MapNum3 = im2double(imread(['Cropped_Map\',num2str(picture_num-yStepNum+1),'.png']));
end

if (picture_num-1)<1 || (picture_num-1)>xyStepNum
    MapNum4 =zeros([CropSize,CropSize,3]);
else
    MapNum4 = im2double(imread(['Cropped_Map\',num2str(picture_num-1),'.png']));
end

if (picture_num)<1 || (picture_num)>xyStepNum
    MapNum5 =zeros([CropSize,CropSize,3]);
else
    MapNum5 = im2double(imread(['Cropped_Map\',num2str(picture_num),'.png']));
end

if (picture_num+1)<1 || (picture_num+1)>xyStepNum
    MapNum6 =zeros([CropSize,CropSize,3]);
else
    MapNum6 = im2double(imread(['Cropped_Map\',num2str(picture_num+1),'.png']));
end

if (picture_num+yStepNum-1)<1 || (picture_num+yStepNum-1)>xyStepNum
    MapNum7 =zeros([CropSize,CropSize,3]);
else    
    MapNum7 = im2double(imread(['Cropped_Map\',num2str(picture_num+yStepNum-1),'.png']));
end

if (picture_num+yStepNum)<1 || (picture_num+yStepNum)>xyStepNum
    MapNum8 =zeros([CropSize,CropSize,3]);
else    
    MapNum8 = im2double(imread(['Cropped_Map\',num2str(picture_num+yStepNum),'.png']));
end

if (picture_num+yStepNum+1)<1 || (picture_num+yStepNum+1)>xyStepNum
    MapNum9 =zeros([CropSize,CropSize,3]);
else
    MapNum9 = im2double(imread(['Cropped_Map\',num2str(picture_num+yStepNum+1),'.png']));
end

% 避免上下邊界的誤判
for  count = 1:1:xStepNum % 第一排和最後一排分別有11個數字要判斷  
    % 若 picture_num 的值位於第一排
    if picture_num == (count-1)*yStepNum+1
        % 就把 MapNum1、MapNum4 和 MapNum7 用零矩陣替換掉
        MapNum1 = zeros([CropSize,CropSize,3]);  
        MapNum4 = zeros([CropSize,CropSize,3]);   
        MapNum7 = zeros([CropSize,CropSize,3]);
    end
    % 若 picture_num 的值位於最末排
    if picture_num == (count)*yStepNum
        % 就把 MapNum3、MapNum6 和 MapNum9 用零矩陣替換掉
        MapNum3 = zeros([CropSize,CropSize,3]);  
        MapNum6 = zeros([CropSize,CropSize,3]);   
        MapNum9 = zeros([CropSize,CropSize,3]);
    end
end

% 再把九張圖片照順序拼回去
Puzzle_c =[MapNum1 MapNum4 MapNum7; MapNum2 MapNum5 MapNum8; MapNum3 MapNum6 MapNum9]; 

% 把拼好的圖片顯示出來
figure('Name','拼接完後的圖片[3*3]');
imshow(double(Puzzle_c));

% 並把圖片儲存起來
PictureName = strcat('小地圖拼接/Puzzle99.png'); % 將合成圖片儲存至想要的位置並命名
imwrite(double(Puzzle_c),PictureName);

%% 圖片轉灰階
  % 調整圖片
  Gray_img1 = rgb2gray(Puzzle_c);

%   Gray_img1 = histeq(Gray_img1);
%   Gray_img1 = imadjust(Gray_img1, [0 1], [0.1 0.9]);
%   Gray_img1 = adapthisteq(Gray_img1);
%   Gray_img1 = Gray_img1 - 0.12;
  Gray_img2 = rgb2gray(img2);
  Gray_img2 = histeq(Gray_img2);
%   Gray_img2 = imhistmatch(Gray_img2, Gray_img1); % 光度均一化
  
  % Gray_img2 = adapthisteq(Gray_img2);

%% 判斷是否需調整亮度、對比
% 計算兩張圖片的亮度（平均值）和對比度（標準差）
Gray_img2 = double(Gray_img2) / 255;  % 將uint8圖片轉換成雙精度
mean1 = mean(Gray_img1(:));
mean2 = mean(Gray_img2(:));
std1 = std(Gray_img1(:));
std2 = std(Gray_img2(:));

% 顯示亮度與對比度的差異
disp(['Image 1 亮度: ', num2str(mean1), ', 對比度: ', num2str(std1)]);
disp(['Image 2 亮度: ', num2str(mean2), ', 對比度: ', num2str(std2)]);

%% 
% % 圖像對比度篩選(看是否需要做均衡化)
%   contrastValue_img1 = std(double(Gray_img1(:)));
%   contrastValue_img2 = std(double(Gray_img2(:)));
%   threshold = 20;
% %   if contrastValue_img1 < threshold
% %       Gray_img1 = histeq(Gray_img1);      % 調整對比
% %   end
% %   
% %   if contrastValue_img2 < threshold
% %       Gray_img2 = histeq(Gray_img2);
% %   end
%   
%   
%   % 计算灰度图的均值和方差
%   mean1 = mean(Gray_img1(:));  % (:)是將矩陣展開變一維向量
%   mean2 = mean(Gray_img2(:));
%   std1 = std(double(Gray_img1(:)));
% %   std2 = std(double(Gray_img2(:)));
% 
%   % 设定均值和方差的阈值
%   meanThreshold = 30;  % 均值差异的阈值
%   stdThreshold = 30;   % 方差差异的阈值
% 
%   % 判断是否进行直方图匹配
% %   if abs(mean1 - mean2) > meanThreshold || abs(std1 - std2) > stdThreshold
% %     Gray_img2 = imhistmatch(Gray_img2, Gray_img1);
% %   end
  
  
  % 顯示結果
  figure('Name','轉成灰階圖片')
  subplot(2,1,1)
  imshow(Gray_img1) 
  title("地圖")
  subplot(2,1,2)
  imshow(Gray_img2)
  title("無人機實拍地圖")

%% 特徵檢測
% OEB
% points_img1 = detectORBFeatures(Gray_img1);
% points_img2 = detectORBFeatures(Gray_img2);

% KAZE
  points_img1 = detectKAZEFeatures(Gray_img1);
  points_img2 = detectKAZEFeatures(Gray_img2);

% SURF
%   points_img1 = detectSURFFeatures(Gray_img1);
%   points_img2 = detectSURFFeatures(Gray_img2);

  points_img1 = points_img1.selectStrongest(15000);
  points_img2 = points_img2.selectStrongest(15000);

% 顯示結果
% figure('Name','對圖像分別進行特徵點提取')
% subplot(2,1,1)
% imshow(Gray_img1); hold on;
% %title("The feature points of Google map")
% [MapSize, m] = size(points_img1.Location);
% plot(points_img1(1:2:MapSize,:),'showOrientation',true);
% subplot(2,1,2)
% imshow(Gray_img2); hold on;
% %title("The feature points Microsoft map")
% [ImgSize, n] = size(points_img2.Location);
% plot(points_img2(1:2:ImgSize,:),'showOrientation',true); 

%% 計算描述符
  [f1, vpts1] = extractFeatures(Gray_img1, points_img1);
  [f2, vpts2] = extractFeatures(Gray_img2, points_img2);

%% 進行匹配
  indexPairs = matchFeatures(f1, f2, 'MatchThreshold', 20, 'MaxRatio', 0.7) ; 
  matched_pts1 = vpts1(indexPairs(:, 1));
  matched_pts2 = vpts2(indexPairs(:, 2));

%% 去除錯誤點
  [tform, inlierimg2Points, inlierimg1Points] = estimateGeometricTransform(matched_pts2, matched_pts1, 'similarity', 'MaxNumTrials', 5000, 'MaxDistance', 10);
  % MaxDistance可試著條大一些
  toc
  
% 可視化匹配結果
  figure('Name','對左右圖像分別進行特徵點提取')
  subplot(2,1,1)
  showMatchedFeatures(Gray_img1,Gray_img2,matched_pts1,matched_pts2,'montage');
  legend('matched points 1','matched points 2'); 
  title("特徵值匹配結果")
  subplot(2,1,2)
  showMatchedFeatures(Gray_img1,Gray_img2,inlierimg1Points,inlierimg2Points,'montage');
  legend('matched points 1','matched points 2');
  title("去除錯誤的匹配點後的特徵值匹配結果")

%% 疊圖 (Xmdn與Ymdn的計算，因現在適用Puzzle_c做計算，經緯度計算會有問題)
 Rfixed = imref2d(size(Puzzle_c));
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
  imshowpair(Puzzle_c,Rfixed,registered2,Rregistered,'blend');
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
  
%% 計算匹配照片經緯度
 % 計算實拍照片經緯度
  % 計算3*3拼接地圖最左上角那個點的像素點(須想一下邊界要怎麼處理)
  CropMap_TL_Pix = [CropSize*(block_X-2), CropSize*(block_Y-2)];   % [X, Y]
  % 計算特徵匹配照片的經緯度
  Estimated_Lat = TL_coor(2) - (CropMap_TL_Pix(2) + Ymdn) * Lat_per_pix;  % 因Xmdn/Ymdn是在局部(Puzzle_c)座標下計算出，須加上那張切割地圖的偏移量
  Estimated_Lon = TL_coor(1) + (CropMap_TL_Pix(1) + Xmdn) * Lon_per_pix; 
  Estimated_Pos = [Estimated_Lat Estimated_Lon];
  disp(['演算法計算之經緯度 -> 經度: ', num2str(Estimated_Lon, '%.8f'), ', 緯度: ', num2str(Estimated_Lat, '%.8f')]);
  disp(['照片實際之經緯度 -> 經度: ', num2str(AerialP_Lon, '%.8f'), ', 緯度: ', num2str(AerialP_Lat, '%.8f')]);
  
 % 將估算經緯度與實際經緯度畫出來
  figure();
  imshow(img1)
  hold on
  plot(CropMap_TL_Pix(1) + Xmdn , CropMap_TL_Pix(2) + Ymdn , 'bo', 'MarkerSize', 12, 'LineWidth', 2.5)
  
  X_shot_point = (AerialP_Lon - TL_coor(1)) / Lon_per_pix;
  Y_shot_point = (TL_coor(2) - AerialP_Lat) / Lat_per_pix;
  hold on
  plot(X_shot_point, Y_shot_point, 'ro', 'MarkerSize', 12, 'LineWidth', 2.5)
  hold off
  legend('Calculated Position', 'Real Position', 'FontSize', 16)
  