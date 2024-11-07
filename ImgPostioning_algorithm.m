%% 影像定位演算法
% 自動將照片與GPS、iMU資料同步
% 自動讀取照片並跑演算法
% 應用此演算法前需準備的資料:
% 1、 圖資資料
% 2、 空拍照片
% 3、 空拍照片對應之位置姿態EXCEL表


clc; clear; close all;

%% 讀取區域地圖 
  Map = imread('G_FullMap5.png');
  GrayMap = rgb2gray(Map);
  MapSize = size(GrayMap);
% 找地圖特徵點
  pointsMap = detectSURFFeatures(GrayMap, 'MetricThreshold', 500);
  [Msize, m] = size(pointsMap.Location);
%   figure(1)
%   imshow(GrayMap); hold on;
%   plot(pointsMap(1:2:Msize,:),'showOrientation',true);
%   hold off;
% 計算地圖描述向量
  [fMap, vptsMap] = extractFeatures(GrayMap, pointsMap);
% 將地圖轉換成2為座標
  Rfixed = imref2d(size(Map));
%% 開發階段，每一秒讀取一張照片
% 指定資料夾
  folderPath = 'IMG_1s1p';

% 列出資料夾中所有照片
  imageFiles = dir(fullfile(folderPath, '*.jpg'));
  
% 檢查是否有照片
  if isempty(imageFiles)
    error('資料夾中沒有找到照片');
  end
  
% 每1秒讀取一張照片
  for i = 1:length(imageFiles)
    try  % 使用try-catch 捕捉錯誤，當有錯誤時顯示錯誤並跳至下一個迴圈
           
    % 建立照片路徑索引
      imagePath = fullfile(folderPath, imageFiles(i).name);

    % 讀取照片
      img = imread(imagePath);
      img = imresize(img, 0.2);


    %% 讀取照片姿態 (這邊先跳過，待自動同步GPS/ATT程式完成) 


    %% 圖片轉灰階
      GrayImg = rgb2gray(img);

    %% 找圖片特徵點
      pointsImg = detectSURFFeatures(GrayImg); 
      [Isize, n] = size(pointsImg.Location);
      figure(2)
      imshow(GrayImg); hold on;
      plot(pointsImg(1:2:Isize,:),'showOrientation',true); 
      hold off

    % 計算描述向量
      [fImg, vptsImg] = extractFeatures(GrayImg, pointsImg);

    %% 提取特徵點位置
%       indexPairs = matchFeatures(fMap, fImg, 'Prenormalized', true) ; 
      indexPairs = matchFeatures(fMap, fImg, 'Method', 'Approximate') ;
      matched_pts1 = vptsMap(indexPairs(:, 1));
      matched_pts2 = vptsImg(indexPairs(:, 2));

    %% 用MSAC演算法剔除錯誤匹配點
      [tform, inlierimg2Points, inlierimg1Points] = estimateGeometricTransform(matched_pts2, matched_pts1, 'similarity', 'MaxDistance', 4);
      figure(3)
      subplot(2,1,1)
      showMatchedFeatures(GrayMap,GrayImg,matched_pts1,matched_pts2,'montage');
      legend('matched points 1','matched points 2'); 
      title("特徵值匹配結果")

      subplot(2,1,2)
      showMatchedFeatures(GrayMap,GrayImg,inlierimg1Points,inlierimg2Points,'montage');
      legend('matched points 1','matched points 2');
      title("去除錯誤的匹配點後的特徵值匹配結果")
    %% 進行圖像合併
      [registered2, registered1] = imwarp(img, tform);
      boxPolygon = [1, 1;... % 左上
        size(img, 2), 1; ... % 右上
        size(img, 2), size(img, 1); ... % 右下
        1, size(img, 1); ... % 左下
        1, 1]; % 重複左下，才可得到一個閉區間的多邊形
      % 將多邊形變換到目標圖片上，變換的結果表示了物體的位置
        newBoxPolygon = transformPointsForward(tform, boxPolygon);
      % 顯示被檢測到的物體
        figure(4)
        imshowpair(Map,Rfixed,registered2,registered1,'blend');
    

      % 暫停1秒，實際上要看演算法跑多久
        pause(1);
     
       % 清除不再需要的變量，釋放內存 (這步很重要，要注意)
        clear img GrayImg pointsImg fImg vptsImg indexPairs matched_pts1 matched_pts2 tform inlierimg2Points inlierimg1Points;
        close(2)
        close(3)
        close(4)
        
      catch ME
        % 捕獲錯誤並顯示錯誤信息
        disp(['發生錯誤: ', ME.message]);
        % 繼續執行下一個迴圈
      continue;
    end
    
    % 監控內存
    memoryInfo = memory;
    disp(['Iteration ', num2str(i), ': Memory in use: ', num2str(memoryInfo.MemUsedMATLAB / 1e6), ' MB']);
  end
  
close all;
  