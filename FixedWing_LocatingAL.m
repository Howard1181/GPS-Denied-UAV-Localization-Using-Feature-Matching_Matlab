%% 以定翼機包含速度與方向的方式運行影像定位演算法
clear; clc; close all;

%% 讀取圖片
 % 提取TIF照片經緯度資訊(這裡在實際運用中應事前處理，先儲存在嵌入式系統內存中)
 % 逢甲
%  GEOinfo = geotiffinfo('離線地圖/逢甲/999.tif');
%  TL_coor = [GEOinfo.BoundingBox(1), GEOinfo.BoundingBox(4)]; % 地圖四個角經緯度，經度先
%  TR_coor = [GEOinfo.BoundingBox(2), GEOinfo.BoundingBox(4)];
%  BL_coor = [GEOinfo.BoundingBox(1), GEOinfo.BoundingBox(3)];
%  BR_coor = [GEOinfo.BoundingBox(2), GEOinfo.BoundingBox(3)];
 % 豐原
 TL_coor = [120.7157564, 24.2630274]; % 地圖四個角經緯度，經度先
 TR_coor = [120.7223255, 24.2630537];
 BL_coor = [120.7158070, 24.2573680];
 BR_coor = [120.7223298, 24.2573723];
 
 % 畫飛機
 Scale_factor_plane = 350;
 PlaneX = [ -0.1   0.1   0.05  0.05  0.3  0.05  0    -0.05  -0.3   -0.05  -0.05 ] .* Scale_factor_plane ;
 PlaneY = [ -0.1  -0.1  -0.05  0     0    0.05  0.3   0.05   0      0     -0.05 ] .* Scale_factor_plane;
 
%  fill(PlaneX, PlaneY,'b')
 
 % 定義地圖資訊
 
 % 逢甲
%  Map_Width = 3677;
%  Map_Height = 4142;
 
 % 豐原
 Map_Width = 4103; 
 Map_Height = 3894;

 CropSize = 300; % 豐原
%  CropSize = 350; % 逢甲
 numrow = ceil(Map_Width / CropSize) * CropSize; % 將numrow往上取到100的倍數
 numcol = ceil(Map_Height / CropSize) * CropSize;
 % 計算一個像素是多少經緯度
 Lon_per_pix = (round(TR_coor(1),8) - round(TL_coor(1),8)) / Map_Width;  % 這些應要都是已知值
 Lat_per_pix = (round(TL_coor(2),8) - round(BL_coor(2),8)) / Map_Height;
 
 % 讀取空拍照片
 folderPath_camera = 'photo';
 % 列出資料夾中所有照片
 imageFiles_camera = dir(fullfile(folderPath_camera, '*.jpg'));

 % 檢查是否有照片
 if isempty(imageFiles_camera)
     error('空拍資料夾中沒有找到照片');
 end
 
 %  讀取照片詳細資料
 PhoInfo_excelFile = 'fdata\FongYuanLog2.csv';
 data = readtable(PhoInfo_excelFile); % 使用 readtable 讀取 Excel
 Excel_dataName_cam = data{:, 1}; % 提取第一列的文件名

%% 建立無人機二維移動模型
 % 根據實際點兩點直線當作無人機移動軌跡
 % 根據實際點當作拍照地點並回傳圖片給演算法
 
 % 打開離線地圖並將拍攝座標標上去
 figure(1)
 folderPath_map = ('map\8KUHD_FongYuan2.png'); % 豐原
 img1 = imread(folderPath_map);
 imshow(img1)
 data = readtable(PhoInfo_excelFile); % 使用 readtable 讀取 Excel
 AerialP_Lon = data{:, 18};  % 用()提取table,用{}提取實際值
 AerialP_Lat = data{:, 17};
 Flight_Yaw = data{:, 9};
 Flight_Pitch = data{:, 10};
 AerialP_pixX = abs((round(AerialP_Lon,8) - round(BL_coor(1),8))) / Lon_per_pix; 
 AerialP_pixY = abs((round(AerialP_Lat,8) - round(BL_coor(2),8))) / Lat_per_pix;
 px = AerialP_pixX - 1;   
 py = Map_Height - AerialP_pixY - 1;
 hold on 
 plot(px(:), py(:), 'r-o', 'LineWidth', 2.5)
 hold on
 
 % 初始化參數
 k = 0; % 初始化重跑演算法旗標
 discard_cn = 0;  % 初始化放棄匹配次數旗標
 t_end = 100000; 
 switch_flag = 0; % 拍攝點到達判斷，因為目前是模擬，上機則是用失去GPS的時機判斷
 i = 1;
 uav_V = 12; 
 uav_px = 1;
 uav_py = 1;
 init_uavp = [uav_px, uav_py];
 Flight_Psi = atan2d(py(1) - init_uavp(2), px(1) - init_uavp(2));
 
 % 飛機根據psi轉向
 R = [cosd(Flight_Psi-90) -sind(Flight_Psi-90); sind(Flight_Psi-90) cosd(Flight_Psi-90)];
 RotatedPlane = R * [PlaneX;  PlaneY];
 PlaneX = RotatedPlane(1,:);
 PlaneY = RotatedPlane(2,:);
 
%  h = plot(init_uavp(1), init_uavp(2), 'bo', 'LineWidth', 3)
 h = fill(init_uavp(1) + PlaneX, init_uavp(2) + PlaneY,'y');
 

 
 % 在每點拍攝點間分成小等分，並以一質點代表飛機
 for dt = 1:1:t_end  % 模擬版本 - 以迴圈實現無人機動畫移動，並在中間用if判斷到達拍攝點並運行演算法

     Vx = uav_V * cosd(Flight_Psi); % AL_psi = +- 180
     Vy = uav_V * sind(Flight_Psi);
     
     uav_px = uav_px + Vx;
     uav_py = uav_py + Vy;
      
     set(h, 'XData', uav_px + PlaneX, 'YData', uav_py + PlaneY);
     drawnow
    
     if sqrt( (px(i)-uav_px)^2 + (py(i)-uav_py)^2 ) < 30
         disp('此處為拍攝點')
         switch_flag = 1;
     end
     
     if switch_flag == 1
         disp('更新航向')
         temp_Psi = Flight_Psi;
         Flight_Psi = atan2d(py(i+1)-py(i), px(i+1)-px(i));
         % 飛機根據psi轉向
         R = [cosd(Flight_Psi-temp_Psi) -sind(Flight_Psi-temp_Psi); sind(Flight_Psi-temp_Psi) cosd(Flight_Psi-temp_Psi)];
         RotatedPlane = R * [PlaneX;  PlaneY];
         PlaneX = RotatedPlane(1,:);
         PlaneY = RotatedPlane(2,:);
%          switch_flag = 0;
        % 迴圈控制
        i = i + 1;
        
     end

    %%
     % 地圖圖片裁切(若有換底圖的話)
     %  MapCrop()

 %%
     % 判斷進入演算法
     % 用if才比較符合實際情況，無人機繼續飛但演算法重算，後續匹配位置依無人機速度方向推算
     % while 接收到相機傳送之照片資料，立即進行處理，若拍下一張時，本張還沒匹配完的機制?
     while switch_flag == 1 
        
        j = i - 1;
        % 執行影像定位演算法
        % 逐張讀取空拍照片
        disp(['目前照片: ', num2str(j), '.jpg'])

        tic
        img2Path = fullfile(folderPath_camera, imageFiles_camera(j).name);
        img2 = imread(img2Path);

        % 讀取照片資訊(對答案用)
        % 使用 fileparts 提取文件夾路徑、文件名與副檔名
%         [fileFolder, imageFileName, ext] = fileparts(img2Path);
%         img2Name = [imageFileName , ext];
% 
%         idx = find(strcmp(Excel_dataName_cam, img2Name)); % 找到文件名的位置
% 
%         if ~isempty(idx)
%           AerialP_Lon = data{idx, 10};  % 用()提取table,用{}提取實際值
%           AerialP_Lat = data{idx, 9};
%           Flight_Yaw = data{idx, 12};
%         else
%           disp('未找到對應文件名');
%         end
% 
%         % 根據無人機失去GPS前的選取週遭離線地圖
%         % 計算pixel對應的圖片(轉換到imshow()座標)
%         px = abs(AerialP_pixX - 1);   
%         py = abs(Map_Height - AerialP_pixY - 1);
       
        % 這裡改用switch_flag=1時的無人機位置當作拼接中心
        block_X = ceil(abs(uav_px) / CropSize);
        block_Y = ceil(abs(uav_py) / CropSize);
        picture_num = (numcol/CropSize)*(block_X-1) + block_Y; % 座標點會落在picture_num這張小圖片上，圖片編號由上而下、由左至右編號
        % 改變picture_num的方法須改成用預測的
        
        
        % 從資料夾中讀取picture_num對應的小圖
%         CropImg_folderPath = 'Cropped_Map'; % 逢甲
        CropImg_folderPath = 'cmap'; % 豐原
        CropimgName = sprintf('%d.png', picture_num);
        CropimgPath = fullfile(CropImg_folderPath, CropimgName);
        if exist(CropimgPath, 'file')
            selectedImg = imread(CropimgPath);
        else
            disp('指定的圖片檔案不存在');
        end

        % 拼接小圖週遭地圖
        xStepNum = floor((numrow-CropSize)/CropSize+1); % 朝負無窮方向取整，寬度方向block移動的次數
        yStepNum = floor((numcol-CropSize)/CropSize+1); % 朝負無窮方向取整，長度方向block移動的次數
        xyStepNum = xStepNum*yStepNum;

        % 拼 3*3
        % 將小圖片及周遭的8張一併讀取出來，考量可能會抓到邊界或是計算到的圖片編號不在範圍內的情況，若發生這種現象則用零取代
        if (picture_num-yStepNum-1)<1 || (picture_num-yStepNum-1)>xyStepNum
            MapNum1 = zeros([CropSize,CropSize,3]);
        else
            MapNum1 = im2double(imread(['cmap\',num2str(picture_num-yStepNum-1),'.png']));
        end

        if (picture_num-yStepNum)<1 || (picture_num-yStepNum)>xyStepNum
            MapNum2 =zeros([CropSize,CropSize,3]);
        else
            MapNum2 = im2double(imread(['cmap\',num2str(picture_num-yStepNum),'.png']));
        end

        if (picture_num-yStepNum+1)<1 || (picture_num-yStepNum+1)>xyStepNum
            MapNum3 =zeros([CropSize,CropSize,3]);
        else
            MapNum3 = im2double(imread(['cmap\',num2str(picture_num-yStepNum+1),'.png']));
        end

        if (picture_num-1)<1 || (picture_num-1)>xyStepNum
            MapNum4 =zeros([CropSize,CropSize,3]);
        else
            MapNum4 = im2double(imread(['cmap\',num2str(picture_num-1),'.png']));
        end

        if (picture_num)<1 || (picture_num)>xyStepNum
            MapNum5 =zeros([CropSize,CropSize,3]);
        else
            MapNum5 = im2double(imread(['cmap\',num2str(picture_num),'.png']));
        end

        if (picture_num+1)<1 || (picture_num+1)>xyStepNum
            MapNum6 =zeros([CropSize,CropSize,3]);
        else
            MapNum6 = im2double(imread(['cmap\',num2str(picture_num+1),'.png']));
        end

        if (picture_num+yStepNum-1)<1 || (picture_num+yStepNum-1)>xyStepNum
            MapNum7 =zeros([CropSize,CropSize,3]);
        else    
            MapNum7 = im2double(imread(['cmap\',num2str(picture_num+yStepNum-1),'.png']));
        end

        if (picture_num+yStepNum)<1 || (picture_num+yStepNum)>xyStepNum
            MapNum8 =zeros([CropSize,CropSize,3]);
        else    
            MapNum8 = im2double(imread(['cmap\',num2str(picture_num+yStepNum),'.png']));
        end

        if (picture_num+yStepNum+1)<1 || (picture_num+yStepNum+1)>xyStepNum
            MapNum9 =zeros([CropSize,CropSize,3]);
        else
            MapNum9 = im2double(imread(['cmap\',num2str(picture_num+yStepNum+1),'.png']));
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
        Puzzle_c = [MapNum1 MapNum4 MapNum7; MapNum2 MapNum5 MapNum8; MapNum3 MapNum6 MapNum9]; 

        % 把拼好的圖片顯示出來
        figure(2);
        subplot(2,2,1);
        imshow(double(Puzzle_c));
        title('拼接完後的圖片[3*3]')
        subplot(2,2,3);
        imshow(img2);
        title('無人機空拍圖')

        % 並把圖片儲存起來
        PictureName = strcat('puzzle/Puzzle99.png'); % 將合成圖片儲存至想要的位置並命名
        imwrite(double(Puzzle_c),PictureName);

        % 調整圖片
        % 調整img解析度
        sfactor2 = 0.3; %豐原 
%         sfactor2 = 0.1; % 逢甲 
        if k == 1 && discard_cn >= 3
            sfactor2 = 0.2;
            disp('調整img2解析度');
        end

        img2 = imresize(img2, sfactor2);
        img2 = imrotate(img2, -Flight_Yaw(j));

        % 圖片轉灰階
        Gray_img1 = rgb2gray(Puzzle_c);
        Gray_img2 = rgb2gray(img2);
        Gray_img2 = histeq(Gray_img2);

        % 判斷是否需調整亮度/對比
        Gray_img2_db = double(Gray_img2) / 255;  % 將uint8圖片轉換成雙精度
        mean1 = mean(Gray_img1(:));
        mean2 = mean(Gray_img2_db(:));
        std1 = std(Gray_img1(:));
        std2 = std(Gray_img2_db(:));
        disp(['Image 1 亮度: ', num2str(mean1), ', 對比度: ', num2str(std1)]);
        disp(['Image 2 亮度: ', num2str(mean2), ', 對比度: ', num2str(std2)]);

        if (mean1 > 0.2 && abs(mean1-mean2) > 0.2 && abs(mean1-mean2) < 0.3) || (mean1 > 0.2 && abs(std1-std2) > 0.15 && abs(std1-std2) < 0.25) || discard_cn >= 2
            Gray_img2 = imhistmatch(Gray_img2, Gray_img1); % 光度均一化
            Gray_img2_db = double(Gray_img2) / 255;
            mean1 = mean(Gray_img1(:));
            mean2 = mean(Gray_img2_db(:));
            std1 = std(Gray_img1(:));
            std2 = std(Gray_img2_db(:));
            disp('亮度/對比差異大，執行光度均一');
            disp(['Image 1 均一化亮度: ', num2str(mean1), ', 均一化對比度: ', num2str(std1)]);
            disp(['Image 2 均一化亮度: ', num2str(mean2), ', 均一化對比度: ', num2str(std2)]);
        end

        % 特徵檢測
        points_img1 = detectKAZEFeatures(Gray_img1);
        points_img2 = detectKAZEFeatures(Gray_img2);
        points_img1 = points_img1.selectStrongest(15000);
        points_img2 = points_img2.selectStrongest(15000);

        % 描述符(可嘗試如何實作 Contextual similarity)
        [f1, vpts1] = extractFeatures(Gray_img1, points_img1);
        [f2, vpts2] = extractFeatures(Gray_img2, points_img2);

        % 進行匹配(參數要隨匹配狀況調整?)
        indexPairs = matchFeatures(f1, f2, 'MatchThreshold', 20, 'MaxRatio', 0.7) ; 
        matched_pts1 = vpts1(indexPairs(:, 1));
        matched_pts2 = vpts2(indexPairs(:, 2));

        % RANSAC演算法
        [tform, inlierimg2Points, inlierimg1Points] = estimateGeometricTransform(matched_pts2, matched_pts1, 'similarity', 'MaxNumTrials', 6000, 'MaxDistance', 10);
        toc

        % 匹配結果
%         figure('Name','對左右圖像分別進行特徵點提取')
%         figure(3)
        subplot(2,2,2);
        showMatchedFeatures(Gray_img1,Gray_img2,matched_pts1,matched_pts2,'montage');
        legend('matched points 1','matched points 2'); 
        title("特徵值匹配結果")
        subplot(2,2,4);
        showMatchedFeatures(Gray_img1,Gray_img2,inlierimg1Points,inlierimg2Points,'montage');
        legend('matched points 1','matched points 2');
        title("去除錯誤的匹配點後的特徵值匹配結果")

        % 重新匹配機制
        if (max(abs(tform.T(:))) > 2000 || inlierimg2Points.Count < 6) && discard_cn < 4
%             close (2);
%             close (3);

            % 清空圖內容但不關掉視窗
            figure(2)
            clf;

            k = 1; % 重新匹配旗標
            discard_cn = discard_cn +1;
            disp('匹配程度不佳，調整參數重新匹配');
            continue;
        elseif discard_cn >= 4
            % 將經緯度 Pixel 更新為演算法計算出來的 ( 目前作弊用實際經緯度代替，等穩定後換成 estimated_Lon/Lat )
%             AerialP_Lon = data{idx+1, 10};  % 用下一張照片位置進行拼接(也是作弊)，上機可用速度向量預測下一次迴圈位置
%             AerialP_Lat = data{idx+1, 9};
%             AerialP_pixX = abs((round(AerialP_Lon,8) - round(BL_coor(1),8))) / Lon_per_pix;  % 加abs目的為防止邊界誤判
%             AerialP_pixY = abs((round(AerialP_Lat,8) - round(BL_coor(2),8))) / Lat_per_pix;
            k = 0;
            discard_cn = 0;
            switch_flag = 0; % 重設拍攝點旗標

            % 清理內存
            figure(2)
            clf;
            clear points_img1 points_img2 f1 vpts1 f2 vpts2 tform inlierimg2Points inlierimg1Points registered2 Registered
            disp('捨棄本張照片，匹配下一張');
            break; % 跳出 while loop
        else
            % 疊圖
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
            figure(3)
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
            hold off
        end 

        % 計算匹配結果照片經緯度
        CropMap_TL_Pix = [CropSize*(block_X-2), CropSize*(block_Y-2)]; 

        % 計算特徵匹配照片的經緯度
        Estimated_Lat = TL_coor(2) - (CropMap_TL_Pix(2) + Ymdn ) * Lat_per_pix;  % -100 -150作弊用不該存在
        Estimated_Lon = TL_coor(1) + (CropMap_TL_Pix(1) + Xmdn ) * Lon_per_pix; 
%         Estimated_Pos = [Estimated_Lat Estimated_Lon];
        
        % 鏡頭角度校正
        uav_altitude = 120; % m
        gimbal_pitch = 8.5; % degree
        if j >= 46
            gimbal_pitch = 0; % 豐原拍攝的時候在第46張照片後鏡頭皆調整為正下
        end
        total_pitch = gimbal_pitch + Flight_Pitch(j);
        delta_d = 120 * tand(total_pitch);
        delta_x = delta_d * sind(Flight_Yaw(j)); % m
        delta_y = -delta_d * cosd(Flight_Yaw(j)); % m
        delta_y2lat = 0.00898 * delta_y/1000; % 1000公尺影響 0.00898緯度
        delta_x2lon = (1000/(111320*cosd(24.26))) * delta_x/1000; % 1000公尺影響 (1000/(111320*cosd(當地緯度))) 的經度，葫蘆墩公園緯度為24.26
        Estimated_Pos = [Estimated_Lat+delta_y2lat  Estimated_Lon-delta_x2lon];

        % 計算實際經緯度與估測經緯度誤差距離
        deltaLon = Estimated_Pos(2) - AerialP_Lon(j);
        deltaLat = Estimated_Pos(1) - AerialP_Lat(j);

        % Haversine 經緯度轉距離
        a = sind(deltaLat / 2)^2 + cosd(Estimated_Pos(1)) * cosd(AerialP_Lat(j)) * sind(deltaLon / 2)^2;
        c = 2 * atan2(sqrt(a), sqrt(1 - a));
        R = 6371; % km
        d = 1000 * R * c; % m 

        % 判斷d大小，若太大重跑回圈(此為開發階段用)，重新匹配時若有新的照片進來該如何解決(告訴相機先不要照?)，預測位置邏輯又該如何解決
        if d > 30 && discard_cn < 4
           k = 1;
           
           figure(2)
           clf;
           figure(3)
           clf;
           discard_cn = discard_cn +1;
           disp(['距離誤差過大: ', num2str(d, '%.2f')]);
           continue;
        end

        Bias(j) = d;

        % 顯示
        disp(['演算法計算之經緯度 -> 經度: ', num2str(Estimated_Pos(2), '%.8f'), ', 緯度: ', num2str(Estimated_Pos(1), '%.8f')]);
        disp(['照片實際之經緯度 -> 經度: ', num2str(AerialP_Lon(j), '%.8f'), ', 緯度: ', num2str(AerialP_Lat(j), '%.8f')]);
        disp(['估測距離差距: ', num2str(d, '%.2f'), ' 公尺']);

        % 將估算經緯度與實際經緯度畫出來
        X_shot_point = (AerialP_Lon(j) - TL_coor(1)) / Lon_per_pix;
        Y_shot_point = (TL_coor(2) - AerialP_Lat(j)) / Lat_per_pix;
        
        % 將演算法計算之拍攝點畫在figure1
        Estimated_pixX = (round(Estimated_Pos(2),8) - round(BL_coor(1),8)) / Lon_per_pix; 
        Estimated_pixY = Map_Height - (round(Estimated_Pos(1),8) - round(BL_coor(2),8)) / Lat_per_pix; 
        
        figure(1)
        plot(Estimated_pixX, Estimated_pixY, 'xy', 'LineWidth', 2.5) % -100 -150作弊用不該存在
        drawnow;
        
        % 儲存實際與估算畫素值，以便執行後畫圖
        Real_points(j,:) = [X_shot_point, Y_shot_point];
        Estimated_points(j,:) = [Estimated_pixX, Estimated_pixY];

%         hold on
%         plot(X_shot_point, Y_shot_point, 'ro', 'MarkerSize', 12, 'LineWidth', 2.5)
%         hold off
%         legend('Calculated Position', 'Real Position', 'FontSize', 16)

        % 將經緯度 Pixel 更新為演算法計算出來的 ( 目前作弊用實際經緯度代替，等穩定後換成 estimated_Lon/Lat )
%         AerialP_Lon = data{idx+1, 10};  % 用下一張照片位置進行拼接(也是作弊)，上機可用速度向量預測下一次迴圈位置
%         AerialP_Lat = data{idx+1, 9};
%         AerialP_pixX = abs((round(AerialP_Lon,8) - round(BL_coor(1),8))) / Lon_per_pix; 
%         AerialP_pixY = abs((round(AerialP_Lat,8) - round(BL_coor(2),8))) / Lat_per_pix;

        % 迴圈結束清理內存
        pause(2)
        figure(2)
        clf;
        figure(3)
        clf;

        clear points_img1 points_img2 f1 vpts1 f2 vpts2 tform inlierimg2Points inlierimg1Points registered2 Registered

        % 監控內存
       memoryInfo = memory;
       disp(['Iteration ', num2str(j), ': Memory in use: ', num2str(memoryInfo.MemUsedMATLAB / 1e6), ' MB ']);
 
       % 迴圈控制

       discard_cn = 0;  % 重設丟棄旗標
       switch_flag = 0; % 重設拍攝點旗標
       
     end
     
     
     
     
 end
 
 
 
 
 
 
 
 
 