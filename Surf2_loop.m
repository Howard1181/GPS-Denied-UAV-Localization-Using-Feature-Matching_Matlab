%% 無GPS環境平面影像定位演算法 Matlab版本
% 1. 讀取大地圖(無裁切)
% 2. 讀取大地圖地理資訊，四個角的確切經緯度
% 3. 手工鍵入大地圖圖像資訊，長 X 寬，與裁切地圖的大小，並計算一個畫素是多少經緯度
% 4. 大地圖裁切(一開始就要做，後續有換大地圖再跑)，高斯濾雜訊後匹配較好
% 5. 讀取首張照片之GPS位置當作初始值，上機版本需獲取無人機在失去GPS訊號前的GPS位置
% 6. 依序讀取資料夾中拍攝好之空拍照片
% 7. 調整讀取圖片之參數，解析度/亮度/對比等等。需尋找能自動依據圖片狀況調整之方法!或需在辨識效果不好時(閥值)重新調整解析度或其他參數
% 8. 根據第5.之GPS座標選取相對應之裁切地圖，並拼接成5*5之地圖，後續需依據估測位置重新拼接5*5地圖
% 9. 5*5地圖與空拍圖灰階化，並將空拍圖調整亮度、對比與地圖一致(這裡有時調比較好不調比較好，也需想調整機制)
% 續9. 在邊界時地圖亮度與對比會受影響，這點需考慮
% 10. 進行兩張灰階化的圖片的特徵檢測，SURF較快效果較差，KAZE較慢效果較好
% 11. 描述子定義，特徵匹配(參數調整需盡量固定，或有自適應的方法)
% 12. 錯誤匹配檢測RANSAC，參數調整需盡量固定
% 13. 計算匹配照片中心點並計算對應畫素再轉成經緯度
% 14. 提取Excel表空照圖資訊(應在迴圈中完成)，此步驟為對答案用，在實飛上用不到
% 15. 比對實際照片經緯度與演算法計算經緯度
% 16. 根據匹配結果(tform, inlierimg1Points, Bias等數據)判斷此次匹配正確性，並決定是否調參(根據旗標k)重新匹配
% 17. 若匹配次數達上限則放棄此次迭代的照片
% --------------------- Python/上機 注意事項 --------------------------------- %
% 如何選擇拼接照片中心picture_num，啟用演算法第一筆為失去GPS前的位置，後續依據演算法計算結果or根據無人機飛向速度向量預測軌跡?
% 在飛行中若匹配失敗的處理機制，連續匹配三次失敗跳下一張，預測位置該如何選取? 飛試該如何知道匹配錯誤三次，以甚麼資訊?

%% 以不同特徵提取方法嘗試
clear; clc; close all;

%% 讀取圖片
 img1 = imread('離線地圖/8K_FJU.png'); % 這裡在實際嵌入式運用中應讀取的是裁切完的照片
 
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
 % 計算一個像素是多少經緯度
 Lon_per_pix = (round(TR_coor(1),8) - round(TL_coor(1),8)) / Map_Width;  % 這些應要都是已知值
 Lat_per_pix = (round(TL_coor(2),8) - round(BL_coor(2),8)) / Map_Height;
 
%%  
% 地圖圖片裁切(若有換底圖的話)
%  MapCrop()

 % 計算丟失GPS前經緯度在2離線地圖哪個Pixel，第一次迭代之後需更新成演算法計算出來的照片的位置
 % 此為初始化，[120.6464858, 24.18187389] 在飛試中需判斷失去GPS並開始執行演算法，接收失去GPS前的GPS座標當作第一次匹配的目標 

 AerialP_pixX = abs((round(120.6464858, 8) - round(BL_coor(1),8))) / Lon_per_pix;
 AerialP_pixY = abs((round(24.18187389, 8) - round(BL_coor(2),8))) / Lat_per_pix;
 
%%   
 % 指定資料夾
  folderPath = '逢甲空拍';

 % 列出資料夾中所有照片
  imageFiles = dir(fullfile(folderPath, '*.jpg'));
  
 % 檢查是否有照片
  if isempty(imageFiles)
    error('資料夾中沒有找到照片');
  end 
  
%%  讀取照片詳細資料
 PhoInfo_excelFile = 'Log.csv';
 data = readtable(PhoInfo_excelFile); % 使用 readtable 讀取 Excel
 Excel_dataName = data{:, 1}; % 提取第一列的文件名
 
%% 影像定位演算法迴圈
i = 1;   % 初始化i，逐張讀去照片用，上機則讀取當下拍的照片
k = 0; % 初始化重跑演算法旗標
discard_cn = 0;  % 初始化放棄匹配次數旗標

while i <= length(imageFiles)
%% 逐張讀取空拍照片    
 disp(['目前照片: ', num2str(i), '.jpg'])
 tic
 img2Path = fullfile(folderPath, imageFiles(i).name);
 img2 = imread(img2Path);
 
 % 讀取照片資訊(對答案用)
 % 使用 fileparts 提取文件夾路徑、文件名與副檔名
 [fileFolder, imageFileName, ext] = fileparts(img2Path);
 img2Name = [imageFileName , ext];
 
 idx = find(strcmp(Excel_dataName, img2Name)); % 找到文件名的位置
  
 if ~isempty(idx)
   AerialP_Lon = data{idx, 10};  % 用()提取table,用{}提取實際值
   AerialP_Lat = data{idx, 9};
   Flight_Yaw = data{idx, 12};
 else
   disp('未找到對應文件名');
 end

%% 根據無人機失去GPS前的選取週遭離線地圖
 % 根據座標進行週遭地圖拼接
 % 計算pixel對應的圖片(轉換到imshow()座標)
 px = abs(AerialP_pixX - 1);   
 py = abs(Map_Height - AerialP_pixY - 1);
 block_X = ceil(px / CropSize);
 block_Y = ceil(py / CropSize);
%   block_X = floor(px / CropSize);
%   block_Y = floor(py / CropSize);
  
 picture_num = (numcol/CropSize)*(block_X-1) + block_Y; % 座標點會落在picture_num這張小圖片上，圖片編號由上而下、由左至右編號
%   hold on
%   plot(AerialP_pixX-1, Map_Height-AerialP_pixY-1, 'ro', 'MarkerSize', 15)
%   hold off
  
 % 從資料夾中讀取picture_num對應的小圖
 CropImg_folderPath = 'Cropped_Map';
 CropimgName = sprintf('%d.png', picture_num);
 CropimgPath = fullfile(CropImg_folderPath, CropimgName);
 if exist(CropimgPath, 'file')
    % 讀取圖片
   selectedImg = imread(CropimgPath);  % 好像不用這個? 
   
 % 顯示圖片
%     figure()
%     imshow(selectedImg);
 else
   disp('指定的圖片檔案不存在');
 end
 
 % 拼接小圖週遭地圖
 xStepNum = floor((numrow-CropSize)/CropSize+1); % 朝負無窮方向取整，寬度方向block移動的次數
 yStepNum = floor((numcol-CropSize)/CropSize+1); % 朝負無窮方向取整，長度方向block移動的次數
 xyStepNum = xStepNum*yStepNum;
%% 拼5*5
%  % 將小圖片及他周遭的24張一併讀取出來，考量可能會抓到邊界或是計算到的圖片編號不在範圍內的情況，若發生這種現象則用零取代
%  if (picture_num-2*yStepNum-2)<1 || (picture_num-2*yStepNum-2)>xyStepNum
%      MapNum1 = zeros([CropSize,CropSize,3]);
%  else
%      MapNum1 = im2double(imread(['Cropped_Map\',num2str(picture_num-2*yStepNum-2),'.png']));
%  end
% 
%  if (picture_num-2*yStepNum-1)<1 || (picture_num-2*yStepNum-1)>xyStepNum
%      MapNum2 = zeros([CropSize,CropSize,3]);
%  else
%      MapNum2 = im2double(imread(['Cropped_Map\',num2str(picture_num-2*yStepNum-1),'.png']));
%  end
% 
%  if (picture_num-2*yStepNum)<1 || (picture_num-2*yStepNum)>xyStepNum
%      MapNum3 = zeros([CropSize,CropSize,3]);
%  else
%      MapNum3 = im2double(imread(['Cropped_Map\',num2str(picture_num-2*yStepNum),'.png']));
%  end
% 
%  if (picture_num-2*yStepNum+1)<1 || (picture_num-2*yStepNum+1)>xyStepNum
%      MapNum4 = zeros([CropSize,CropSize,3]);
%  else
%      MapNum4 = im2double(imread(['Cropped_Map\',num2str(picture_num-2*yStepNum+1),'.png']));
%  end
% 
%  if (picture_num-2*yStepNum+2)<1 || (picture_num-2*yStepNum+2)>xyStepNum
%      MapNum5 = zeros([CropSize,CropSize,3]);
%  else
%      MapNum5 = im2double(imread(['Cropped_Map\',num2str(picture_num-2*yStepNum+2),'.png']));
%  end
% 
%  if (picture_num-yStepNum-2)<1 || (picture_num-yStepNum-2)>xyStepNum
%      MapNum6 = zeros([CropSize,CropSize,3]);
%  else
%      MapNum6 = im2double(imread(['Cropped_Map\',num2str(picture_num-yStepNum-2),'.png']));
%  end
% 
%  if (picture_num-yStepNum-1)<1 || (picture_num-yStepNum-1)>xyStepNum
%      MapNum7 = zeros([CropSize,CropSize,3]);
%  else    
%      MapNum7 = im2double(imread(['Cropped_Map\',num2str(picture_num-yStepNum-1),'.png']));
%  end
% 
%  if (picture_num-yStepNum)<1 || (picture_num-yStepNum)>xyStepNum
%      MapNum8 = zeros([CropSize,CropSize,3]);
%  else    
%      MapNum8 = im2double(imread(['Cropped_Map\',num2str(picture_num-yStepNum),'.png']));
%  end
% 
%  if (picture_num-yStepNum+1)<1 || (picture_num-yStepNum+1)>xyStepNum
%      MapNum9 = zeros([CropSize,CropSize,3]);
%  else
%      MapNum9 = im2double(imread(['Cropped_Map\',num2str(picture_num-yStepNum+1),'.png']));
%  end
% 
%  if (picture_num-yStepNum+2)<1 || (picture_num-yStepNum+2)>xyStepNum
%      MapNum10 = zeros([CropSize,CropSize,3]);
%  else
%      MapNum10 = im2double(imread(['Cropped_Map\',num2str(picture_num-yStepNum+2),'.png']));
%  end
% 
%  if (picture_num-2)<1 || (picture_num-2)>xyStepNum
%      MapNum11 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum11 = im2double(imread(['Cropped_Map\',num2str(picture_num-2),'.png']));
%  end
% 
%  if (picture_num-1)<1 || (picture_num-1)>xyStepNum
%      MapNum12 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum12 = im2double(imread(['Cropped_Map\',num2str(picture_num-1),'.png']));
%  end
% 
%  if (picture_num)<1 || (picture_num)>xyStepNum
%      MapNum13 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum13 = im2double(imread(['Cropped_Map\',num2str(picture_num),'.png']));
%  end
% 
%  if (picture_num+1)<1 || (picture_num+1)>xyStepNum
%      MapNum14 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum14 = im2double(imread(['Cropped_Map\',num2str(picture_num+1),'.png']));
%  end
% 
%  if (picture_num+2)<1 || (picture_num+2)>xyStepNum
%      MapNum15 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum15 = im2double(imread(['Cropped_Map\',num2str(picture_num+2),'.png']));
%  end
% 
%  if (picture_num+yStepNum-2)<1 || (picture_num+yStepNum-2)>xyStepNum
%      MapNum16 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum16 = im2double(imread(['Cropped_Map\',num2str(picture_num+yStepNum-2),'.png']));
%  end
% 
%  if (picture_num+yStepNum-1)<1 || (picture_num+yStepNum-1)>xyStepNum
%      MapNum17 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum17 = im2double(imread(['Cropped_Map\',num2str(picture_num+yStepNum-1),'.png']));
%  end
% 
%  if (picture_num+yStepNum)<1 || (picture_num+yStepNum)>xyStepNum
%      MapNum18 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum18 = im2double(imread(['Cropped_Map\',num2str(picture_num+yStepNum),'.png']));
%  end
% 
%  if (picture_num+yStepNum+1)<1 || (picture_num+yStepNum+1)>xyStepNum
%      MapNum19 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum19 = im2double(imread(['Cropped_Map\',num2str(picture_num+yStepNum+1),'.png']));
%  end
% 
%  if (picture_num+yStepNum+2)<1 || (picture_num+yStepNum+2)>xyStepNum
%      MapNum20 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum20 = im2double(imread(['Cropped_Map\',num2str(picture_num+yStepNum+2),'.png']));
%  end
% 
%  if (picture_num+2*yStepNum-2)<1 || (picture_num+2*yStepNum-2)>xyStepNum
%      MapNum21 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum21 = im2double(imread(['Cropped_Map\',num2str(picture_num+2*yStepNum-2),'.png']));
%  end
% 
%  if (picture_num+2*yStepNum-1)<1 || (picture_num+2*yStepNum-1)>xyStepNum
%      MapNum22 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum22 = im2double(imread(['Cropped_Map\',num2str(picture_num+2*yStepNum-1),'.png']));
%  end
% 
%  if (picture_num+2*yStepNum)<1 || (picture_num+2*yStepNum)>xyStepNum
%      MapNum23 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum23 = im2double(imread(['Cropped_Map\',num2str(picture_num+2*yStepNum),'.png']));
%  end
% 
%  if (picture_num+2*yStepNum+1)<1 || (picture_num+2*yStepNum+1)>xyStepNum
%      MapNum24 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum24 = im2double(imread(['Cropped_Map\',num2str(picture_num+2*yStepNum+1),'.png']));
%  end
% 
%  if (picture_num+2*yStepNum+2)<1 || (picture_num+2*yStepNum+2)>xyStepNum
%      MapNum25 =zeros([CropSize,CropSize,3]);
%  else
%      MapNum25 = im2double(imread(['Cropped_Map\',num2str(picture_num+2*yStepNum+2),'.png']));
%  end
% 
%  % 避免上下邊界的誤判
%  for  count = 1:1:xStepNum  
%     % 若 picture_num 的值位於第一排
%     if picture_num == (count-1)*yStepNum + 1
%         % 就把 MapNum1、MapNum6、MapNum11、MapNum16、MapNum21 用零矩陣替換掉
%         MapNum1 = zeros([CropSize,CropSize,3]);  
%         MapNum6 = zeros([CropSize,CropSize,3]);   
%         MapNum11 = zeros([CropSize,CropSize,3]);
%         MapNum16 = zeros([CropSize,CropSize,3]);
%         MapNum21 = zeros([CropSize,CropSize,3]);
%         
%         MapNum2 = zeros([CropSize,CropSize,3]);  
%         MapNum7 = zeros([CropSize,CropSize,3]);   
%         MapNum12 = zeros([CropSize,CropSize,3]);
%         MapNum17 = zeros([CropSize,CropSize,3]);
%         MapNum22 = zeros([CropSize,CropSize,3]);
%     end
%     
%     % 若 picture_num 的值位於第二排
%     if picture_num == (count-1)*yStepNum + 2
%         % 就把 MapNum1、MapNum6、MapNum11、MapNum16、MapNum21 用零矩陣替換掉
%         MapNum1 = zeros([CropSize,CropSize,3]);  
%         MapNum6 = zeros([CropSize,CropSize,3]);   
%         MapNum11 = zeros([CropSize,CropSize,3]);
%         MapNum16 = zeros([CropSize,CropSize,3]);
%         MapNum21 = zeros([CropSize,CropSize,3]);
%     end
%     
%     % 若 picture_num 的值位於最末排
%     if picture_num == (count)*yStepNum
%         % 就把 MapNum5、MapNum10、MapNum15、MapNum20、MapNum25 用零矩陣替換掉
%         MapNum5 = zeros([CropSize,CropSize,3]);  
%         MapNum10 = zeros([CropSize,CropSize,3]);   
%         MapNum15 = zeros([CropSize,CropSize,3]);
%         MapNum20 = zeros([CropSize,CropSize,3]);
%         MapNum25 = zeros([CropSize,CropSize,3]);
%         
%         MapNum4 = zeros([CropSize,CropSize,3]);  
%         MapNum9 = zeros([CropSize,CropSize,3]);   
%         MapNum14 = zeros([CropSize,CropSize,3]);
%         MapNum19 = zeros([CropSize,CropSize,3]);
%         MapNum24 = zeros([CropSize,CropSize,3]);
%     end
%     
%     % 若 picture_num 的值位於倒數第二排
%     if picture_num == (count)*yStepNum - 1
%         % 就把 MapNum5、MapNum10、MapNum15、MapNum20、MapNum25 用零矩陣替換掉
%         MapNum5 = zeros([CropSize,CropSize,3]);  
%         MapNum10 = zeros([CropSize,CropSize,3]);   
%         MapNum15 = zeros([CropSize,CropSize,3]);
%         MapNum20 = zeros([CropSize,CropSize,3]);
%         MapNum25 = zeros([CropSize,CropSize,3]);
%     end 
%  end
% 
%  % 再把25張圖片照順序拼回去
%  Puzzle_c =[MapNum1 MapNum6 MapNum11 MapNum16 MapNum21; 
%            MapNum2 MapNum7 MapNum12 MapNum17 MapNum22; 
%            MapNum3 MapNum8 MapNum13 MapNum18 MapNum23;
%            MapNum4 MapNum9 MapNum14 MapNum19 MapNum24;
%            MapNum5 MapNum10 MapNum15 MapNum20 MapNum25]; 
% 
%  % 把拼好的圖片顯示出來
%  figure('Name','拼接完後的圖片[5*5]');
%  imshow(double(Puzzle_c));
% 
%  % 並把圖片儲存起來
%  PictureName = strcat('小地圖拼接/Puzzle55.png'); % 將合成圖片儲存至想要的位置並命名
%  imwrite(double(Puzzle_c),PictureName);

%% 拼 3*3
% 將小圖片及他周遭的8張一併讀取出來，考量可能會抓到邊界或是計算到的圖片編號不在範圍內的情況，若發生這種現象則用零取代
if (picture_num-yStepNum-1)<1 || (picture_num-yStepNum-1)>xyStepNum
    MapNum1 = zeros([CropSize,CropSize,3]);
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
Puzzle_c = [MapNum1 MapNum4 MapNum7; MapNum2 MapNum5 MapNum8; MapNum3 MapNum6 MapNum9]; 

% 把拼好的圖片顯示出來
figure('Name','拼接完後的圖片[3*3]');
imshow(double(Puzzle_c));

% 並把圖片儲存起來
PictureName = strcat('小地圖拼接/Puzzle99.png'); % 將合成圖片儲存至想要的位置並命名
imwrite(double(Puzzle_c),PictureName);

%% 調整圖片
 % 調整img解析度
 sfactor2 = 0.1;  % 想想如何根據兩張照片性質自動調整scale?
 
 if k == 1 && discard_cn >= 3
    sfactor2 = 0.2;
    disp('調整img2解析度');
 end
 
 img2 = imresize(img2, sfactor2);
 % 調整圖片方向
 img2 = imrotate(img2, -Flight_Yaw); % 根據excel yaw值旋轉圖片，imrotate(img, angle)->逆時針轉angle度，上機理論上也可透過接收Yaw值做到

 % 顯示結果
 % figure('Name','讀取彩色圖片')
 % subplot(2,1,1)
 % imshow(img1)
 % title("離線地圖")
 % subplot(2,1,2)
 % imshow(img2)
 % title("無人機實拍畫面") 

%% 圖片轉灰階
  % 調整圖片
  Gray_img1 = rgb2gray(Puzzle_c);
  
%   Gray_img1 = histeq(Gray_img1);
%   Gray_img1 = imadjust(Gray_img1, [0 1], [0.1 0.9]);
%   Gray_img2 = adapthisteq(Gray_img2);

%   Gray_img1 = Gray_img1 - 0.12;
  Gray_img2 = rgb2gray(img2);
  Gray_img2 = histeq(Gray_img2);
  
  % Gray_img2 = adapthisteq(Gray_img2);
  
%% 判斷是否需調整亮度、對比
 % 計算兩張圖片的亮度（平均值）和對比度（標準差）
 hist_Gray_img1 = Gray_img1(Gray_img1 > 0);  % 似乎不太妥當
 Gray_img2_db = double(Gray_img2) / 255;  % 將uint8圖片轉換成雙精度
 mean1 = mean(Gray_img1(:));
 mean2 = mean(Gray_img2_db(:));
 std1 = std(Gray_img1(:));
 std2 = std(Gray_img2_db(:));

 % 顯示亮度與對比度的差異
 disp(['Image 1 亮度: ', num2str(mean1), ', 對比度: ', num2str(std1)]);
 disp(['Image 2 亮度: ', num2str(mean2), ', 對比度: ', num2str(std2)]);
 
 if (mean1 > 0.2 && abs(mean1-mean2) > 0.1 && abs(mean1-mean2) < 0.25) || (mean1 > 0.2 && abs(std1-std2) > 0.1 && abs(std1-std2) < 0.25) || discard_cn >= 2
    Gray_img2 = imhistmatch(Gray_img2, Gray_img1); % 光度均一化
    Gray_img2_db = double(Gray_img2) / 255;
    mean1 = mean(Gray_img1(:));
    mean2 = mean(Gray_img2_db(:));
    std1 = std(Gray_img1(:));
    std2 = std(Gray_img2_db(:));
    disp('光度均一化');
    disp(['Image 1 均一化亮度: ', num2str(mean1), ', 均一化對比度: ', num2str(std1)]);
    disp(['Image 2 均一化亮度: ', num2str(mean2), ', 均一化對比度: ', num2str(std2)]);
 end
 
  % 顯示結果
%   figure('Name','轉成灰階圖片')
%   subplot(2,1,1)
%   imshow(Gray_img1) 
%   title("地圖")
%   subplot(2,1,2)
%   imshow(Gray_img2)
%   title("無人機實拍地圖")

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

% 計算描述符
  [f1, vpts1] = extractFeatures(Gray_img1, points_img1);
  [f2, vpts2] = extractFeatures(Gray_img2, points_img2);

% 進行匹配
  indexPairs = matchFeatures(f1, f2, 'MatchThreshold', 20, 'MaxRatio', 0.7) ; 
  matched_pts1 = vpts1(indexPairs(:, 1));
  matched_pts2 = vpts2(indexPairs(:, 2));

% 去除錯誤點
  [tform, inlierimg2Points, inlierimg1Points] = estimateGeometricTransform(matched_pts2, matched_pts1, 'similarity', 'MaxNumTrials', 6000, 'MaxDistance', 10);
  % MaxDistance可試著調大一些
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

  
%% 判斷匹配結果，設置旗標

% 若特徵匹配點少於6個則當前照片不變重跑演算法
 if (max(abs(tform.T(:))) > 2000 || inlierimg2Points.Count < 6) && discard_cn < 4
    close(1)
    close(2)
    k = 1;
    discard_cn = discard_cn +1;
    disp('匹配程度不佳，調整參數重新匹配');
    continue;
 elseif discard_cn >= 4
    % 將經緯度 Pixel 更新為演算法計算出來的 ( 目前作弊用實際經緯度代替，等穩定後換成 estimated_Lon/Lat )
    AerialP_Lon = data{idx+1, 10};  % 用下一張照片位置進行拼接(也是作弊)，上機可用速度向量預測下一次迴圈位置
    AerialP_Lat = data{idx+1, 9};
    AerialP_pixX = abs((round(AerialP_Lon,8) - round(BL_coor(1),8))) / Lon_per_pix;  % 加abs目的為防止邊界誤判
    AerialP_pixY = abs((round(AerialP_Lat,8) - round(BL_coor(2),8))) / Lat_per_pix;
    k = 0;
    discard_cn = 0;
    i = i+1;
    % 清理內存
    close all;
    clear points_img1 points_img2 f1 vpts1 f2 vpts2 tform inlierimg2Points inlierimg1Points registered2 Registered
    disp('捨棄本張照片，匹配下一張');
    continue;
 else
%% 疊圖 (Xmdn與Ymdn的計算)
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
 end

%% 計算匹配照片經緯度
 % 計算實拍照片經緯度
  % 計算5*5拼接地圖最左上角那個點的像素點(需想一下邊界要怎麼處理)
  CropMap_TL_Pix = [CropSize*(block_X-2), CropSize*(block_Y-2)];   % [X, Y] 這邊需根據拼3*3or5*5決定，3*3 -2，5*5 -3
%   if block_X == 2 || block_X == 1   % -> 這邊有問題需處理
%       CropMap_TL_Pix(1) = 0;    
%   end
%   
%   if block_Y == 2 || block_Y == 1
%       CropMap_TL_Pix(2) = 0;    
%   end
  
  % 計算特徵匹配照片的經緯度
  Estimated_Lat = TL_coor(2) - (CropMap_TL_Pix(2) + Ymdn) * Lat_per_pix;  % 因 Xmdn/Ymdn 是在局部(Puzzle_c)座標下計算出，須加上那張切割地圖的偏移量
  Estimated_Lon = TL_coor(1) + (CropMap_TL_Pix(1) + Xmdn) * Lon_per_pix; 
  Estimated_Pos = [Estimated_Lat Estimated_Lon];
  
  R = 6371; % km
  % 計算實際經緯度與估測經緯度誤差距離
  
  deltaLon = Estimated_Lon - AerialP_Lon;
  deltaLat = Estimated_Lat - AerialP_Lat;
  
  % Haversine
  a = sind(deltaLat / 2)^2 + cosd(Estimated_Lat) * cosd(AerialP_Lat) * sind(deltaLon / 2)^2;
  c = 2 * atan2(sqrt(a), sqrt(1 - a));
  
  d = 1000 * R * c; % m
  
  % 判斷d大小，若太大重跑回圈(此為開發階段用)
  if d > 30 && discard_cn < 4
      k = 1;
      close(1)
      close(2)
      close(3)
      discard_cn = discard_cn +1;
      disp(['距離誤差過大: ', num2str(d, '%.2f')]);
      continue;
  end
  
  % 存進矩陣
  Bias(i) = d;
  
  % 顯示
  disp(['演算法計算之經緯度 -> 經度: ', num2str(Estimated_Lon, '%.8f'), ', 緯度: ', num2str(Estimated_Lat, '%.8f')]);
  disp(['照片實際之經緯度 -> 經度: ', num2str(AerialP_Lon, '%.8f'), ', 緯度: ', num2str(AerialP_Lat, '%.8f')]);
  disp(['估測距離差距: ', num2str(d, '%.2f'), ' 公尺']);
  
 % 將估算經緯度與實際經緯度畫出來
  figure();
  imshow(img1)
  hold on
  plot(CropMap_TL_Pix(1) + Xmdn , CropMap_TL_Pix(2) + Ymdn , 'bo', 'MarkerSize', 12, 'LineWidth', 2.5)
  
  X_shot_point = (AerialP_Lon - TL_coor(1)) / Lon_per_pix;
  Y_shot_point = (TL_coor(2) - AerialP_Lat) / Lat_per_pix;
  
  Real_points(i,:) = [X_shot_point, Y_shot_point];
  Estimated_points(i,:) = [CropMap_TL_Pix(1) + Xmdn, CropMap_TL_Pix(2) + Ymdn];
  
  hold on
  plot(X_shot_point, Y_shot_point, 'ro', 'MarkerSize', 12, 'LineWidth', 2.5)
  hold off
  legend('Calculated Position', 'Real Position', 'FontSize', 16)
 
 % 將經緯度 Pixel 更新為演算法計算出來的 ( 目前作弊用實際經緯度代替，等穩定後換成 estimated_Lon/Lat )
  AerialP_Lon = data{idx+1, 10};  % 用下一張照片位置進行拼接(也是作弊)，上機可用速度向量預測下一次迴圈位置
  AerialP_Lat = data{idx+1, 9};
  AerialP_pixX = abs((round(AerialP_Lon,8) - round(BL_coor(1),8))) / Lon_per_pix;  % 加abs目的為防止邊界誤判
  AerialP_pixY = abs((round(AerialP_Lat,8) - round(BL_coor(2),8))) / Lat_per_pix;
  
  pause(1);
  
 % 清理內存
  close all;
  clear points_img1 points_img2 f1 vpts1 f2 vpts2 tform inlierimg2Points inlierimg1Points registered2 Registered
  
 % 監控內存
  memoryInfo = memory;
  disp(['Iteration ', num2str(i), ': Memory in use: ', num2str(memoryInfo.MemUsedMATLAB / 1e6), ' MB ']);
 
  % 迴圈控制
  i = i+1;
  discard_cn = 0;  % 重設丟棄旗標
end  

% 將估測值與實際值畫在圖1
figure()
imshow(img1)
hold on
for k = 1:length(Real_points)

    plot(Real_points(k,1),  Real_points(k,2), 'ro', 'MarkerSize', 12, 'LineWidth', 2.5)
    text(Real_points(k,1), Real_points(k,2), num2str(k), 'Color', 'y', 'FontSize', 15)
   
    plot(Estimated_points(k,1),  Estimated_points(k,2), 'bo', 'MarkerSize', 12, 'LineWidth', 2.5)
    text(Estimated_points(k,1), Estimated_points(k,2), num2str(k), 'Color', 'y', 'FontSize', 15)
end
hold off

% 計算平均誤差
avg_bias = sum(Bias(1:end)) / length(Bias(1:end));
