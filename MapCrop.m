function MapCrop()
    img1 = imread('離線地圖/豐原/8KUHD_FongYuan2.png');
    
    % 圖像調整
    img1 = histeq(img1);
    img1 = imgaussfilt(img1, 2);  
    img1 = cat(3, img1(:,:,1), img1(:,:,2), img1(:,:,3));
    % h = fspecial('unsharp'); 
%     img1 = img1 + 30; 
    
    % 調整解析度
    sfactor1 = 1;
    img1 = imresize(img1, sfactor1);
    
    % 切割
    % 把矩陣大小補零，這樣在做切割時才不會有問題
    numrow = size(img1, 1); % 取img1的row數
    numcol = size(img1, 2);
    numrow = ceil(numrow / 300) * 300; % 將numrow往上取到size的倍數
    numcol = ceil(numcol / 300) * 300;
    img1(numrow,numcol,3) = 0;
    
    figure()
    imshow(img1)
    SkipStep = 300; % 每個切割後圖片塊的大小
    M_img1 = 300; % 圖片塊的長
    N_img1 = 300; % 圖片塊的寬
    n = 0; % 圖片塊的編號
    Double_img1 = im2double(img1); % 把圖片轉成雙精值
    [H,W,t] = size(Double_img1); % 得到矩陣的大小
    xStepNum = floor((W-N_img1)/SkipStep+1); % 朝負無窮方向取整，寬度方向block移動的次數
    yStepNum = floor((H-M_img1)/SkipStep+1); % 朝負無窮方向取整，長度方向block移動的次數
    
    for i = 1:xStepNum
     for j = 1:yStepNum
        n = n + 1;
        P_img1 = Double_img1((j-1)*SkipStep+1:(j-1)*SkipStep+M_img1, (i-1)*SkipStep+1:(i-1)*SkipStep+N_img1, :); % 分割圖像
        Cut_a = strcat('Cropped_Map_FU\',num2str(n),'.png'); % 儲存的圖片位置及每幅圖片塊的命名
        imwrite(double(P_img1),Cut_a);        
     end
    end

end