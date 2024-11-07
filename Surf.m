%%  Surf �t��k
%   �� Surf ���S�x���� �äǰt�Ϲ�

% 1.  ���J���չϹ�
% 2.  ��Ƕ�����O�˴�SURF�S�x�I
% 3.  ���O����SURF�y�z�l�A�Y�S�x�V�q
% 4.  �ǰt��ӹϹ����S�x�V�q
% 5.  �Q�Τǰt���G�p���̤�����transform���Y
% 6.  �z�LMSAC�N���~���S�x�I�簣
% 7.  �N�쥻���p�ϫ����ܤj�Ϥ�
% 8.  �ھڤp�Ϧb�j�ϤW����m�P�ܴ����Ytform�A�b�j�ϤW�إX�p�Ϧ�m
% 9.  �p��p�Ϫ������I��m(pixel)=>�⭫��
% 10. �إ߸g�n�ׯ��޼��ҡA�A�z�L�p�Ϧ�m(pixel)�����u��a�Ϫ��g�n�סA�ӥB�i���٭n�g�Ӥ�����S���g�n�׮y�Ъ��I
% ******************************** �ثe�i�ר��10�� *************************************
% 11. �κϤO�p(�o�ä��ǡA���S��k�A�]���ڭ̨S��GPS)�B���׭p�M����(�������A��+���Y�w�˨�)���X�����Z���A
%     �A�M�e����X����m�����v�i�o������ڦ�m
% PS. �Q�@�U�Ӧp���@�s��w��(�Y�ɩʦp��)�A�ݭn���u�ìP�a�ϩM�ũ��Ƥ~�����...

clc; clear; close all;

%% Ū���Ϥ�
rgbImage = imread('G_FullMap_FJU5.png'); % Google �a��(���� X �e�� X 3)�A���۰�Ū���ۤ��A�}�o���q���]1��Ū���@��
% rgbImage = histeq(rgbImage);
% img1 = imgaussfilt(img1, 2);

% rgbImage = cat(3, img1(:,:,1), img1(:,:,2), img1(:,:,3));


img2 = imread('�{�Ҫũ�/094.jpg'); % �L�H������
% img2 = imrotate(img2, -180);
% img2 = histeq(img2);
% img2 = imgaussfilt(img2, 2);

% �վ�img�ѪR��
sfactor1 = 1;
sfactor2 = 0.1;
rgbImage = imresize(rgbImage, sfactor1);
% rgbImage = rgbImage + 10;    % �����G��
img2 = imresize(img2, sfactor2);
% img2 = img2 - 30;    % ���C�G��
% imwrite(img2, 'Resize_pic.jpg')
% figure;
% imshow(img2);
% title('Resized Image');
% img3 = imread('twdtm_asterV2_30m.tif'); % �K�������x�W�a��(����)

% ��ܵ��G
figure('Name','Ū���m��Ϥ�')
subplot(2,1,1)
imshow(rgbImage) % ��3����size = 3�A�O�]��rgb???
title("���u�a��")
subplot(2,1,2)
imshow(img2)
title("�L�H�����e��")  

%% �Ϥ����A�ե�
roll = -2.69;   % ���۰�Ū��EXCEL��
pitch = 9.32;
yaw = 310;

% roll = 16.65;
% pitch = 9.05;
% yaw = 297.45;

roll = deg2rad(roll);
pitch = deg2rad(pitch);
yaw = deg2rad(yaw);
% totalAngle = sqrt(pitch^2 + roll^2);

% ����x�}
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

% �p��Ϲ����|�Ө��I
[h, w, ~] = size(img2);
corners = [0 0 1; w 0 1; w h 1; 0 h 1]';

% �N���I�i������ܴ�
new_corners = R * corners;

% �p��s���Ϲ����(?
min_x = min(new_corners(1, :));
max_x = max(new_corners(1, :));
min_y = min(new_corners(2, :));
max_y = max(new_corners(2, :));
min_z = min(new_corners(3, :));

% �p��s���Ϲ��ؤo
new_width = ceil(max_x - min_x);   %% ceil �O�����L���j��V���̤p���
new_height = ceil(max_y - min_y);

% �p���ܴ��x�}
tform = fitgeotrans(corners(1:2, :)', new_corners(1:2, :)', 'projective');

% ���γz���ܴ�
outputView = imref2d([new_height new_width], [min_x max_x], [min_y max_y]);
AttCorrected_img2 = imwarp(img2, tform, 'OutputView', outputView);

% �ե��᪺�Ϲ�
% figure()
% subplot(2,1,1)
% imshow(img2)
% title('���')
% subplot(2,1,2)
% imshow(AttCorrected_img2);
% title('���A�ե��᪺�Ϲ�');

%% �Ϥ�����
% New_img1 = img1; % ���Q�л\�쥻���x�}�A�ҥH�t�~�s�@�ӯx�}
% % ��x�}�j�p�ɹs�A�o�˦b�����ήɤ~���|�����D
% numrow = size(img1, 1); % ��img1��row��
% numcol = size(img1, 2);
% numrow = ceil(numrow / 100) * 100; % �Nnumrow���W����100������
% numcol = ceil(numcol / 100) * 100;
% 
% New_img1(numrow,numcol,3)=0;  % Gb1.JPG �o�ӭn
% figure;
% imshow(New_img1)
% SkipStep = 100; % �C�Ӥ��Ϋ�Ϥ������j�p
% M_img1 = 100; % �Ϥ�������
% N_img1 = 100; % �Ϥ������e
% n = 0; % �Ϥ������s��
% Double_img1 = im2double(New_img1); % ��Ϥ��ন�����
% [H,W,t] = size(Double_img1); % �o��x�}���j�p
% xStepNum = floor((W-N_img1)/SkipStep+1); % �­t�L�a��V����A�e�פ�Vblock���ʪ�����
% yStepNum = floor((H-M_img1)/SkipStep+1); % �­t�L�a��V����A���פ�Vblock���ʪ�����
% for i = 1:xStepNum
%     for j = 1:yStepNum
%         n = n + 1;
%         P_img1 = Double_img1((j-1)*SkipStep+1:(j-1)*SkipStep+M_img1, (i-1)*SkipStep+1:(i-1)*SkipStep+N_img1, :); % ���ιϹ�
%         Cut_a = strcat('Cut_GFig\',num2str(n),'.jpg'); % �x�s���Ϥ���m�ΨC�T�Ϥ������R�W
%         imwrite(double(P_img1),Cut_a);        
%     end
% end
% 
% %% �Ϥ�����(ALL)�G�ثe�O�N�������a�ϸ��J�A�i�����
% Puzzle_img1 = zeros(numrow,numcol,3); % �ɹs�᪺��Ϥ��j�p�� ()
% m = 0;
% for i = 1:xStepNum
%     for j = 1:yStepNum
%         m = m + 1;
%         % �����Ϲ�
%         Puzzle_a = strcat('Cut_GFig\',num2str(m),'.jpg'); % �x�s���Ϥ���m�ΨC�T�Ϥ������R�W
%         Puzzle_b = im2double(imread(Puzzle_a));     
%         Puzzle_img1((j-1)*SkipStep+1:(j-1)*SkipStep+M_img1, (i-1)*SkipStep+1:(i-1)*SkipStep+N_img1, :) = Puzzle_b;
%     end
% end
% w = strcat('Puzzle.jpg'); % �N�X���Ϥ��x�s�ܷQ�n����m�éR�W
% figure('Name','�������᪺�Ϥ�[�`��]')
% imshow(double(Puzzle_img1));
% imwrite(double(Puzzle_img1),w);
% 
% %% �Ϥ�����(3*3)�G�ھڸ��J���ϸ����X�P��8����
% % % �]���{�b�u�O�Ϥ��A�ҥH�u�ઽ�����s�����m�A�������ӬO�ھڸg�n�ץh��ڭ̭n�����i�Ϥ��A�o��������A�g
% % promp1 = '�п�J�Ϥ��s��(1 ~ ';
% % promp2 = num2str(xyStepNum);
% % promp3 = ' �䤤�@��):';
% % promp = append(promp1,promp2,promp3);
% % % ��ܭn���a�Ͻs��
% % p = input(promp);
% 
% % ��ܭn���a�Ͻs��
% xyStepNum = xStepNum*yStepNum;
% p = input(append('�п�J�Ϥ��s��(1 ~ ',num2str(xyStepNum),' �䤤�@��):'));
% 
% % �A�N���i�Ϥ��ΥL�P�D��8�i�@��Ū���X�ӡA�Ҷq�i��|�����ɩάO�p��쪺�Ϥ��s�����b�d�򤺪����p�A�Y�o�ͳo�ز{�H�h�ιs���N
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
% % �קK�W�U��ɪ��~�P
% for  count = 1:1:xStepNum % �Ĥ@�ƩM�̫�@�Ƥ��O��11�ӼƦr�n�P�_  
%     % �Y p ���Ȧ��Ĥ@��
%     if p == (count-1)*yStepNum+1
%         % �N�� MapNum1�BMapNum4 �M MapNum7 �ιs�x�}������
%         MapNum1 = zeros([M_img1,N_img1,3]);  
%         MapNum4 = zeros([M_img1,N_img1,3]);   
%         MapNum7 = zeros([M_img1,N_img1,3]);
%     end
%     % �Y p ���Ȧ��̥���
%     if p == (count)*yStepNum
%         % �N�� MapNum3�BMapNum6 �M MapNum9 �ιs�x�}������
%         MapNum3 = zeros([M_img1,N_img1,3]);  
%         MapNum6 = zeros([M_img1,N_img1,3]);   
%         MapNum9 = zeros([M_img1,N_img1,3]);
%     end
% end
% 
% % �A��E�i�Ϥ��Ӷ��ǫ��^�h
% Puzzle_c =[MapNum1 MapNum4 MapNum7; MapNum2 MapNum5 MapNum8; MapNum3 MapNum6 MapNum9]; 
% 
% % ����n���Ϥ���ܥX��
% figure('Name','�������᪺�Ϥ�[3*3]');
% imshow(double(Puzzle_c));
% 
% % �ç�Ϥ��x�s�_��
% PictureName = strcat('Puzzle99.jpg'); % �N�X���Ϥ��x�s�ܷQ�n����m�éR�W
% imwrite(double(Puzzle_c),PictureName);

%% �Ϥ���Ƕ�
% img1 = imread('Puzzle9.jpg');  % �N9�c��a�ϭ��s���Nimg1�A���ھڵL�H����m��ܤE�c�椤���I(�L�H���ثe��mŪ�����P�D�E�c���) -> �ݽT�{
Gray_img1 = rgb2gray(rgbImage);
% Gray_img1 = histeq(Gray_img1);

Gray_img2 = rgb2gray(img2);
% Gray_img2 = histeq(Gray_img2);
% Gray_img2 = imgaussfilt(Gray_img2, 2);

% Gray_img3 = rgb2gray(AttCorrected_img2);
imageSize = size(rgbImage);
[M,N] = size(Gray_img2);
% ��ܵ��G
figure('Name','�ন�Ƕ��Ϥ�')
subplot(2,1,1)
imshow(Gray_img1) 
title("Google�a��")
subplot(2,1,2)
imshow(Gray_img2)
title("�L�H�����a��")

%% �Ϥ��G�Ȥ�(SURF�u�ݭn�ǫ׹ϡA�ҥH�o�����Τ���)
% Bi_mg1 = im2bw(Gray_img1,0.7); 
% Bi_mg2 = im2bw(Gray_img2,0.7);

%% ��Ϲ����S�x�I
tic
points1 = detectSURFFeatures(Gray_img1);
points2 = detectSURFFeatures(Gray_img2);  % �i�H�յ۽� Threshold
% points3 = detectSURFFeatures(Gray_img3); % �g���A�ե����ũ�Ӥ�

% ��ܵ��G
figure('Name','�索�k�Ϲ����O�i��S�x�I����')
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

% figure('Name', '��i�櫺�A�ե��᪺�϶i��S�x�I����')
% imshow(Gray_img3); hold on;
% [Asize, l] = size(points3.Location);
% plot(points3(1:2:Asize, :), 'showOrientation', true);

%% �v���ե��G�]���۾��|���Ө���(�w�˨�+�������A)�A�ҥH�v���|�ܧ�(�������)�A�]���ݭn��v���Ԧ^������v�����A
% �g�L���յo�{SURF�t��k�������X��ե����\��A�ҥH�ثe���������μg�o����
% Theta = 0:180; % 0:179;
% R = radon(Gray_img2,Theta);
% % �D�X�Ϲ������I����ɪ��Z��
% L = round(sqrt((M/2)^2+(N/2)^2));
% [C,angle] = max(R(L,:)); % �������ӧڭ̦ۤv���|����n
% % angle ���Ϲ��ɱר���
% angle = angle - 1;
% % �N�Ϥ����ե�
% A = imrotate(img2,angle,'nearest');
% 
% % ��ܮե��᪺�Ϥ�
% figure;imshow(A);

%% �p��y�z�V�q
[f1, vpts1] = extractFeatures(Gray_img1, points1);
[f2, vpts2] = extractFeatures(Gray_img2, points2);
% [f3, vpts3] = extractFeatures(Gray_img3, points3);

%% �i��ǰt�F�����S�x�I��m�ASURF �S�x�V�q���������W�ƳB�z
% �g���A�ե��e���ũ��
% indexPairs1 = matchFeatures(f1, f2) ;
indexPairs1 = matchFeatures(f1, f2, 'Prenormalized', true) ; 
% indexPairs1 = matchFeatures(f1, f2, 'Method', 'Approximate') ;
% indexPairs1 = matchFeatures(f1, f2, 'Unique', true) ;
matched_pts1 = vpts1(indexPairs1(:, 1));
matched_pts2 = vpts2(indexPairs1(:, 2));

% �g���A�ե��᪺�ũ��
% % indexPairs2 = matchFeatures(f1, f2,'MatchThreshold',2) ;
% indexPairs2 = matchFeatures(f1, f3, 'Prenormalized', true) ; 
% % indexPairs2 = matchFeatures(f1, f2, 'Method', 'Approximate') ;
% matched_pts3 = vpts1(indexPairs2(:, 1));
% matched_pts4 = vpts3(indexPairs2(:, 2));

%% ��MSAC�t��k�h�����~���ǰt�I
% �q�L�S�x�I�ǰt�ٱo��F�ĤG�T�Ϫ��ܴ��x�}tform�A�ĤG�T�ϭn�g�L�ܴ��x�}�ܦ��M�Ĥ@�T�Ϫ��y�Ф@�P
[tform, inlierimg2Points, inlierimg1Points] = estimateGeometricTransform(matched_pts2, matched_pts1, 'similarity'); %�g�v�ܴ��Atform�M�g�I��1���I���I��2���I
% �Ө�ƨϥ��H���˥��@�P�ʡ]RANSAC�ARandom Sample Consensus�^�t��k������MSAC�t��k��{�A�h���~�ǰt�I
% ��^���X��M�g�x�}�M�g�Ĥ@�ѼƤ��I��ĤG�ѼƤ��I
toc
% ��ܹ�ǰt���G�A�i�H�ݨ����٦��@�ǲ��`�ȡA���g�LMSAC�ᦳ����ﵽ
  figure('Name','�索�k�Ϲ����O�i��S�x�I����')
  subplot(2,1,1)
  showMatchedFeatures(Gray_img1,Gray_img2,matched_pts1,matched_pts2,'montage');
  legend('matched points 1','matched points 2'); 
  title("�S�x�Ȥǰt���G")
  subplot(2,1,2)
  showMatchedFeatures(Gray_img1,Gray_img2,inlierimg1Points,inlierimg2Points,'montage');
  legend('matched points 1','matched points 2');
  title("�h�����~���ǰt�I�᪺�S�x�Ȥǰt���G")
  
% �P�˪��ƱN�a�ϻP���A�ե��L��ϦA���@��
%   [tform_attC, inlierimg2Points, inlierimg1Points] = estimateGeometricTransform(matched_pts4, matched_pts3, 'similarity', 'MaxDistance', 4); %�g�v�ܴ��Atform�M�g�I��1���I���I��2���I
%   figure('Name', '�g���A�ե���索�k�ϯS�x�I����')
%   subplot(2,1,1)
%   showMatchedFeatures(Gray_img1,Gray_img3,matched_pts3,matched_pts4,'montage');
%   title("�S�x�ǰt���G")
%   subplot(2,1,2)
%   showMatchedFeatures(Gray_img1,Gray_img3,inlierimg1Points,inlierimg2Points,'montage');
%   title("�h�����~�I�᪺�S�x�ǰt���G")
    
%% �i��Ϲ��X�֡Atform�O�ܴ��x�}�A�H�Ĥ@�T�Ϲ�����Ǯy�СA�ĤG�T�ϭn�i���ܴ��P�����
% Rfixed ���Ĥ@�T�Ϫ��@�ɤG���y��
  Rfixed = imref2d(size(rgbImage));
  [registered2, Rregistered] = imwarp(img2, tform);
  boxPolygon = [1, 1;... % ���W
    size(img2, 2), 1; ... % �k�W
    size(img2, 2), size(img2, 1); ... % �k�U
    1, size(img2, 1); ... % ���U
    1, 1]; % ���ƥ��U�A�~�i�o��@�ӳ��϶����h���
% �N�h����ܴ���ؼйϤ��W�A�ܴ������G��ܤF���骺��m
  newBoxPolygon = transformPointsForward(tform, boxPolygon);
% ��ܳQ�˴��쪺����
  figure()
  imshowpair(rgbImage,Rfixed,registered2,Rregistered,'blend');
%%% �p������Ϥ������ߦ�m(pixel)�G����ھڳo�ӥh����������ޭ�
% �n��쥻�����٬O�I�����ءH���I�����اڭn�Q�@�U����
  hold on;
% [xlim, ylim] = outputLimits(tform, [1 imageSize(1)], [1 imageSize(2)]);
% �Ĥ@�غ⭫�ߪ��覡�G��N����
  Xmdn = (newBoxPolygon(1, 1)+newBoxPolygon(2, 1)+newBoxPolygon(3, 1)+newBoxPolygon(4, 1))/4;
  Ymdn = (newBoxPolygon(1, 2)+newBoxPolygon(2, 2)+newBoxPolygon(3, 2)+newBoxPolygon(4, 2))/4;
  line(newBoxPolygon(:, 1), newBoxPolygon(:, 2), 'Color', 'y');
  hold on;
  InterestPoint1 = scatter(Xmdn, Ymdn, '*b'); % �Mmatlab�⪺���ߦ��ǷL�~�t
% �ĤG�غ⭫�ߪ��覡�G�� matlab ���
  x1 = newBoxPolygon(1:4, 1)';
  y1 = newBoxPolygon(1:4, 2)';
  polyin = polyshape(x1,y1); 
  [CenterX,CenterY] = centroid(polyin);% matlab�۱a�p�⭫�ߪ����
  InterestPoint2 = scatter(CenterX, CenterY, '*g'); 
  title('�Ϲ��t��');

%% ���Ϋ��A�ե��L�᪺�ϦA���@��
%   [registered_attC, Rregistered_attC] = imwarp(AttCorrected_img2, tform_attC);
%   boxPolygon_attC = [1, 1;... % ���W
%     size(AttCorrected_img2, 2), 1; ... % �k�W
%     size(AttCorrected_img2, 2), size(AttCorrected_img2, 1); ... % �k�U
%     1, size(AttCorrected_img2, 1); ... % ���U
%     1, 1]; % ���ƥ��U�A�~�i�o��@�ӳ��϶����h���
%   newBoxPolygon_attC = transformPointsForward(tform_attC, boxPolygon_attC);
%   figure()
%   imshowpair(img1,Rfixed,registered_attC,Rregistered_attC,'blend');
% % �n��쥻�����٬O�I�����ءH���I�����اڭn�Q�@�U����
%   hold on;
% % [xlim, ylim] = outputLimits(tform, [1 imageSize(1)], [1 imageSize(2)]);
% % �Ĥ@�غ⭫�ߪ��覡�G��N����
%   Xmdn = (newBoxPolygon_attC(1, 1)+newBoxPolygon_attC(2, 1)+newBoxPolygon_attC(3, 1)+newBoxPolygon_attC(4, 1))/4;
%   Ymdn = (newBoxPolygon_attC(1, 2)+newBoxPolygon_attC(2, 2)+newBoxPolygon_attC(3, 2)+newBoxPolygon_attC(4, 2))/4;
%   line(newBoxPolygon_attC(:, 1), newBoxPolygon_attC(:, 2), 'Color', 'y');
%   hold on;
%   InterestPoint1 = scatter(Xmdn, Ymdn, '*b'); % �Mmatlab�⪺���ߦ��ǷL�~�t
% % �ĤG�غ⭫�ߪ��覡�G�� matlab ���
%   x1 = newBoxPolygon_attC(1:4, 1)';
%   y1 = newBoxPolygon_attC(1:4, 2)';
%   polyin = polyshape(x1,y1); 
%   [CenterX,CenterY] = centroid(polyin);% matlab�۱a�p�⭫�ߪ����
%   InterestPoint2 = scatter(CenterX, CenterY, '*g'); 
%   title('�Ϲ��t��');
 
%%% �p������Ϥ������ߦ�m(pixel)�G����ھڳo�ӥh����������ޭ�
% �n��쥻�����٬O�I�����ءH���I�����اڭn�Q�@�U����
%   hold on;
% % [xlim, ylim] = outputLimits(tform, [1 imageSize(1)], [1 imageSize(2)]);
% % �Ĥ@�غ⭫�ߪ��覡�G��N����
%   Xmdn = (newBoxPolygon_attC(1, 1)+newBoxPolygon_attC(2, 1)+newBoxPolygon_attC(3, 1)+newBoxPolygon_attC(4, 1))/4;
%   Ymdn = (newBoxPolygon_attC(1, 2)+newBoxPolygon_attC(2, 2)+newBoxPolygon_attC(3, 2)+newBoxPolygon_attC(4, 2))/4;
%   line(newBoxPolygon_attC(:, 1), newBoxPolygon_attC(:, 2), 'Color', 'y');
%   hold on;
%   InterestPoint1 = scatter(Xmdn, Ymdn, '*b'); % �Mmatlab�⪺���ߦ��ǷL�~�t
% % �ĤG�غ⭫�ߪ��覡�G�� matlab ���
%   x1 = newBoxPolygon_attC(1:4, 1)';
%   y1 = newBoxPolygon_attC(1:4, 2)';
%   polyin = polyshape(x1,y1); 
%   [CenterX,CenterY] = centroid(polyin);% matlab�۱a�p�⭫�ߪ����
%   InterestPoint2 = scatter(CenterX, CenterY, '*g'); 
%   title('�Ϲ��t��');
  
 %% �H�Ӥ������I�y�Эp��L�H����ڸg�n��
 % ���u�a�ϥ|�Ө����g�n��
%   Pos_Map_UL = [24.265632, 120.815443];
%   Pos_Map_LL = [24.263959, 120.815445];
%   Pos_Map_UR = [24.265644, 120.817518];
%   Pos_Map_LR = [24.264007, 120.817536];
%  % �p�� X/Y �b�C�ӹ����O�h�ָg�n��(�g�n�ץ�x10e6)
%   Lon_per_pixel = (Pos_Map_UR(2) * 10^6 - Pos_Map_UL(2) * 10^6) / size(img1, 2);
%   Lat_per_pixel = (Pos_Map_UL(1) * 10^6 - Pos_Map_LL(1) * 10^6) / size(img1, 1);
%  % �p����Ӥ��g�n��
%   Estimated_Lat = Pos_Map_UL(1) * 10^6 - Ymdn*Lat_per_pixel;
%   Estimated_Lon = Pos_Map_UL(2) * 10^6 + Xmdn*Lon_per_pixel;
%   Estimated_Pos = [Estimated_Lat Estimated_Lon] / 10^6;
%  % �N����g�n�׻P��ڸg�n�׵e�X��(���i�ण�ǡA�����a�ϥ|�Ө��O�Τ��I��)
%   hold on
%   plot(Xmdn, Ymdn, 'ro', 'MarkerSize', 12)
%   
%   X_shot_point = (120.8169532 - Pos_Map_UL(2))*10^6 / Lon_per_pixel;
%   Y_shot_point = (Pos_Map_UL(1) - 24.2643635)*10^6 / Lat_per_pixel;
%   hold on
%   plot(X_shot_point, Y_shot_point, 'bo', 'MarkerSize', 12)
%   hold off
  
 % ��� [24.2643635000000,120.816953200000]
 % ���� [24.264338,120.81686]
 % �~�t����10����
 
%% �H�Ϥ������I�y�Эp���ڭ�����m(�o�̸�T������A�ٵL�k����)
% % �N������m������g�n��(�ڲ{�b�٨S��index�A�ҥH����쥻����)
% PictureLon = CenterX;
% PictureLat = CenterY;
% % ���Y�w�˨��׬O20�סA��ڥ���
% CameraAngle = 20;
% % ����������
% Xsens_Pitch = Data(:,12); % ������
% %  Madgwick_Pitch = Data(:,15); % ������
% CPAngle = 90-(CameraAngle + Xsens_Pitch);
% % �������סAm(NED)�A���F�S��LLA���Шt������!!!
% Air_Altitude = Data(:,19); 
% LLA_Altitude = �L�k���olla��������;% ���׭n��LLA���Шt�����M�|��@@!!!
% 
% % ��k�@�G����k
% �p��׶Z�G�o�̥Ϊ���k����]�O�]�����ݭn�y���ഫ�A�o�˥i�H���C�p��ɶ��A�S�[�W�Z���ܵu�A�ҥH�i�H������έp��
%  ARC = 6371.393*1000; % �����b�|(m)
%  TrueLon = PictureLon + LLA_Altitude/(ARC*cos(PictureLon)*2*pi/360);
%  TrueLat = PictureLat +(LLA_Altitude/cos(CPAngle))/ (ARC *2*pi/360);
% 
% % ��k�G�G�y���ഫ
% % LLA �y�Шt �� NED �y�Шt [xNorth,yEast,zDown] = geodetic2ned(lat,lon,h,lat0,lon0,h0,spheroid) 
% % �_����m�G
%   LonS_GPS = Data(1,8); % GPS LLA �y�ФU�� �g�� ��m
%   LatS_GPS = Data(1,9); % GPS LLA �y�ФU�� �n�� ��m
%   HS_GPS = Data(1,10); % GPS LLA �y�ФU�� ���� ��m
% % �ϤO�p
%   MagX = Data(:,8); % X �b
%   MagY = Data(:,9); % Y �b 
% % LLA �y�Шt �� NED �y�Шt [xNorth,yEast,zDown] = geodetic2ned(lat,lon,h,lat0,lon0,h0,spheroid)�A�{�b�٨S����� 
%  [GPS_N_S,GPS_E_S,GPS_D_S] = geodetic2ned(PictureLon,PictureLat,�S���צn�{�e,LonS_GPS,LatS_GPS,HS_GPS,wgs84Ellipsoid);
% %  True_N_S = (MagX/MagY) *True_E_S;
% %  (True_N_S - GPS_N_S)^2 + (True_E_S - GPS_E_S)^2 = (Air_Altitude/cos(CPAngle))^2;
% % �@���G����{��(�p��Z��)
%  syms True_N_S True_E_S;
%  f = 5*True_N_S+(-2*True_N_S*GPS_N_S-4*True_N_S*GPS_E_S)-(Air_Altitude/cos(CPAngle))^2+GPS_N_S^2+GPS_E_S^2;
%  % �D�Ѥ�{��
%  True_N_S =solve(f==0);
%  True_E_S = 2*True_N_S;
%  % NED �y�Шt �� LLA �y�Шt
%  [True_LAT_S, True_LON_S, True_H] = ned2geodetic(True_N_S,True_E_S,�S�S����,LonS_GPS,LatS_GPS,HS_GPS,wgs84Ellipsoid); 
 

%%  *************************************************** �H�U�O�ۤv���p��s�A�M�ץΤ��� ***************************************************
%% �T�w��T�Ϲ�ڤj�p�Φ�m 
% tic
% %��X�y�нd��
% [xlim, ylim] = outputLimits(tform, [1 imageSize(1)], [1 imageSize(2)]);
% 
% % ����X�Ŷ�����̤j�̤p��
% xMin = min([1; xlim(:)]);
% xMax = max([imageSize(2); xlim(:)]); 
% yMin = min([1; ylim(:)]);
% yMax = max([imageSize(1); ylim(:)]);
%  
% % �����Ϫ��e��
% width  = round(xMax - xMin);
% height = round(yMax - yMin); 
%  
% %�Ы�2D�Ŷ��ѦҪ���w�q�����Ϥؤo
% xLimits = [xMin xMax];
% yLimits = [yMin yMax];
% panoramaView = imref2d([height width ], xLimits, yLimits);
%  
% % �ܴ��Ϥ��������
% unwarpedImage = imwarp(img1,projective2d, 'OutputView', panoramaView);
% warpedImage = imwarp(img2, tform, 'OutputView', panoramaView);
% 
% %% ���J���X�ĦX�G������i�Ϥ�������t�A���䧹���ĤJ
% % �ҿ׺��J���X�N�O�N��T�ϭ��X���ϰ���ӶZ����T�Ϫ��Z�����Ӥ@�w���v�����s���t���X�����ϵe���T����v��
% % ��p�̤������N�O0.5 0.5����ҡC�U�@�B�N�O���L�̭��|�ϰ�A�]�N�O�ۦP���Ұ�
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
% maskA = (warpedImage(:,:,1)>0 |warpedImage(:,:,2)>0 | warpedImage(:,:,3)>0);%�ܴ��Ϲ�����
% mask1 = (newImage(:,:,1)>0 | newImage(:,:,2)>0 | newImage(:,:,3)>0);%�D�ܴ��Ϲ�����
% mask1 = and(maskA, mask1);%���|�ϱ���
% % figure,imshow(mask1)
% [row,col] = find(mask1==1);
% left = min(col);
% right = max(col);%��o���|�ϥ��k�d��
% up=min(row);
% down=max(row);
% mask = ones(size(mask1));
% % figure()
% % imshow(mask)
% %mask(:,left:right) = repmat(linspace(0,1,right-left+1),size(mask,1),1);%�ƻs���Q�x�}
% mask(up:down,:) = repmat(linspace(1,0,down-up+1)',1,size(mask,2));%�ƻs���Q�x�}
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
% mask(up:down,:) = repmat(linspace(0,1,down-up+1)',1,size(mask,2));%�ƻs���Q�x�}
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
% title('���J���X�ĦX');
% toc


