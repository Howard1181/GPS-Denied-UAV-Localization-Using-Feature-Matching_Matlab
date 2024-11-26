#pragma once
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <filesystem>
#include <opencv2/opencv.hpp>

#define M_PI 3.1415926

using namespace std;
namespace fs = std::filesystem;

class ImgLocate
{


	// �a�ϥ|���g�n��


	// �Ӥ�Ū����T
	string folderPath_camera = "E:/�L�H���v���w��/CY�v���w���@/�׭�ũ�_2";
	string PhoInfo_excelFile = "../../FongYuanLog2.csv";

	// �t��k�Ѽ�
	int xStepNum = 0;
	int yStepNum = 0;
	int xyStepNum = 0;
	int discard_count = 0; // �˱󦸼�


	// �B��禡
	cv::Mat im2double(const cv::Mat& img);
	cv::Mat histMatch(cv::Mat&, cv::Mat&);


public:
	vector<string> getCamImgFiles();
	vector<vector<string>> getCamImgLog();
	string CroppedMapPath = "E:/�L�H���v���w��/CY�v���w���@/Cropped_Map_FU"; // ������|�]�i�H
	string folderPath_map = "E:/VS/ImgLocation/Map/8KUHD_FongYuan2.png";
	cv::Mat  puzzleMap(int picture_num);
	void setResolution(cv::Mat&, double);
	void setAngle(cv::Mat&, double);
	void grayImg(cv::Mat&);
	void adjustBright_Contrast(cv::Mat&, cv::Mat&, int);
	void PresetPhoto(cv::Mat&, cv::Mat&, double, double);
	cv::Point2d locatingAlgorithm(cv::Mat&, cv::Mat&, double);
	cv::Point2d PositionCalculation(int, int, int, double, double, cv::Point2d);


	// �ŧi�a��/�ũ�򥻸�T(�p��ʸ�?)(�i�Ϋغc�l?�άO�g�@�禡)
	double Map_Width = 4103; // �׭�
	double Map_Height = 3894;
	double CropSize = 300;
	int numrow = ceil(Map_Width / CropSize) * CropSize; // �`�NMap_Width/Height, CropSize�n�]double���G�~�|��
	int numcol = ceil(Map_Height / CropSize) * CropSize;

	double TL_coor[2] = { 120.7157564, 24.2630274 };
	double TR_coor[2] = { 120.7223255, 24.2630537 };
	double BL_coor[2] = { 120.7158070, 24.2573680 };
	double BR_coor[2] = { 120.7223298, 24.2573723 };
	double Lon_per_pix = (TR_coor[0] - TL_coor[0]) / Map_Width;
	double Lat_per_pix = (TL_coor[1] - BL_coor[1]) / Map_Height;

	// �t��k�Ѽ�
	bool k = false; // �t��k���]�X��

};

// �N�Ӥ��ഫ��double���禡
cv::Mat ImgLocate::im2double(const cv::Mat& img)
{
	cv::Mat img_double;
	img.convertTo(img_double, CV_64F, 1.0 / 255.0);
	return img_double;
}

// ���o�ũ�Ӥ����|string
vector<string> ImgLocate::getCamImgFiles()
{
	vector<string> imageFiles_camera;
	if (fs::exists(folderPath_camera) && fs::is_directory(folderPath_camera))
	{
		for (const auto& entry : fs::directory_iterator(folderPath_camera))
		{
			imageFiles_camera.push_back(entry.path().string());
		}
	}
	else
	{
		cerr << "��Ƨ����s�b�Τ��O�@�Ӧ��Ī���Ƨ�" << folderPath_camera << endl;
	}
	return imageFiles_camera;
}

// ���o�ũ�Ӥ�log��
vector<vector<string>> ImgLocate::getCamImgLog()
{
	ifstream file(PhoInfo_excelFile);
	string line;
	vector<vector<string>> data;
	if (!file.is_open())
	{
		cerr << "�L�k���}�ɮ�" << PhoInfo_excelFile << endl;
		return data;
	}

	// �v��Ū��csv
	while (getline(file, line))
	{
		stringstream ss(line);
		string cell;
		vector<string> row; // �

		while (getline(ss, cell, ','))
		{
			row.push_back(cell);
		}

		data.push_back(row);
	}
	file.close();
	return data;
}

// ������e�y�зӤ�
cv::Mat ImgLocate::puzzleMap(int picture_num)
{
	cv::Mat MapNum1, MapNum2, MapNum3, MapNum4, MapNum5, MapNum6, MapNum7, MapNum8, MapNum9 = cv::Mat::zeros(CropSize, CropSize, CV_8UC3);
	cv::Mat Puzzle_map33(CropSize * 3, CropSize * 3, CV_8UC3);
	xStepNum = floor((numrow - CropSize) / CropSize + 1);
	yStepNum = floor((numcol - CropSize) / CropSize + 1);
	xyStepNum = xStepNum * yStepNum;

	if ((picture_num - yStepNum - 1) < 1 || (picture_num - yStepNum - 1) > xyStepNum)
	{
		MapNum1.setTo(cv::Scalar(0, 0, 0)); // �¦�
	}
	else {
		string full_path = CroppedMapPath + "/" + to_string(picture_num - yStepNum - 1) + ".png";
		MapNum1 = cv::imread(full_path);
	}

	if ((picture_num - yStepNum) < 1 || (picture_num - yStepNum) > xyStepNum)
	{
		MapNum2.setTo(cv::Scalar(0, 0, 0)); // �¦�
	}
	else {
		string full_path = CroppedMapPath + "/" + to_string(picture_num - yStepNum) + ".png";
		MapNum2 = cv::imread(full_path);
	}

	if ((picture_num - yStepNum + 1) < 1 || (picture_num - yStepNum + 1) > xyStepNum)
	{
		MapNum3.setTo(cv::Scalar(0, 0, 0)); // �¦�
	}
	else {
		string full_path = CroppedMapPath + "/" + to_string(picture_num - yStepNum + 1) + ".png";
		MapNum3 = cv::imread(full_path);
	}

	if ((picture_num - 1) < 1 || (picture_num - 1) > xyStepNum)
	{
		MapNum4.setTo(cv::Scalar(0, 0, 0)); // �¦�
	}
	else {
		string full_path = CroppedMapPath + "/" + to_string(picture_num - 1) + ".png";
		MapNum4 = cv::imread(full_path);
	}

	if ((picture_num) < 1 || (picture_num) > xyStepNum)
	{
		MapNum5.setTo(cv::Scalar(0, 0, 0)); // �¦�
	}
	else {
		string full_path = CroppedMapPath + "/" + to_string(picture_num) + ".png";
		MapNum5 = cv::imread(full_path);
	}

	if ((picture_num + 1) < 1 || (picture_num + 1) > xyStepNum)
	{
		MapNum6.setTo(cv::Scalar(0, 0, 0)); // �¦�
	}
	else {
		string full_path = CroppedMapPath + "/" + to_string(picture_num + 1) + ".png";
		MapNum6 = cv::imread(full_path);
	}

	if ((picture_num + yStepNum - 1) < 1 || (picture_num + yStepNum - 1) > xyStepNum)
	{
		MapNum7.setTo(cv::Scalar(0, 0, 0)); // �¦�
	}
	else {
		string full_path = CroppedMapPath + "/" + to_string(picture_num + yStepNum - 1) + ".png";
		MapNum7 = cv::imread(full_path);
	}

	if ((picture_num + yStepNum) < 1 || (picture_num + yStepNum) > xyStepNum)
	{
		MapNum8.setTo(cv::Scalar(0, 0, 0)); // �¦�
	}
	else {
		string full_path = CroppedMapPath + "/" + to_string(picture_num + yStepNum) + ".png";
		MapNum8 = cv::imread(full_path);
	}

	if ((picture_num + yStepNum + 1) < 1 || (picture_num + yStepNum + 1) > xyStepNum)
	{
		MapNum9.setTo(cv::Scalar(0, 0, 0)); // �¦�
	}
	else {
		string full_path = CroppedMapPath + "/" + to_string(picture_num + yStepNum + 1) + ".png";
		MapNum9 = cv::imread(full_path);
	}

	// �קK�W�U��ɻ~�P
	for (int count = 0; count < xStepNum; count++)
	{
		if (picture_num == (count - 1) * yStepNum + 1)
		{
			MapNum1.setTo(cv::Scalar(0, 0, 0));
			MapNum4.setTo(cv::Scalar(0, 0, 0));
			MapNum7.setTo(cv::Scalar(0, 0, 0));
		}
		if (picture_num == count * xyStepNum)
		{
			MapNum3.setTo(cv::Scalar(0, 0, 0));
			MapNum6.setTo(cv::Scalar(0, 0, 0));
			MapNum9.setTo(cv::Scalar(0, 0, 0));
		}
	}

	// ������ 3*3 ���p�a��
	cv::Mat row1, row2, row3;
	cv::hconcat(vector<cv::Mat>{MapNum1, MapNum4, MapNum7}, row1); // �٥i�H�o�˥�vector?
	cv::hconcat(vector<cv::Mat>{MapNum2, MapNum5, MapNum8}, row2);
	cv::hconcat(vector<cv::Mat>{MapNum3, MapNum6, MapNum9}, row3);

	cv::vconcat(vector<cv::Mat>{row1, row2, row3}, Puzzle_map33);

	cv::imshow("FinalImg", Puzzle_map33);
	cv::waitKey(800);
	cv::destroyWindow("FinalImg");

	// �N�����n���Ӥ��g�J��Ƨ�
	string puzzlePath = "E:/�L�H���v���w��/CY�v���w���@/�p�a�ϫ���/Puzzle33_cpp.png";
	cv::imwrite(puzzlePath, Puzzle_map33);
	return Puzzle_map33;
}

// �վ�Ӥ��ѪR�ר��
void ImgLocate::setResolution(cv::Mat& img2, double sFactor)
{
	cv::resize(img2, img2, cv::Size(), sFactor, sFactor);
}

// �վ�Ӥ�����
void ImgLocate::setAngle(cv::Mat& img2, double angle)
{
	cv::Point2f center(img2.cols / 2.0F, img2.rows / 2.0F); // �p��Ϲ�����
	cv::Mat rotateMatrix = cv::getRotationMatrix2D(center, angle, 1.0); // �Ыر���x�}
	cv::Rect boundingbox = cv::RotatedRect(center, img2.size(), angle).boundingRect(); // �p������Ϲ��j�p
	rotateMatrix.at<double>(0, 2) += boundingbox.width / 2.0 - center.x; // �����H�ŦX���
	rotateMatrix.at<double>(1, 2) += boundingbox.height / 2.0 - center.y;

	cv::warpAffine(img2, img2, rotateMatrix, boundingbox.size());
	//cv::namedWindow("rotated_img", cv::WINDOW_NORMAL);
	//cv::imshow("rotated_img", rotated_img);
	//cv::waitKey(0);
}

// �Ƕ���
void ImgLocate::grayImg(cv::Mat& img)
{
	cv::cvtColor(img, img, cv::COLOR_BGR2GRAY);
	//cv::namedWindow("rotated_img", cv::WINDOW_NORMAL);
	//cv::imshow("rotated_img", grayImg);
	//cv::waitKey(0);
}

// �P�_�ýվ�G�׹��
void ImgLocate::adjustBright_Contrast(cv::Mat& img1, cv::Mat& img2, int dicard_cn)
{
	cv::Mat adjusted_img2;
	cv::Mat temp_img1, temp_img2;
	cv::Scalar mean1, mean2, std1, std2; // �|��double���o���O�A���O��BGR�P�z����

	temp_img1 = img1;
	temp_img2 = img2;
	temp_img1.convertTo(temp_img1, CV_32F); // �N�Ϲ��ন�B�I�ƫ��O
	temp_img2.convertTo(temp_img2, CV_32F);
	cv::normalize(temp_img1, temp_img1, 0, 255, cv::NORM_MINMAX); // ���W�ƨ�[0 255]����
	cv::normalize(temp_img2, temp_img2, 0, 255, cv::NORM_MINMAX);

	cv::meanStdDev(temp_img1, mean1, std1);
	cv::meanStdDev(temp_img2, mean2, std2);

	cout << "������:" << mean1 << "," << mean2 << endl;
	cout << "�зǮt:" << std1 << "," << std2 << endl;

	// �̾ڥ����ȼзǮt�P�_�O�_�n���ק��@��
	if (abs(mean1[0] - mean2[0]) >= 10 && abs(mean1[0] - mean2[0] <= 25) || abs(std1[0] - std2[0]) >= 10)
	{
		if (discard_count <= 1)
		{
			cout << "������ק��@��" << endl;
			adjusted_img2 = histMatch(img1, img2);
		}
		else if (discard_count > 1)
		{
			cout << "���~�F2���A��������ק��@��" << endl;
			adjusted_img2 = img2;
		}
		temp_img2 = adjusted_img2;
		cv::normalize(temp_img2, temp_img2, 0, 255, cv::NORM_MINMAX);
		cv::meanStdDev(temp_img2, mean2, std2);
		cout << "�վ�ᥭ����:" << mean1 << "," << mean2 << endl;
		cout << "�վ��зǮt:" << std1 << "," << std2 << endl;
	}
	else
	{
		if (discard_count <= 1)
		{
			cout << "��i�Ӥ��t���p�ιL�j�A��������ק��@��" << endl;
			adjusted_img2 = img2;
		}
		else if (discard_count > 1)
		{
			adjusted_img2 = histMatch(img1, img2);
			cout << "���~�F2���A������ק��@��" << endl;
			temp_img2 = adjusted_img2;
			cv::normalize(temp_img2, temp_img2, 0, 255, cv::NORM_MINMAX);
			cv::meanStdDev(temp_img2, mean2, std2);
			cout << "�վ�ᥭ����:" << mean1 << "," << mean2 << endl;
			cout << "�վ��зǮt:" << std1 << "," << std2 << endl;
		}

	}
}

// ���ק��@��
cv::Mat ImgLocate::histMatch(cv::Mat& map, cv::Mat& photo)
{
	if (photo.empty() || map.empty())
	{
		cerr << "Error: one or both img empty" << endl;
	}

	cv::Mat mapHist, photoHist;
	int histSize = 256;
	float range[] = { 0, 256 };
	const float* histRange = { range };

	// �p�⪽���
	// calcHist(��J�Ϲ�, ��J�Ϲ��ƶq, �ϥΪ��q�D����, ���ϥα��X(�i�w�ﳡ���ϰ�p��Ȥ��), ��X�����, ����Ϻ���, �C�Ӻ��ת����Ƚd��)
	cv::calcHist(&photo, 1, 0, cv::Mat(), photoHist, 1, &histSize, &histRange);
	cv::calcHist(&map, 1, 0, cv::Mat(), mapHist, 1, &histSize, &histRange);

	// �p��CDF(�ֿn�������) - �����Ȫ��ֿn���v����
	photoHist /= photo.total();
	mapHist /= map.total();
	for (int i = 1; i < histSize; i++)
	{
		photoHist.at<float>(i) += photoHist.at<float>(i - 1); // ??
		mapHist.at<float>(i) += mapHist.at<float>(i - 1);
	}

	// �d���(LUT) - 1D�ƲաA�Ω��x�s�C�ӹ����Ȫ��M�g���Y
	cv::Mat lut(1, 256, CV_8U); // �@�w�nCV_8U?
	int refidx = 0;
	for (int srcidx = 0; srcidx < histSize; ++srcidx)
	{
		while (refidx < histSize && photoHist.at<float>(srcidx) > mapHist.at<float>(refidx))
		{
			refidx++;
		}
		lut.at<uchar>(srcidx) = refidx;
	}

	cv::Mat matched;
	cv::LUT(photo, lut, matched); // �ǹL�Ӫ��Ӥ��n��cv::Mat�����A������g�Lconver_to

	return matched;
}

// �Ӥ��e�B�z
void ImgLocate::PresetPhoto(cv::Mat& Map, cv::Mat& photo, double Flight_yaw, double s_factor)
{
	// �վ�Ӥ��ѪR��
	setResolution(photo, s_factor);

	// �L�H���ũ�Ӥ�/�a����Ƕ�
	grayImg(Map);
	grayImg(photo);
	//cv::equalizeHist(photo, photo); // �W�j�ũ�Ӥ��G�׹��A�Pmatlab�p��y���P

	// �վ�G�׹��
	adjustBright_Contrast(Map, photo, discard_count);

	// �վ�Ӥ����ײŦX�L�H��Yaw(�n���B�z�Ӥ��A��A���M�৹������|�v�T������/�зǮt)
	setAngle(photo, -Flight_yaw);
}

// �S�x�˴��P�y�z��
cv::Point2d ImgLocate::locatingAlgorithm(cv::Mat& Map, cv::Mat& photo, double yaw)
{
	// �Ыذ������Ѽ�
	bool extended_descriptor = true;			  // �X�i�y�z�l
	bool rotate = false;					   	  // �Ҽ{����
	float detect_threshold = 0.002f;			  // �˴����e
	int nOctaves = 3;							  // ���r��h��
	int nOctaveLayers = 4;						  // �C�h�����r��h�ơA�v�T�˴��ӽo�{��

	do
	{
		cv::Mat temp_Map = Map;
		cv::Mat temp_photo = photo;

		if (discard_count == 0)
		{
			PresetPhoto(temp_Map, temp_photo, yaw, 0.2);
			printf("�ũ�Ӥ��ѪR�׽վ㬰'0.2'��\n");
		}
		else if (discard_count == 1)
		{
			PresetPhoto(temp_Map, temp_photo, yaw, 0.2);
			detect_threshold = 0.001f;
			rotate = true;
			printf("�ũ�Ӥ��ѪR�׽վ㬰'0.2'�� �˴����e���C �Ҽ{����\n");
		}
		else if (discard_count == 2)
		{
			PresetPhoto(temp_Map, temp_photo, yaw, 0.3);
			detect_threshold = 0.001f;
			rotate = true;
			nOctaves = 4;
			printf("�ũ�Ӥ��ѪR�׽վ㬰'0.3'�� �˴����e���C �Ҽ{���� ���r��h�ƼW�[\n");
		}

		// KAZE�S�x�˴���
		cv::Ptr<cv::KAZE> kaze = cv::KAZE::create(
			extended_descriptor,		 // �O�_�ϥ��X�i�y�z�l
			rotate,						 // �O�_�Ҽ{����
			detect_threshold,			 // �˴����e�A�V�j�V�֯S�x�I
			nOctaves,					 // ���r��h�ơA�h�ثפ��R�A�ȶV�j�˴��d��V�s
			nOctaveLayers,				 // �C�h�����r��h�ơA�v�T�˴��ӽo�{��
			cv::KAZE::DIFF_CHARBONNIER	 // �X���ҫ�(DIFF_PM_G1: �����n�Ϲ��A���˴��L�z�S�x)(DIFF_PM_G2: �Ϲ���q���ι�ʯ��׭n�D������)
			// (DIFF_WEICKERT: ���c�Ƴ����A�p�ؿv��D�a��)(DIFF_CHARBONNIER: í�w�ʭn�D���������A�j�d��ǰt)
		);

		// �S�x�I/�y�z��
		vector<cv::KeyPoint> keyPoints1, keyPoints2;
		cv::Mat descriptor1, descriptor2;
		kaze->detect(temp_Map, keyPoints1);
		kaze->detect(temp_photo, keyPoints2);

		// ���Z���ƧǡA�U�p�ǰt�U�n�A�i��֭p��q�P��ֿ��~
		sort(keyPoints1.begin(), keyPoints1.end(), [](const cv::KeyPoint& a, const cv::KeyPoint& b) {return a.response > b.response; });
		sort(keyPoints2.begin(), keyPoints2.end(), [](const cv::KeyPoint& a, const cv::KeyPoint& b) {return a.response > b.response; });
		keyPoints1.resize(1800);
		if (keyPoints2.size() > 1800) { keyPoints2.resize(1800); } // �o�ˬO�_�N���|�X�{ Assertion Failed ���~

		//cv::KeyPointsFilter::retainBest(keyPoints1, 2000);
		//cv::KeyPointsFilter::retainBest(keyPoints2, 2000);

		// �y�z��
		kaze->compute(temp_Map, keyPoints1, descriptor1); // �ݭ���S�x�I�ƶq�A���M�|�]���e�Τj�q���s�ɭP�O����X��
		kaze->compute(temp_photo, keyPoints2, descriptor2);

		// �ǰt
		cv::BFMatcher matcher(cv::NORM_L2, true);
		vector<cv::DMatch> match_pts;
		matcher.match(descriptor1, descriptor2, match_pts);


		cv::Mat inliers, tform, H;
		vector<cv::Point2f> points1, points2;
		double averageDistance = 0.0;
		// �����ǰt�I
		for (const auto& matcher : match_pts)
		{
			points1.push_back(keyPoints1[matcher.queryIdx].pt);
			points2.push_back(keyPoints2[matcher.trainIdx].pt);
			averageDistance += matcher.distance;
		}

		// �ϥ�RANSAC�t��k
		//tform = cv::findHomography(points2, points1, cv::RANSAC, 5.0, inliers);  // 3.0���֭ȡA���j���֭ȥi�e�\��h�����I�Ctform���z���ܴ��x�}
																				 // (srcPoints, desPoints, method, ransacReprojThreshold, mask, maxIters, confidence)
		H = cv::estimateAffine2D(points2, points1, inliers, cv::RANSAC, 5.0); // 2x3���x�}�A�����ۦ��ܴ�
		tform = cv::Mat::eye(3, 3, CV_64F);
		H.copyTo(tform(cv::Rect(0, 0, 3, 2))); // �N�ۦ��ܴ��x�}�ɦ�3x3

		double tform_min, tform_max;
		cv::minMaxLoc(tform, &tform_min, &tform_max);
		cout << "�̤j��: " << tform_max << " ��C��: " << cv::determinant(tform) << endl;

		// �Nkeypoints�z�אּRANSAC�t��k�p��X���Ī��ǰt�I
		vector<cv::DMatch> inlier_matche_pts;
		for (size_t i = 0; i < match_pts.size(); i++)
		{
			if (inliers.at<uchar>(i))  // �B�n(0��outliers; 1��inliers)
			{
				inlier_matche_pts.push_back(match_pts[i]);
			}
		}

		// �P�_�ǰt���G(�o��ӤӦh�귽)
		averageDistance /= match_pts.size();
		cout << "�����Z��" << averageDistance << endl;

		int matchedThreshold = 8;
		if (inlier_matche_pts.size() < matchedThreshold && discard_count < 3) // �į���Ӹ��p
		{
			k = true;
			discard_count++;
			if (discard_count == 3)
			{
				cout << "���~�F�T���A�˱�" << endl;
				k = false;
				discard_count = 0;
				//tform = cv::Mat::zeros(3, 3, CV_8U);
				cv::Point CenterPoints(0, 0);
				// �M�z���s
				temp_Map.release();
				temp_photo.release();
				descriptor1.release();
				descriptor2.release();
				keyPoints1.clear();
				keyPoints2.clear();
				points1.clear();
				points2.clear();
				inliers.release();
				return CenterPoints;
			}
			cout << "�S�x�I�ƶq�L�֭��s�ǰt" << endl;
			continue;
		}
		else if (tform.empty() || tform_max > 1500 || cv::determinant(tform) < 1e-5) // �į���Ӹ��j
		{
			k = true;
			discard_count++;
			if (discard_count == 3)
			{
				cout << "���~�F�T���A�˱�" << endl;
				k = false;
				discard_count = 0;
				//tform = cv::Mat::zeros(3, 3, CV_8U);
				cv::Point CenterPoints(0, 0);
				// �M�z���s
				temp_Map.release();
				temp_photo.release();
				descriptor1.release();
				descriptor2.release();
				keyPoints1.clear();
				keyPoints2.clear();
				points1.clear();
				points2.clear();
				inliers.release();
				return CenterPoints;
			}
			cout << "�S���ܴ��x�}���ũα���ƹL�j�A�ǰt��í�w" << endl;
			continue;
		}
		else
		{
			// �ǰt���G
			cv::cvtColor(temp_Map, temp_Map, cv::COLOR_GRAY2BGR); // �ন3�q�D
			cv::cvtColor(temp_photo, temp_photo, cv::COLOR_GRAY2BGR);
			cv::Mat img_matches;
			cv::drawMatches(temp_Map, keyPoints1, temp_photo, keyPoints2, match_pts, img_matches);
			cv::namedWindow("Kaze Matches", cv::WINDOW_NORMAL);
			cv::imshow("Kaze Matches", img_matches);
			cv::waitKey(800);

			cv::drawMatches(temp_Map, keyPoints1, temp_photo, keyPoints2, inlier_matche_pts, img_matches);
			cv::namedWindow("Inlier Kaze Matches", cv::WINDOW_NORMAL);
			cv::imshow("Inlier Kaze Matches", img_matches);
			cv::waitKey(800);

			// �|��
			// 1.�p���ܴ��᪺�Ϲ����
			vector<cv::Point2d> newphotoCorners = { cv::Point2d(0,0), cv::Point2d(temp_photo.cols,0),
													cv::Point2d(temp_photo.cols, temp_photo.rows), cv::Point2d(0, temp_photo.rows) };
			vector<cv::Point2d> transformedCorner;
			cv::perspectiveTransform(newphotoCorners, transformedCorner, tform); // �z���ܴ���Ϲ��|�Ө��I�y�СA��ܦbmap�W���s�y��
			// 2.�p��s����ɮ�
			// �z���ܴ��᪺photo�i��W�Xmap��ɮت��j�p�A���ˬd�sphoto�|�Ө��I����ɡAmin�ˬd�O�_�W�X�t��V�Amax�ˬd�O�_�W�X����V
			double x_min = min({ transformedCorner[0].x, transformedCorner[1].x, transformedCorner[2].x, transformedCorner[3].x, 0.0 });
			double x_max = max({ transformedCorner[0].x, transformedCorner[1].x, transformedCorner[2].x, transformedCorner[3].x, (double)temp_Map.cols });
			double y_min = min({ transformedCorner[0].y, transformedCorner[1].y, transformedCorner[2].y, transformedCorner[3].y, 0.0 });
			double y_max = max({ transformedCorner[0].y, transformedCorner[1].y, transformedCorner[2].y, transformedCorner[3].y, (double)temp_Map.rows });
			// 3.�p��s�������ؤo
			int newWidth = static_cast<int>(x_max - x_min);
			int newHeight = static_cast<int>(y_max - y_min);
			// 4.�����ܴ�
			// ������Q�X�i�ɡA���I(0,0)�i��ݭn��۰���
			cv::Mat translation = (cv::Mat_<double>(3, 3) << 1, 0, -x_min,
				0, 1, -y_min,
				0, 0, 1);
			tform = translation * tform;
			// 5.�z���ܴ�
			cv::Mat warped_img;
			cv::warpPerspective(temp_photo, warped_img, tform, cv::Size(newWidth, newHeight));
			// 6.�N�ؼйϹ��K��s�Ϲ���
			cv::Mat stitchedImg = cv::Mat::zeros(newHeight, newWidth, temp_photo.type());
			cv::Mat roi = stitchedImg(cv::Rect(-x_min, -y_min, temp_Map.cols, temp_Map.rows));
			temp_Map.copyTo(roi);
			// 7.�|�[
			double alpha = 0.5; // �z����
			cv::addWeighted(stitchedImg, alpha, warped_img, 1 - alpha, 0.0, stitchedImg);

			// �p���|�Ϥ����I
			double Xmdn = ((transformedCorner[0].x + transformedCorner[1].x + transformedCorner[2].x + transformedCorner[3].x) / 4.0);
			double Ymdn = ((transformedCorner[0].y + transformedCorner[1].y + transformedCorner[2].y + transformedCorner[3].y) / 4.0);
			cv::Point2d CenterPoints(Xmdn, Ymdn);
			cv::circle(stitchedImg, CenterPoints, 5, cv::Scalar(0, 255, 255), -1);

			cv::namedWindow("blended", cv::WINDOW_NORMAL);
			cv::imshow("blended", stitchedImg);
			cv::waitKey(800);

			// �M�z���s
			cv::destroyWindow("Kaze Matches");
			cv::destroyWindow("Inlier Kaze Matches");
			cv::destroyWindow("blended");
			temp_Map.release();
			temp_photo.release();
			descriptor1.release();
			descriptor2.release();
			keyPoints1.clear();
			keyPoints2.clear();
			points1.clear();
			points2.clear();
			img_matches.release();
			inliers.release();
			warped_img.release();
			stitchedImg.release();

			// ���]�P�_�Ѽ�
			discard_count = 0;
			k = false;

			return CenterPoints; // �Y�����ǰt���G�h��^�|�Ϥ����I
		}

	} while (k == true);

}

// ����g�n�סA�P��ڶZ���t��
cv::Point2d ImgLocate::PositionCalculation(int blockX, int blockY, int j, double Flight_pitch, double Flight_yaw, cv::Point2d CenterPoints)
{
	cv::Point2d CropMap_TL_Pix(CropSize * (blockX - 2), CropSize * (blockY - 2)); // Puzzle�a�Ϧ^�_��j�a�ϭp��
	double Estimated_Lat, Estimated_Lon;
	Estimated_Lat = TL_coor[1] - (CropMap_TL_Pix.y + CenterPoints.y) * Lat_per_pix;
	Estimated_Lon = TL_coor[0] + (CropMap_TL_Pix.x + CenterPoints.x) * Lon_per_pix;

	// ���Y���׮ե�
	double uav_alt = 120; // �׭�
	double gimbal_pitch = 8.5;
	if (j >= 46) { gimbal_pitch = 0; }

	//double uav_alt = 60; // �{��
	//double gimbal_pitch = 0;

	double total_pitch = gimbal_pitch + Flight_pitch;
	double delta_d = uav_alt * tan(total_pitch * M_PI / 180);
	double delta[2] = { delta_d * sin(Flight_yaw * M_PI / 180), -delta_d * cos(Flight_yaw * M_PI / 180) }; // m

	double deltay2lat = 0.00898 * delta[1] / 1000; // 1000���ؼv�T 0.00898 �n��

	// 1000���ؼv�T (1000/(111320*cosd(��a�n��))) ���g�סA��Ī�[����n�׬�24.26
	double deltax2lon = (1000 / (111320 * cos(24.26 * M_PI / 180))) * delta[0] / 1000;
	//double deltax2lon = (1000 / (111320 * cos(24.17 * M_PI / 180))) * delta[0] / 1000; // �{��

	// �p��g�ե��᪺��m
	cv::Point2d Estimated_pos(Estimated_Lon - deltax2lon, Estimated_Lat + deltay2lat);

	return Estimated_pos;
}