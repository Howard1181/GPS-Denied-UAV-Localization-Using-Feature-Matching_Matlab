#pragma once
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <filesystem>
#include <opencv2/opencv.hpp>


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

	// �B��禡
	cv::Mat im2double(const cv::Mat& img);


public:
	vector<string> getCamImgFiles();
	vector<vector<string>> getCamImgLog();
	string CroppedMapPath = "E:/�L�H���v���w��/CY�v���w���@/Cropped_Map_FU"; // ������|�]�i�H
	string folderPath_map = "E:/VS/ImgLocation/Map/8KUHD_FongYuan2.png";
	void ImgLicatingAL_init();
	void puzzleMap(int picture_num);
	cv::Mat setResolution(cv::Mat, double);
	cv::Mat setAngle(cv::Mat, double);



	// �ŧi�a��/�ũ�򥻸�T(�p��ʸ�?)
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


};

// �N�Ӥ��ഫ��double���禡
cv::Mat ImgLocate::im2double(const cv::Mat& img)
{
	cv::Mat img_double;
	img.convertTo(img_double, CV_64F, 1.0 / 255.0);
	return img_double;
}

// ��l�ƺt��k�Ѽ�(�b�禡��l���N����F)
void ImgLocate::ImgLicatingAL_init()
{

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
void ImgLocate::puzzleMap(int picture_num)
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
	cv::waitKey(500);
	cv::destroyWindow("FinalImg");

	// �N�����n���Ӥ��g�J��Ƨ�
	string puzzlePath = "E:/�L�H���v���w��/CY�v���w���@/�p�a�ϫ���/Puzzle33_cpp.png";
	cv::imwrite(puzzlePath, Puzzle_map33);
}

// �վ�Ӥ��ѪR�ר��
cv::Mat ImgLocate::setResolution(cv::Mat img2, double sFactor)
{
	cv::Mat resized_img;
	cv::resize(img2, resized_img, cv::Size(), sFactor, sFactor);
	return resized_img;
}

// �վ�Ӥ�����
cv::Mat ImgLocate::setAngle(cv::Mat img2, double angle)
{
	cv::Mat rotated_img;
	cv::Point2f center(img2.cols / 2.0F, img2.rows / 2.0F); // �p��Ϲ�����
	cv::Mat rotateMatrix = cv::getRotationMatrix2D(center, angle, 1.0); // �Ыر���x�}
	cv::Rect boundingbox = cv::RotatedRect(center, img2.size(), angle).boundingRect(); // �p������Ϲ��j�p
	rotateMatrix.at<double>(0, 2) += boundingbox.width / 2.0 - center.x; // �����H�ŦX���
	rotateMatrix.at<double>(1, 2) += boundingbox.height / 2.0 - center.y;


	cv::warpAffine(img2, rotated_img, rotateMatrix, boundingbox.size());
	cv::namedWindow("rotated_img", cv::WINDOW_NORMAL);
	cv::imshow("rotated_img", rotated_img);
	cv::waitKey(0);
	return rotated_img;
}