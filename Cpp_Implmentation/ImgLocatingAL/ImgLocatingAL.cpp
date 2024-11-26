
#include "Header.h"

using namespace std;


int main()
{
    // 儲存資料到記憶體
    ImgLocate IL;
    vector<string> ImgName = IL.getCamImgFiles();
    vector<vector<string>> ImgLog = IL.getCamImgLog();

    // 演算法參數初始化
    double AerialP_Lon[104] = { 0 };
    double AerialP_Lat[104] = { 0 };
    double Aerial_pixX = 0;
    double Aerial_pixY = 0;
    double Flight_Yaw[104] = { 0 };
    double Flight_Pitch[104] = { 0 };
    int activate_flag = 1; // 執行演算法旗標，目前為判斷無人機到拍攝點 d < 30
    int j = 0; // 演算法迴圈
    int idx = 0; // 讀取空拍照片用
    int block_X = 0;
    int block_Y = 0;
    int picture_num = 0;
    vector<string> log_FirstColumn;
    string cropped_name;
    cv::Point2d Center_points, Estimated_Points;
    cv::Point2d Aerial_pix0[104];


    // 提取空拍照片log內容
    size_t LengthLog = ImgLog.size();
    for (int i = 1; i < LengthLog; i++)
    {
        AerialP_Lon[i - 1] = stod(ImgLog[i][17]);
        AerialP_Lat[i - 1] = stod(ImgLog[i][16]);
        Flight_Yaw[i - 1] = stod(ImgLog[i][8]);
        Flight_Pitch[i - 1] = stod(ImgLog[i][9]);
        log_FirstColumn.push_back(ImgLog[i][0]);
        Aerial_pix0[i - 1].x = abs(AerialP_Lon[i - 1] - IL.BL_coor[0]) / IL.Lon_per_pix;
        Aerial_pix0[i - 1].y = IL.Map_Height - (abs(AerialP_Lat[i - 1] - IL.BL_coor[1]) / IL.Lat_per_pix);
    }

    // 讀取地圖(不一定要)
    cv::Mat img1 = cv::imread(IL.folderPath_map);
    if (img1.empty())
    {
        cerr << "無法讀取地圖" << endl;
        return -1;
    }

    // 在地圖上標註實際拍攝點位置
    for (int i = 0; i < LengthLog - 1; i++)
    {
        cv::Point2d drawp(Aerial_pix0[i].x, Aerial_pix0[i].y);
        cv::circle(img1, drawp, 22, cv::Scalar(0, 0, 255), 10);
        if (i < (LengthLog - 2)) { cv::line(img1, Aerial_pix0[i], Aerial_pix0[i + 1], cv::Scalar(0, 0, 255), 10); }
    }

    cv::namedWindow("Map", cv::WINDOW_NORMAL); // 可自己調整大小的視窗
    cv::imshow("Map", img1);
    cv::waitKey(1000); // 若()為0則是無限等待使用者輸入，1000則為1秒以此類推

    // 演算法執行，要想要用甚麼觸發演算法執行，傳照片過來? flag要打開(=1)多久?
    while (j <= ImgName.size())
    {
        printf("目前照片: %d.png\n", j + 1);
        // 讀取img2(上機要換成空拍照片)
        cv::Mat img2 = cv::imread(ImgName[j]);
        fs::path full_path = ImgName[j];
        fs::path photoName = full_path.filename();

        // 根據空拍照片檔案名字找到log檔對應位子
        string s1 = photoName.string();
        auto it = find(log_FirstColumn.begin(), log_FirstColumn.end(), s1);
        if (it != log_FirstColumn.end())
        {
            idx = distance(log_FirstColumn.begin(), it);
        }
        else
        {
            cout << "Target: " << s1 << " not found." << endl;
        }

        // 計算經緯度對應的pixel(c++坐標系是甚麼?)
        Aerial_pixX = abs(AerialP_Lon[idx] - IL.BL_coor[0]) / IL.Lon_per_pix;
        Aerial_pixY = IL.Map_Height - (abs(AerialP_Lat[idx] - IL.BL_coor[1]) / IL.Lat_per_pix);
        block_X = ceil(Aerial_pixX / IL.CropSize);
        block_Y = ceil(Aerial_pixY / IL.CropSize);
        picture_num = (IL.numcol / IL.CropSize) * (block_X - 1) + block_Y;

        // 選取離線地圖對應小圖
        cropped_name = to_string(picture_num) + ".png";
        fs::path fullCropPath = IL.CroppedMapPath + "/" + cropped_name;

        // 拼接對應 3*3 小地圖並存入資料夾
        cv::Mat Puzzle33 = IL.puzzleMap(picture_num);
        if (Puzzle33.empty() || img2.empty())
        {
            cout << "照片為空" << endl;
            break;
        }

        // 執行影像定位演算法(包含照片處理)，回傳疊圖中心點
        Center_points = IL.locatingAlgorithm(Puzzle33, img2, Flight_Yaw[j]);

        // 估算演算法計算之經緯度
        Estimated_Points = IL.PositionCalculation(block_X, block_Y, j, Flight_Pitch[j], Flight_Yaw[j], Center_points);

        cout << "Estimated lat: " << Estimated_Points.y << " ,Estimated lon: " << Estimated_Points.x << endl;
        cout << "Actual lat: " << AerialP_Lat[j] << " ,Actual lon: " << AerialP_Lon[j] << endl;

        // 畫出估計經緯度與實際經緯度
        cv::Point2d Estimated_Points2pix((Estimated_Points.x - IL.BL_coor[0]) / IL.Lon_per_pix,
            IL.Map_Height - (Estimated_Points.y - IL.BL_coor[1]) / IL.Lat_per_pix);
        int size = 20;
        cv::Point2d TL(Estimated_Points2pix.x - size, Estimated_Points2pix.y - size);
        cv::Point2d TR(Estimated_Points2pix.x + size, Estimated_Points2pix.y - size);
        cv::Point2d BL(Estimated_Points2pix.x - size, Estimated_Points2pix.y + size);
        cv::Point2d BR(Estimated_Points2pix.x + size, Estimated_Points2pix.y + size);
        cv::line(img1, TL, BR, cv::Scalar(0, 255, 255), 10);
        cv::line(img1, TR, BL, cv::Scalar(0, 255, 255), 10);

        cv::imshow("Map", img1);
        cv::waitKey(800);

        // 清理內存
        Puzzle33.release();
        img2.release();


        // 迴圈控制
        j++;
    }





    return 0;

}