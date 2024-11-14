
#include "Header.h"

using namespace std;


int main()
{

    // 儲存資料到記憶體
    ImgLocate IL;
    IL.ImgLicatingAL_init();
    vector<string> ImgName = IL.getCamImgFiles();
    vector<vector<string>> ImgLog = IL.getCamImgLog();


    // 演算法參數初始化
    double AerialP_Lon[104] = { 0 };
    double AerialP_Lat[104] = { 0 };
    double Aerial_pixX = 0;
    double Aerial_pixY = 0;
    double Flight_Yaw[104] = { 0 };
    double Flight_Pitch[104] = { 0 };
    int k = 0; // 演算法重跑旗標
    int discard_count = 0; // 捨棄匹配次數
    int activate_flag = 1; // 執行演算法旗標，目前為判斷無人機到拍攝點 d < 30
    int j = 0; // 演算法迴圈
    int idx = 0; // 讀取空拍照片用
    int block_X = 0;
    int block_Y = 0;
    int picture_num = 0;
    vector<string> log_FirstColumn;
    string cropped_name;




    // 提取空拍照片log內容
    size_t LengthLog = ImgLog.size();
    for (int i = 1; i < LengthLog; i++)
    {
        AerialP_Lon[i - 1] = stod(ImgLog[i][17]);
        AerialP_Lat[i - 1] = stod(ImgLog[i][16]);
        Flight_Yaw[i - 1] = stod(ImgLog[i][8]);
        Flight_Pitch[i - 1] = stod(ImgLog[i][9]);
        log_FirstColumn.push_back(ImgLog[i][0]);
    }

    // 讀取地圖
    cv::Mat img1 = cv::imread(IL.folderPath_map);
    if (img1.empty())
    {
        cerr << "無法讀取地圖" << endl;
        return -1;
    }


    cv::namedWindow("Map", cv::WINDOW_NORMAL); // 可自己調整大小的視窗
    cv::imshow("Map", img1);
    cv::waitKey(1000); // 若()為0則是無限等待使用者輸入，1000則為1秒以此類推


    // 演算法執行，要想要用甚麼觸發演算法執行，傳照片過來? flag要打開(=1)多久?
    for (int j = 0; j <= ImgName.size(); j++)
    {
        printf("目前照片: %d.png\n", j + 1);
        // 讀取img2(上機要換成空拍照片)
        cv::Mat img2 = cv::imread(ImgName[j]);
        fs::path full_path = ImgName[j];
        fs::path photoName = full_path.filename();

        // 根據空拍照片檔案名字找到log檔對應位子
        string s1 = photoName.string();
        auto it = find(log_FirstColumn.begin(), log_FirstColumn.end(), s1); // 這行有問題
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
        IL.puzzleMap(picture_num);

        // 調整空拍照片解析度
        if (k == 1 && discard_count >= 3)
        {
            img2 = IL.setResolution(img2, 0.2);
            printf("空拍照片解析度調整為'0.2'倍\n");
        }
        else
        {
            img2 = IL.setResolution(img2, 0.3);
            printf("空拍照片解析度調整為'0.3'倍\n");
        }

        // 調整照片角度符合無人機Yaw
        img2 = IL.setAngle(img2, -Flight_Yaw[j]);

        // 無人機空拍照片轉灰階


    }





    return 0;

}


