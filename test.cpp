#include "utils/file.hpp"
#include "RepeatFind.hpp"
#include <stdint.h>
#include <time.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <iostream>
#include <algorithm>

// argc是参数个数
// arcv是具体参数，argv[0]是可执行exe文件的路径，argv[1]是输入文件的路径
int main(int argc, char *argv[])
{
    clock_t start, end;
    char *data_path; // 文件路径
    std::string output_path; // 输出路径
    std::string line;     // 用于读取文件的中间变量
    uint32_t k = 0;  // 读取的文件中序列的数量

    data_path = argv[1]; // 数据路径
    output_path = argv[2]; // 输出文件路径


    ifstream f(data_path);

    std::cout << "The input data path is " << data_path << std::endl; // 打印文件路径
    std::cout << "The ouput file path is " << output_path << std::endl;

    if (f.good())
    {                                                       // 判断输入的文件路径文件是否存在
        data = file_load_multiple(data_path, &k, 2, false); // data_path是输入文件的路径，k是返回的序列条数, n是返回的所有序列的碱基数量
        cout << "The input data has been loaded" << endl;
    }
    else
    {
        cout << "ERROR! input data path may not exist" << endl; // 打印文件路径不存在
        return 0;                                               // 退出程序
    }
    std::cout<< "The seqence number is "<< k << std::endl;
    
    std::string merged_data = "";
    // 拼接序列
    for(int i = 0; i<k;i++){
        string tmp = data[i];
        merged_data.append(tmp); 
        merged_data.append("0");
    }

    // 输入的序列确保以 0 作为结尾
    start = clock();

    unsigned int ml = 2000;
   
    RepeatFind(merged_data, ml, output_path);
    
    end = clock();                                                                                       // 结束时间戳
    std::cout << "Total time consume: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << std::endl; // 耗费的时间

       
}