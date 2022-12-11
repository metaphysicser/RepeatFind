#include "saca-k.hpp"
#include "binaryTree.hpp"
#include <stdint.h>
#include <time.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <algorithm>
#include <vector>
/*
本文件是对论文Efficient computation of all perfect repeats in genomic sequences of 
up to half a gigabyte, with a case study on the human genome的复现
Becher, Verónica et al. “Efficient computation of all perfect repeats in genomic sequences of up to half a gigabyte, 
with a case study on the human genome.” Bioinformatics (Oxford, England) vol. 25,14 (2009): 1746-53. doi:10.1093/bioinformatics/btp321
*/

// 比较函数，指定按照pair的first元素进行比较
bool compareByFirst(const std::pair<unsigned int, unsigned int>& a, const std::pair<unsigned int, unsigned int>& b) {   
        return a.first < b.first;
    
}

/*
@param: seq: 输入的拼接数据，末尾必须是0，论文中的w
        ml：min_length, 寻找的repeat的最小长度

本函数是基于后缀数组的repeat寻找算法
*/
void RepeatFind(std::string seq, unsigned int ml, std::string output_path)
{
    if(seq[seq.size()-1] != '0') // 判断是否0结尾
    {
        std::cout<< "The merged data should end with 0, please check your input"<<std::endl;
        exit(0);
    }

    if(ml < 0)
    {
        std::cout<< "The parameter ml should be bigger than 0, please check your input"<< std::endl;
        exit(0);
    }

    ofstream output;
    output.open(output_path, std::ios::out | std::ios::app); //以写入和在文件末尾添加的方式打开.txt文件，没有的话就创建该文件。
    if (!output.is_open())
    {
        std::cout<< "can not create ouput file, please check the output path"<<std::endl;
        exit(0);
    }

    unsigned int i, j, k;

    //--------------------构建后缀数组------------------------------------------------------------
    clock_t start, end;
    start = clock(); // 开始时间戳
    unsigned int n = seq.size(); // 输入序列的长度
    unsigned char * s_ch = (unsigned char *)seq.c_str();  // 输入的数据
    unsigned int * Sa = new unsigned int[n]; // 论文中的r，构建的后缀数组，Sa[i]代表第i个后缀在seq中的位置
    SACA_K(s_ch, Sa, n, 256, n, 0);
    end = clock();                                                                                       // 结束时间戳
    std::cout << "Sort SA time consume: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << std::endl; // 耗费的时间
    //------------------------------------------------------------------------------------------

    start = clock();
    unsigned int * SaRank = new unsigned int[n];  //论文中的p，SA的逆序数组，SaRank[i]代表seq第i个字符对应后缀数组中的位置
    for(i = 0; i < n; i++)
        SaRank[Sa[i]] = i;


    std::vector<std::pair<unsigned int, unsigned int>> LCP_tmp;
    LCP_tmp.resize(n);
    unsigned int * LCP = new unsigned int[n]; // LCP数组，LCP[i]代表第i个后缀和第i+1个后缀的公共前缀长度
   
    // n 是 seq 的长度
    // Sa 是排名为 i 的后缀的位置
    // SaRank 是后缀的排名
    // seq 是要求 LCP 的序列
    // LCP 是存储后缀 i 和后缀 i+1 的 LCP 的数组
    for (int i = 0, k = 0; i < n; ++i)
    {
        // 如果 i 是序列末尾，则跳过它
        if (SaRank[i] == n-1)
            continue;
    
        // 如果 k 不为 0，则减去 1
        if (k)
            --k;
    
        // 找到 LCP
        while (seq[i + k] == seq[Sa[SaRank[i] + 1] + k])
            ++k;
    
        // 存储 LCP
        LCP[SaRank[i]] = k;
        LCP_tmp[SaRank[i]].first = k;
        LCP_tmp[SaRank[i]].second = SaRank[i];
    }
    

     // 构建 n 个节点的二叉树
    TreeNode* S = buildTree(n, 0, LCP, ml); 
    // 设置内部节点的inS属性
    setInS(S);

    // find_index(S, n, n-1);

    std::sort(LCP_tmp.begin(), LCP_tmp.end(), compareByFirst); // 论文中的I数组
    unsigned int t, pi, ni;

    // 论文中的算法1
    for (t = 0; t < n; t++) 
    {
        const auto& item = LCP_tmp[t];
        if(item.first < ml)
            continue;
        i = item.second;
        pi = max_lessthan(S, n, i) + 1;
        ni = min_morethan(S, n, i);
        set_target_inS(S, n, i);
        if((pi==1 || LCP[pi-1]!=LCP[i])&&(ni==n || LCP[ni]!= LCP[i]))
        {
            if(Sa[pi]==0 || Sa[ni]==0 || seq[Sa[pi]-1] != seq[Sa[ni]-1] || labs(SaRank[Sa[ni]-1]-SaRank[Sa[pi]-1]) != ni - pi)
            {
                // 此时我们找到了ni-pi+1个长度为LCP[i]的公共子串，他们在seqence上的位置是Sa[m], pi<=m<=ni
                unsigned int wide = LCP[i];
                output << std::to_string(wide) << std::endl;
                for(unsigned int k = pi; k<=ni; k++)
                {
                    output<< std::to_string(Sa[k]) << " ";
                }
                output << std::endl;
            }
        }
    }

    output.close(); // 关闭文件
}