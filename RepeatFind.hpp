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
#include <queue>
#include <map>


// 比较函数，指定按照pair的first元素进行比较
bool compareByFirst(const std::pair<unsigned int, unsigned int>& a, const std::pair<unsigned int, unsigned int>& b) {   
        return a.first < b.first;
    
}

/*
@param: seq: 输入的拼接数据，末尾必须是0，论文中的w
        ml：min_length, 寻找的repeat的最小长度

本函数是对论文Efficient computation of all perfect repeats in genomic sequences of 
up to half a gigabyte, with a case study on the human genome的复现
Becher, Verónica et al. “Efficient computation of all perfect repeats in genomic sequences of up to half a gigabyte, 
with a case study on the human genome.” Bioinformatics (Oxford, England) vol. 25,14 (2009): 1746-53. doi:10.1093/bioinformatics/btp321
*/
void RepeatFind1(std::string seq, unsigned int ml, std::string output_path)
{
    size_t memory = 0; // 使用的内存大小
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
    output.open(output_path, std::ios::out | std::ios::trunc); //如果文件存在，将其清零，没有的话就创建该文件。
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

    memory += sizeof(char) * n;
    std::cout << "The seqence size is " << memory / pow(2, 20) << " MB" << std::endl;

    unsigned char * s_ch = (unsigned char *)seq.c_str();  // 输入的数据
    unsigned int * Sa = new unsigned int[n]; // 论文中的r，构建的后缀数组，Sa[i]代表第i个后缀在seq中的位置
    

    SACA_K(s_ch, Sa, n, 256, n, 0);
    end = clock();                                                                                       // 结束时间戳
    std::cout << "Sort SA time consume: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << std::endl; // 耗费的时间
    //------------------------------------------------------------------------------------------

    start = clock();
    unsigned int * SaRank = new unsigned int[n];  //论文中的p，SA的逆序数组，SaRank[i]代表seq第i个字符对应后缀数组中的位置
    memory += n * (2 * sizeof(unsigned int)); // SA, SaRank占用大小

    for(i = 0; i < n; i++)
        SaRank[Sa[i]] = i;


    std::vector<std::pair<unsigned int, unsigned int>> LCP_tmp;
    LCP_tmp.resize(n);
    unsigned int * LCP = new unsigned int[n]; // LCP数组，LCP[i]代表第i个后缀和第i+1个后缀的公共前缀长度
    memory += n * (sizeof(unsigned int) + sizeof(std::pair<unsigned int, unsigned int>)); // LCP和LCP_tmp
   
    // n 是 seq 的长度
    // Sa 是排名为 i 的后缀的位置
    // SaRank 是后缀的排名
    // seq 是要求 LCP 的序列
    // LCP 是存储后缀 i 和后缀 i+1 的 LCP 的数组
    for (int i = 0, k = 0; i < n; ++i)
    {
        // 如果 i 是序列末尾，则跳过它
        if (SaRank[i] == n-1)
        {
            LCP[SaRank[i]] = 0;
            LCP_tmp[SaRank[i]].first = 0;
            LCP_tmp[SaRank[i]].second = SaRank[i];
            continue;
        }
            
    
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
    memory += 2 * n * sizeof(TreeNode); 

    std::sort(LCP_tmp.begin(), LCP_tmp.end(), compareByFirst); // 论文中的I数组
    unsigned int t, pi, ni, res_count;
    res_count = 0; // 统计公共子串个数

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
                res_count++;
                output << std::to_string(wide) << std::endl;
                for(unsigned int k = pi; k<=ni; k++)
                {
                    output<< std::to_string(Sa[k]) << " ";
                }
                output << std::endl;
            }
        }
    }
    std::cout<< "The total result number is " << std::to_string(res_count) << std::endl;
    std::cout << "The total used memory is " << memory / pow(2, 20) << " MB" << std::endl;
    /*
    大部分内存都耗在了二叉树那里，大约占据了80%的内存，
    等以后有时间了优化一下，结构体的内存对齐耗费太大了
    */
    

    output.close(); // 关闭文件
    delete Sa;
    delete SaRank;
    delete LCP;
}


/*
@param: seq: 输入的拼接数据，末尾必须是0，论文中的w
        ml：min_length, 寻找的repeat的最小长度

本函数是对论文Space-Efficient Computation of Maximal and Supermaximal Repeats in Genome Sequences的复现
*/
void RepeatFind2(std::string seq, unsigned int ml, std::string output_path)
{
    size_t memory = 0; // 使用的内存大小
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
    output.open(output_path, std::ios::out | std::ios::trunc); //如果文件存在，将其清零，没有的话就创建该文件。
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

    memory += sizeof(char) * n;
    std::cout << "The seqence size is " << memory / pow(2, 20) << " MB" << std::endl;

    unsigned char * s_ch = (unsigned char *)seq.c_str();  // 输入的数据
    unsigned int * Sa = new unsigned int[n]; // 构建的后缀数组，Sa[i]代表第i个后缀在seq中的位置
    memory += sizeof(unsigned int) * n;
    

    SACA_K(s_ch, Sa, n, 256, n, 0);
    
    end = clock();                                                                                       // 结束时间戳
    std::cout << "Sort SA time consume: " << (double)(end - start) / CLOCKS_PER_SEC << "s" << std::endl; // 耗费的时间
    //------------------------------------------------------------------------------------------

    unsigned int * SaRank = new unsigned int[n];  //SA的逆序数组，SaRank[i]代表seq第i个字符对应后缀数组中的位置
    unsigned int * LCP = new unsigned int[n];
    memory += sizeof(unsigned int) * n * 2;
    for(i = 0; i < n; i++)
        SaRank[Sa[i]] = i;

    // 构建LCP数组
    LCP[0] = 0;
    for (i = 0, k = 0; i < n; ++i)
    {
        if (SaRank[i] == 0)
            continue;
        if (k)
            --k;
        while (seq[i + k] == seq[Sa[SaRank[i] - 1] + k])
        {
            ++k;
        }
        LCP[SaRank[i]] = k;
    }

    /*
        0 -- 0 A -- 1 C -- 2 G -- 3 T -- 4
    */


   std::map<char, unsigned int> char2int;
   char2int['0'] = 0;
   char2int['A'] = 1;
   char2int['C'] = 2;
   char2int['G'] = 3;
   char2int['T'] = 4;
   char int2char[] = {'0', 'A', 'C','G','T'};
   int C[5]; // 存储每个字符在SA中第一次出现的位置
   unsigned int left,right,mid;
   
   C[0] = 0; // 寻找AGCT第一次出现的位置
   for(i = 1; i <= 4; i++)
   {
        left = 0;
        right = n-1;
        while(true)
        {           
            if(right == left + 1)
            {
                if(seq[Sa[right]] == int2char[i] && seq[Sa[left] != int2char[i]])
                {
                    C[i] = right; 
                }
                else
                    C[i] = -1; // 该字符没有出现过
                break;
            }
                        
            mid = (left + right) / 2;
            if(seq[Sa[mid]] >= int2char[i])
                right = mid;
            else
                left = mid;
        }
        
          
   }
   

   bool* B = new bool[n+1]; // 论文中的B数组
   memory += sizeof(bool) * (n+1);
   B[0] = 1;
   B[n] = 1;

   bool* B_rank = new bool[n]; // 用于判断一个区间是否是maximal的数组,只要某个区间内存在1，就说明是maximal
   memory += sizeof(bool) * n;
   B_rank[0] = 0;
   for(i = 1; i < n; i++)
   {
        if(Sa[i-1]==0 || Sa[i] == 0 || seq[Sa[i-1]-1] != seq[Sa[i]-1])
            B_rank[i] = 1;
        else
            B_rank[i] = 0;

   }

   std::vector<std::queue<std::pair<unsigned int, unsigned int>>> Q; 
   Q.resize(5); // 存储队列
   std::pair<unsigned int, unsigned int> p;


   j = n-1;
   for(int m = 4; m>=0; m--)
   {
        if(C[m] < 0)
            continue;
        else
        {
            p = make_pair(C[m], j);
            j = C[m] - 1;
            Q[m].push(p);
        }
        
   }
   

   unsigned int l = 0;
   int last_lb = -1;
   int last_idx = -1;
   unsigned int size_c[5];
   unsigned int res_count = 0;

   while(true)
   {
        // 如果所有队列都空了，就退出循环
        bool empty_all = true;
        for(i = 0; i < 5; i++)
        {
            empty_all = empty_all && Q[i].empty();
            size_c[i] = Q[i].size(); // 统计每个队列的长度
        }
        if(empty_all)
            break;

        for(i = 0; i < 5; i++)
        {
            while(size_c[i] > 0)
            {
                p = Q[i].front();
                
                Q[i].pop();
                size_c[i] -= 1;
                left = p.first;
                right = p.second;

                if(B[right+1]==0)
                {
                    B[right+1] = 1;
                    if(last_lb == -1)
                    {
                        last_lb = left;
                    }
                    last_idx = right + 1;

                    // ---------------get_interval------------------------------------
                    std::map<char, std::pair<unsigned int, unsigned int>> new_interval;
                    char cur_char;
                    

                    for(j = left; j <= right; j++)
                    {
                        if(Sa[j] == 0)
                            continue;
                        cur_char = seq[Sa[j]-1];
                        auto it = new_interval.find(cur_char);
                        if(it == new_interval.end())
                        {
                            new_interval[cur_char] = make_pair(n-1,0);
                        }
                        
                        
                        if(SaRank[Sa[j]-1] <= new_interval[cur_char].first)
                            new_interval[cur_char].first = SaRank[Sa[j]-1];
                        if(SaRank[Sa[j]-1] >= new_interval[cur_char].second)
                            new_interval[cur_char].second = SaRank[Sa[j]-1];
                        

                    }

                    for(auto &t : new_interval)
                        Q[char2int[t.first]].push(t.second);
                    
                }
                else if(last_idx == left)
                {
                    // process
                    if(l>=ml)
                    {
                        bool maximal = false;
                        for(k = last_lb+1; k<=right; k++)
                        {
                            if(B_rank[k] == 1)
                            {
                                maximal = true;
                                break;
                            }
                        }
                        if(maximal)
                        {
                            res_count++;
                            output << std::to_string(l) << std::endl;
                            for(k = last_lb; k<=right; k++)
                            {
                                output<< std::to_string(Sa[k]) << " ";
                            }
                            output << std::endl;
                        }
                        
                        
                    }
                        
                    last_lb = -1;
                    last_idx = -1;
                    // ---------------get_interval------------------------------------
                    std::map<char, std::pair<unsigned int, unsigned int>> new_interval;
                    char cur_char;
                    

                    for(j = left; j <= right; j++)
                    {
                        if(Sa[j] == 0)
                            continue;
                        cur_char = seq[Sa[j]-1];
                        auto it = new_interval.find(cur_char);
                        if(it == new_interval.end())
                        {
                            new_interval[cur_char] = make_pair(n-1,0);
                        }
                        
                        if(SaRank[Sa[j]-1] <= new_interval[cur_char].first)
                            new_interval[cur_char].first = SaRank[Sa[j]-1];
                        if(SaRank[Sa[j]-1] >= new_interval[cur_char].second)
                            new_interval[cur_char].second = SaRank[Sa[j]-1];
                        

                    }

                    for(auto &t : new_interval)
                        Q[char2int[t.first]].push(t.second);
                }
            }
            
        }
        l++;
        
   }
   output.close(); // 关闭文件
   std::cout<< "the result number is " << res_count << std::endl;
   std::cout << "The total used memory is " << memory / pow(2, 20) << " MB" << std::endl;

}