# RepeatFind

本项目是对论文"Efficient computation of all perfect repeats in genomic sequences of up to half a gigabyte, with a case study on the human genome"的复现，主要功能是对输入序列中重复片段的寻找。

## 使用方法

核心函数为"RepeatFind.hpp"文件中的void RepeatFind(std::string seq, unsigned int ml, std::string output_path)函数

- seq: string类型，是输入的序列，其末尾必须以"0"结尾
- ml: unsigned int类型，是重复片段的最短长度，至少要大于0
- output_path: string类型，是存放输出结果的文件

运行RepeatFind函数后，可以在output_path的文件存放计算结果，格式如下：

```
长度

位置1 位置2 位置3 ... 位置n

...(不停重复)
```

第一行只有一个长度，长度指的是重复片段的长度，第二行有很多个位置，位置指的是每一个子串在序列seq中的开始位置

两行是一个重复子串，不停循环。

---

当然也可以直接调用test.cpp文件，方便快捷的处理fasta文件

首先编译test.cpp文件

```shell
g++ test.cpp -o test
```

然后调用test可执行程序

```shell
./test input.fasta output.txt
```

如果是window系统，命令中的./test 应当更换为./test.exe

test文件接受两个参数，输入文件的路径和输出文件的路径。目前输入文件只支持fasta格式，会自动将所有序列拼接起来，用"0"作为分隔符，然后输入进RepeatFind函数。ml参数也就是最小阈值需要在test文件中手动修改。

## 引用论文

> Verónica Becher, Alejandro Deymonnaz, Pablo Heiber, Efficient computation of all perfect repeats in genomic sequences of up to half a gigabyte, with a case study on the human genome, *Bioinformatics*, Volume 25, Issue 14, 15 July 2009, Pages 1746–1753, https://doi.org/10.1093/bioinformatics/btp321