#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <istream>
#include <iostream>
#include <string>
#include <algorithm>
#include <climits>
#include <vector>
// #include "../Auxiliary.h"
using namespace std;

vector<std::string> data;

/* Returns the file extension 返回文件扩展名
 */
const char *get_filename_ext(const char *filename)
{

  char *dot = const_cast<char *>(strrchr(filename, '.'));

  if (!dot || dot == filename)
    return "";

  return dot + 1;
}

// read sequences separeted by '>' line
// 读取fasta文件
vector<string> load_multiple_fasta(char *c_file, uint32_t *k)
{

  vector<string> c_buffer; // 分配最后的返回结果内存
  uint32_t i;

  ifstream f_in(c_file); // 打开文件

  if (!f_in)
  { // 读取文件失败
    fprintf(stderr, "file_open of '%s' failed: %s.\n", c_file, strerror(errno));
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < *k; i++)
  { // 每次读取一个序列
    string buf;

    getline(f_in, buf); // 读取每个序列第一行

    if (buf.length() == 0)
    { //如果没有数据，说明全部读取结束
      *k = i;
      break;
    }

    getline(f_in, buf, '>'); // 每次读取全部内容

    buf.erase(std::remove(buf.begin(), buf.end(), '\n'), buf.end());
    buf.erase(std::remove(buf.begin(), buf.end(), ':'), buf.end());

    c_buffer.push_back(buf);
  }

  f_in.close(); // 关闭文件

  return c_buffer;
}

/*******************************************************************/
// 目前只支持fasta文件
vector<string> file_load_multiple(char *c_file, uint32_t *k, int in_type, int ignore)
{

  /* .ext
   * .txt   - strings per line
   * .fasta - strings separated by '>' line
   * .fastq - strings separated by four lines
   * .gz    - use kseq parser to extract sequences
   * intype类型: 1是txt，2是fasta，3是fastq
   * ignore就是如果遇到其他后缀文件，该报错退出，还是返回空指针
   */

  vector<string> c_buffer; // 多个字符串序列
  if (*k == 0)
    *k = ULONG_MAX;

  const char *type = get_filename_ext(c_file); // 获取文件扩展名类型

  if (in_type)
  { // intype 是输入文件的类型

    // if(in_type==1) c_buffer = load_multiple_txt(c_file, k, n);
    if (in_type == 2)
      c_buffer = load_multiple_fasta(c_file, k);
    // else if(in_type==3)	c_buffer = load_multiple_fastq(c_file, k, n);
    else
    {
      exit(EXIT_FAILURE);
    }
  }
  else
  {

    // if(strcmp(type,"txt") == 0)
    //   c_buffer = load_multiple_txt(c_file, k, n);

    if (strcmp(type, "fasta") == 0 || strcmp(type, "fa") == 0 || strcmp(type, "fna") == 0)
      c_buffer = load_multiple_fasta(c_file, k);

    // else if(strcmp(type,"fastq") == 0 || strcmp(type,"fq") == 0)
    //   c_buffer = load_multiple_fastq(c_file, k, n);

    else
    {
      if (ignore)
      {
        *k = 0;
        //*n=0;
        return c_buffer;
      }
      printf("%s\n", c_file);
      printf("ERROR: file not recognized (.txt, .fasta, .fa, .fna, .fastq, .fq)\n");
      exit(EXIT_FAILURE);
    }
  }

  return c_buffer;
}

/*******************************************************************/

// /*******************************************************************/

// // read line by line
// char** load_multiple_txt(char *c_file, int *k, size_t *n) {

//   ifstream f_in;
//   f_in.open(c_file);

//   if(!f_in){
//     fprintf (stderr, "file_open of '%s' failed: %s.\n", c_file, strerror (errno));
//     exit (EXIT_FAILURE);
//   }

//   int n_alloc = N_ALLOC;
//   char **c_buffer = (char**) malloc(n_alloc*sizeof(char*));

//   int i;
//   for(i=0; i<*k; i++){

//     size_t len = 0; c_buffer[i] = NULL;
//     ssize_t size = getline(&c_buffer[i], &len, f_in);
//     if (size == -1){
//       *k = i;
//       break;
//     }
//     c_buffer[i][size-1] = 0;

//     (*n) += size;

//     if(i==n_alloc-1){
//       n_alloc+=N_ALLOC;
//       c_buffer = (char**) realloc(c_buffer, n_alloc*sizeof(char*));
//     }
//   }
//   f_in.close();

//   return c_buffer;
// }

// /*******************************************************************/

// // read sequences separeted by '@' line
// char** load_multiple_fastq(char *c_file, int *k, size_t *n){

//   ifstream f_in;
//   f_in.open(c_file);

//   if(!f_in){
//     fprintf (stderr, "file_open of '%s' failed: %s.\n", c_file, strerror (errno));
//     exit (EXIT_FAILURE);
//   }

//   int n_alloc = N_ALLOC;
//   char **c_buffer = (char**) malloc(n_alloc*sizeof(char*));

//   size_t len = 0;
//   char *buf = NULL;
//   int i;
//   for(i=0; i<*k; i++){

//     len = 0; buf = NULL;
//     ssize_t size = getline(&buf, f_in); // @'s line
//     free(buf);
//     if (size <= 1){
//       *k = i;
//       break;
//     }

//     len = 0; c_buffer[i] = NULL;
//     size = getline(&c_buffer[i], f_in, "\n"); // read line
//     c_buffer[i][size-1] = 0;

//     (*n) += size;

//     len = 0; buf = NULL;
//     size = getline(&buf, &len, f_in); // +'s line
//     free(buf);
//     len = 0; buf = NULL;
//     size = getline(&buf, &len, f_in); // @'s line
//     free(buf);

//     if(i==n_alloc-1){
//       n_alloc+=N_ALLOC;
//       c_buffer = (char**) realloc(c_buffer, n_alloc*sizeof(char*));
//     }
//   }
//   f_in.close();

//   return c_buffer;
// }

// /*******************************************************************/