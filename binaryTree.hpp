#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <istream>

//---------------------------------------
struct node // 用于求取数组I
{
    unsigned int index;
    unsigned int value;
};

//-----------------------------------------
// 定义二叉树的节点结构
struct TreeNode {
    bool isleaf;
    unsigned int index; // 节点的值
    bool inS;        
    TreeNode* left; // 左子树
    TreeNode* right;// 右子树
    TreeNode* parent; // 父亲
    TreeNode()
    {
        isleaf = false;
        index = -1;
        inS = false;
		left = NULL;
		right = NULL;
        parent = NULL;
    }
};


// 构建 1 到 n 的二叉树
TreeNode* buildTree(unsigned int n, unsigned int value, unsigned int * LCP, unsigned int ml) {
    // 如果 n 为 0 或者 1，则直接返回空节点
    if (n <= 0) return NULL;
    // 否则，创建一个根节点，并递归构建左右子树
    TreeNode* root = new TreeNode();
    if (n==1)
    {
        root->isleaf = true;
        root->index = value;
        if(LCP[value] < ml)
            root->inS = true;
        else
            root->inS = false;
        return root;
    } 

    root->left = buildTree(n / 2, value, LCP, ml); // 构建左子树
    if(root->left != NULL)
        root->left->parent = root;
    root->right = buildTree(n - n / 2, value + n / 2, LCP, ml); // 构建右子树
    if(root->right != NULL)
        root->right->parent = root;
    return root;
}

// 设置非叶子节点是否在S中
// 该函数用于遍历整棵树，并判断每个非叶子节点是否在集合 S 中
// 若该节点的左右子节点中至少有一个在 S 中，则该节点也在 S 中
// 该函数返回该节点是否在 S 中
bool setInS(TreeNode* root)
{
    // 如果该节点是叶子节点，则直接返回它在 S 中的状态
    if(root->isleaf)
        return root->inS;

    bool left_res = setInS(root->left);
    bool right_res = setInS(root->right);

    // 如果该节点的左子节点或右子节点在 S 中，则该节点也在 S 中
    if(left_res || right_res)
    {
        root->inS = true;
        return true;
    }
    else
        return false;
}


// find_index: 递归地搜索树，并返回第 n 个索引的树节点。
// root: 树的根节点。
// n: 树的总节点数。
// index: 要搜索的索引。
TreeNode* find_index(TreeNode* root, unsigned int n, unsigned int index)
{
    if(root->isleaf) // 如果当前节点是叶节点，则返回当前节点。
        return root;
    
    if(index < n / 2) // 如果要搜索的索引在左子树中，则递归地搜索左子树。
        return find_index(root->left, n / 2, index);
    else // 否则，递归地搜索右子树。
        return find_index(root->right, n - n / 2, index - n / 2);
}


/**
 * set_target_inS - 设置指定值的inS属性为true
 * @root: 树的根节点
 * @n: 树中节点的总数
 * @value: 要设置的值
 */
void set_target_inS(TreeNode* root, unsigned int n, unsigned int value)
{
    // 设置根节点的inS属性为true
    root->inS = true;
    if(root->isleaf) 
        return;
    
    // 如果目标值小于等于总节点数的一半，则递归设置左子树
    if(value < n / 2)
        set_target_inS(root->left, n / 2, value);
    else
        // 否则，递归设置右子树
        set_target_inS(root->right, n - n / 2, value - n / 2);
} 

// 寻找最大的小于某个值的值, 本算法是论文中的算法2
unsigned int max_lessthan(TreeNode* root, unsigned int n, unsigned int target)
{
    // 先找到target所在的节点
    TreeNode* t = find_index(root, n, target);
    while(!t->inS)
    {
        if(t == t->parent->left) // 如果t是左孩子，就在他的左边寻找最右的节点
        {
            unsigned int count = 0;
            // 分界点就是从左孩子变成右孩子
            while(t->parent != NULL && t == t->parent->left)
            {
                t = t->parent;
                count++;
            }
            

            // 从上一个循环跳出来有两种情况，第一种是t变成了右孩子，第二种是t变成了根节点
            if(t->parent != NULL)
            {   
                // 如果是右孩子，就跳到左边的左孩子上
                t = t->parent;
                t = t->left;
            }
            else
            {   // 如果是根节点，说明找不到更小的了，返回target
                return target;
            }

            // 然后一直往右走
            while(count-- && t->right != NULL)
                t = t->right;
        }
        else
            t = t->parent;
    }

    while(!t->isleaf)
    { 
        if(t->right->inS)
            t = t->right;
        else
            t = t->left;
    }

    return t->index;

}

// 寻找最小的大于某个值的值
unsigned int min_morethan(TreeNode* root, unsigned int n, unsigned int target)
{
    
    TreeNode* t = find_index(root, n, target);
    while(!t->inS)
    {
        if(t == t->parent->right)
        {
            unsigned int count = 0;
            while(t->parent != NULL && t == t->parent->right)
            {
                t = t->parent;
                count++;
            }
            if(t->parent != NULL)
            {
                t = t->parent;
                t = t->right;
            }
            else
            {
                return target;
            }
            

            while(count-- && t->left != NULL)
                t = t->left;
        }
        else
            t = t->parent;
    }

    while(!t->isleaf)
    {
        if(t->left->inS)
            t = t->left;
        else
            t = t->right;
    }

    return t->index;

}
