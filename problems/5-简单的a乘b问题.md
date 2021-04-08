# 简单的a * b问题

## 题目描述

给定两个整数 a 和 b ，请计算出它们的乘积。

## 输入

输入共两行

第一行一个整数 a

第二行一个整数 b

数据范围 1 < a, b < 1e1000000

### 输出

输出一个整数表示整数 a 和 b 的乘积

### 样例输入

```
2
3
```

### 样例输出

```
6
```

## 题解

> Java & Python: 我怀疑你在刁难我.

首先因为数据太大了, 创建 `BitInteger` 肯定超时. 所以就肯定不能使用了.

因此我们需要使用 **FFT _(快速傅立叶变换)_** 算法来解决问题.

抄模板即可.

**注意**: 开数组开太大和开太小都会报错, 考虑到数据大小因此这里开 `3e6 + 10`

### Java 的 BigInteger 分析

`Java8` 对 BigInteger 的乘法进行了更改, 主要加入了两种算法: `Karatsuba` 和 `ToomCook3`, 提高了速度.

但是从字符串创建 BigInteger 对象的速度特别慢 _(大约16秒)_, 导致了超时.

## 代码

```cpp
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
using namespace std;

const double PI = acos(-1.0);
struct complex
{
    double l, r;
    complex(double ll = 0.0, double rr = 0.0) {
        l = ll; r = rr;
    }
    complex operator +(const complex& B) {
        return complex(l + B.l, r + B.r);
    }
    complex operator - (const complex& B) {
        return complex(l - B.l, r - B.r);
    }
    complex operator *(const complex& B) {
        return complex(l * B.l - r * B.r, l * B.r + B.l * r);
    }
};

/*
 * 进行FFT和IFFT前的反转变换。
 * 位置i和j（i二进制反转后位置）互换
 * len必须是2的幂
 */
void change(complex y[], int len) {
    int i, j, k;
    for (int i = 1, j = len / 2; i < len - 1; i++) {
        if (i < j) swap(y[i], y[j]);
        k = len / 2;
        while (j >= k) {
            j -= k;
            k >>= 1;
        }
        if (j < k) j += k;
    }
}
/*
 * 做FFT
 * len必须为2^k形式，
 * on==1时是DFT，on==-1时是IDFT
 */
void fft(complex y[], int len, int on) {
    change(y, len);
    for (int h = 2; h <= len; h <<= 1) {
        complex wn(cos(-on * 2 * PI / h), sin(-on * 2 * PI / h));
        for (int j = 0; j < len; j += h) {
            complex w(1, 0);
            for (int k = j; k < j + h / 2; k++) {
                complex u = y[k];
                complex t = w * y[k + h / 2];
                y[k] = u + t;
                y[k + h / 2] = u - t;
                w = w * wn;
            }
        }
    }
    if (on == -1) {
        for (int i = 0; i < len; i++) {
            y[i].l /= len;
        }
    }
}
const int MAXN = 3e6 + 10;
complex x1[MAXN], x2[MAXN];
char s1[MAXN], s2[MAXN];
int sum[MAXN];
int main()
{
    scanf("%s", s1); scanf("%s", s2);
    int len1 = strlen(s1);
    int len2 = strlen(s2);
    int len = 1;
    while (len < len1 * 2 || len < len2 * 2) len <<= 1;
    for (int i = 0; i < len1; i++) {
        x1[i] = complex(s1[len1 - 1 - i] - '0', 0);
    }
    for (int i = len1; i < len; i++) {
        x1[i] = complex(0, 0);
    }
    for (int i = 0; i < len2; i++) {
        x2[i] = complex(s2[len2 - 1 - i] - '0', 0);
    }
    for (int i = len2; i < len; i++) x2[i] = complex(0, 0);

    fft(x1, len, 1);
    fft(x2, len, 1);
    for (int i = 0; i < len; i++) {
        x1[i] = x1[i] * x2[i];
    }
    fft(x1, len, -1);

    //化简和进位
    for (int i = 0; i < len; i++) {
        sum[i] = (int)(x1[i].l + 0.5);
    }
    for (int i = 0; i < len; i++) {
        sum[i + 1] += (sum[i] / 10);
        sum[i] %= 10;
    }
    len = len1 + len2 - 1;
    while (sum[len] <= 0 && len > 0) len--;
    for (int i = len; i >= 0; i--) {
        printf("%c", sum[i] + '0');
    }
    printf("\n");
}
```
