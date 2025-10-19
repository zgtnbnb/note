## A3.1 Hermite方程与级数解法
$Hermite$方程的标准形式为：  
$$
u'' - 2z u' + (\lambda - 1) u = 0 \tag{1}
$$  
其中 $\lambda$ 为参数。当 $\lambda$ 为非负整数时，方程有**多项式解**（$Hermite$多项式）。


### 1. 级数解法
设 $u(z) = \sum_{k=0}^{\infty} c_k z^k$，代入方程$(1)$得递推关系：  $$u'(z)=\sum^\infty_{k=0}kC_kz^{k-1}$$$$u''(z)=\sum^\infty_{k=0}k(k-1)C_kz^{k-2}=\sum^\infty_{k=0}(k+2)(k+1)C_{k+2}z^k$$
$$
c_{k+2} = \frac{2k - (\lambda - 1)}{(k+1)(k+2)} c_k \tag{3}
$$  
由此可得：  
- 偶次项系数：$c_{2k} = \frac{2^{k} (\lambda - 1)(\lambda - 3) \cdots (\lambda - 2k + 1)}{(2k)!} c_0$  
- 奇次项系数：$c_{2k+1} = \frac{2^{k} \lambda (\lambda - 2) \cdots (\lambda - 2k)}{(2k+1)!} c_1$  


### 2. 多项式解条件
当 $\lambda = 2n$（$n = 0,1,2,\dots$ 为非负整数）时，递推关系终止，解为多项式。此时：  
$$
u_n(z) = c_0 H_n(z) + c_1 \tilde{H}_n(z) \tag{4}
$$  
其中 $H_n(z)$ 和 $\tilde{H}_n(z)$ 分别为**第一类**和**第二类Hermite多项式**。


## A3.2 Hermite多项式的正交性与生成函数
### 1. 正交性
Hermite多项式在区间 $(-\infty, +\infty)$ 上关于权重函数 $e^{-z^2}$ 正交：  
$$
\int_{-\infty}^{\infty} H_m(z) H_n(z) e^{-z^2} dz = \sqrt{\pi} 2^n n! \delta_{mn} \tag{12}
$$  
#### 归一化详细步骤


要推导Hermite多项式的正交关系  
$$
\int_{-\infty}^{+\infty} H_m(\xi) H_n(\xi) e^{-\xi^2} d\xi = \sqrt{\pi} \, 2^n n! \, \delta_{mn},
$$  
我们采用**生成函数法**，逐步展开如下：


##### **步骤1：定义Hermite多项式的生成函数**  
Hermite多项式的生成函数为：  
$$
G(\xi, t) = e^{2t\xi - t^2} = \sum_{n=0}^\infty \frac{H_n(\xi)}{n!} t^n. \tag{1}
$$  
该式将Hermite多项式表示为幂级数的形式，其中 $H_n(\xi)$ 是第 n 阶Hermite多项式。


##### **步骤2：构造两个生成函数的乘积并引入权重函数**  
考虑两个生成函数 $G(\xi, s)$ 和 $G(\xi, t)$ 的乘积，乘以权重函数 $e^{-\xi^2}$ 后积分：  
$$
\int_{-\infty}^{+\infty} G(\xi, s) G(\xi, t) e^{-\xi^2} d\xi = \int_{-\infty}^{+\infty} e^{2s\xi - s^2} e^{2t\xi - t^2} e^{-\xi^2} d\xi. \tag{2}
$$  
合并指数项：  
$$
e^{2s\xi - s^2 + 2t\xi - t^2 - \xi^2} = e^{2(s+t)\xi - (s^2 + t^2 + \xi^2)}. \tag{3}
$$


##### **步骤3：整理指数并完成平方**  
将指数中的二次项配方：  
$$
2(s+t)\xi - \xi^2 = -\left( \xi^2 - 2(s+t)\xi \right) = -\left[ (\xi - (s+t))^2 - (s+t)^2 \right]. \tag{4}
$$  
代入式(3)得：  
$$
e^{2(s+t)\xi - s^2 - t^2 - \xi^2} = e^{-(\xi - (s+t))^2 + (s+t)^2 - s^2 - t^2}. \tag{5}
$$  
计算常数项：  
$$
(s+t)^2 - s^2 - t^2 = s^2 + 2st + t^2 - s^2 - t^2 = 2st. \tag{6}
$$  
因此，式(5)简化为：  
$$
e^{2st} e^{-(\xi - (s+t))^2}. \tag{7}
$$


##### **步骤4：计算高斯积分**  
对式(7)积分，利用高斯积分公式  $\int_{-\infty}^{+\infty} e^{-x^2} dx = \sqrt{\pi}$ ：  
$$
\int_{-\infty}^{+\infty} e^{2st} e^{-(\xi - (s+t))^2} d\xi = e^{2st} \int_{-\infty}^{+\infty} e^{-(\xi - (s+t))^2} d\xi = e^{2st} \cdot \sqrt{\pi}. \tag{8}
$$


##### **步骤5：展开生成函数的级数形式**  
左边生成函数乘积的级数展开为：  
$$
G(\xi, s) G(\xi, t) = \left( \sum_{m=0}^\infty \frac{H_m(\xi)}{m!} s^m \right) \left( \sum_{n=0}^\infty \frac{H_n(\xi)}{n!} t^n \right) = \sum_{m=0}^\infty \sum_{n=0}^\infty \frac{H_m(\xi) H_n(\xi)}{m! n!} s^m t^n. \tag{9}
$$  
乘以权重函数并积分后：  
$$
\int_{-\infty}^{+\infty} G(\xi, s) G(\xi, t) e^{-\xi^2} d\xi = \sum_{m=0}^\infty \sum_{n=0}^\infty \left( \int_{-\infty}^{+\infty} H_m(\xi) H_n(\xi) e^{-\xi^2} d\xi \right) \frac{s^m t^n}{m! n!}. \tag{10}
$$  
右边式(8)的级数展开为：  
$$
e^{2st} \sqrt{\pi} = \sqrt{\pi} \sum_{k=0}^\infty \frac{(2st)^k}{k!} = \sqrt{\pi} \sum_{k=0}^\infty \frac{2^k s^k t^k}{k!}. \tag{11}
$$


##### **步骤6：比较系数确定正交性**  
式(10)和式(11)均为双变量幂级数，根据幂级数唯一性定理，对应项系数必须相等。  
- 当 $m \neq n$ 时，左边 $s^m t^n$ 项的系数为  $\frac{1}{m! n!} \int_{-\infty}^{+\infty} H_m(\xi) H_n(\xi) e^{-\xi^2} d\xi$ ，而右边无此交叉项（仅含 $s^k t^k$ 形式），故系数为零，即：  
  $$
  \int_{-\infty}^{+\infty} H_m(\xi) H_n(\xi) e^{-\xi^2} d\xi = 0 \quad (m \neq n). \tag{12}
  $$  
- 当 $m = n = k$ 时，左边 $s^k t^k$ 项的系数为  $\frac{1}{(k!)^2} \int_{-\infty}^{+\infty} H_k(\xi)^2 e^{-\xi^2} d\xi$ ，右边对应项系数为  $\sqrt{\pi} \frac{2^k}{k!}$ ，故：  
  $$
  \frac{1}{(k!)^2} \int_{-\infty}^{+\infty} H_k(\xi)^2 e^{-\xi^2} d\xi = \sqrt{\pi} \frac{2^k}{k!}, \tag{13}
  $$  
  解得：  
  $$
  \int_{-\infty}^{+\infty} H_k(\xi)^2 e^{-\xi^2} d\xi = \sqrt{\pi} \, 2^k k!. \tag{14}
  $$


##### **步骤7：综合结果**  
结合式(12)和式(14)，正交关系可统一表示为：  
$$
\int_{-\infty}^{+\infty} H_m(\xi) H_n(\xi) e^{-\xi^2} d\xi = \sqrt{\pi} \, 2^n n! \, \delta_{mn},
$$  
其中  $\delta_{mn}$  是Kronecker delta 函数（m = n 时为1，否则为0）。
### 2. 生成函数
Hermite多项式的生成函数为：  
$$
e^{2zt - t^2} = \sum_{n=0}^{\infty} \frac{H_n(z)}{n!} t^n \tag{13}
$$  


## A3.3 Hermite多项式的递推关系与特殊值
### 1. 递推关系
- 微分递推：  
  $$
  H_{n+1}(z) = 2z H_n(z) - 2n H_{n-1}(z) \tag{8}
  $$  
- 积分递推：  
  $$
  \int_{-\infty}^{\infty} H_m(z) H_n(z) e^{-z^2} dz = \sqrt{\pi} 2^n n! \delta_{mn} \tag{12}
  $$  


### 2. 特殊值
- $H_0(z) = 1$  
- $H_1(z) = 2z$  
- $H_2(z) = 4z^2 - 2$  
- $H_3(z) = 8z^3 - 12z$  


## A3.4 Hermite多项式与量子力学
Hermite多项式是**量子谐振子**的本征函数。对于质量为 $\mu$、角频率为$\omega$ 的谐振子，能量本征值为：  
$$
E_n = \left( n + \frac{1}{2} \right) \hbar \omega, \quad n = 0,1,2,\dots
$$  
对应的波函数为：  
$$
\psi_n(x) = \left( \frac{\mu \omega}{\pi \hbar} \right)^{1/4} \frac{1}{\sqrt{2^n n!}} H_n\left( \sqrt{\frac{\mu \omega}{\hbar}} x \right) e^{-\frac{\mu \omega x^2}{2\hbar}}
$$
## hermite图像
![[hermite图像.png]]

