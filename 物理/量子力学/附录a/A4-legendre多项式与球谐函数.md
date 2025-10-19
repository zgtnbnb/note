---
title: A4 Legendre多项式与球谐函数
tags:
  - 量子力学
  - 数学物理方法
date:
  "{ date }":
---


> **核心目标**：求解角动量算符 $\hat{L}^2$ 和 $\hat{L}_z$ 的共同本征函数，即球谐函数 $Y_{lm}(\theta, \varphi)$。这需要先解决球坐标下角度部分的微分方程，引出Legendre多项式和连带Legendre多项式。

## 角动量本征方程

在球坐标系中，轨道角动量算符 $(\hat{L}^2, \hat{L}_z)$ 的共同本征函数的角度部分 $\Theta(\theta)$ 满足以下微分方程（来自教材3.3.2节式(16)）：

$$
\frac{1}{\sin\theta} \frac{d}{d\theta}\left( \sin\theta \frac{d\Theta}{d\theta} \right) + \left( \lambda - \frac{m^2}{\sin^2\theta} \right)\Theta = 0\tag{1}
$$

其中：
- $\lambda$ 是 $\hat{L}^2$ 的本征值除以 $\hbar^2$，即 $\lambda = l(l+1)$。
- $m$ 是 $\hat{L}_z$ 的本征值除以 $\hbar$，即 $m = 0, \pm1, \pm2, \cdots$。

>这个方程是分离变量法的产物。它描述了波函数在球面上的“形状”，完全由角动量决定，与径向距离无关。理解这个方程是掌握原子结构和分子轨道的基础。

为了简化，我们引入变量代换：
令 $x = \cos\theta$，则 $|x| \leq 1$，且 $\Theta(\theta) = y(x)$。

代入后，方程(1)化为 **连带Legendre方程**：

$$
\frac{d}{dx}\left[ (1-x^2) \frac{dy}{dx} \right] + \left( \lambda - \frac{m^2}{1-x^2} \right)y = 0 \tag{2}
$$

当 $m=0$ 时，上式退化为标准的 **Legendre方程**：

$$
(1-x^2) \frac{d^2y}{dx^2} - 2x \frac{dy}{dx} + \lambda y = 0 \tag{3}
$$

>$x=\cos\theta$ 是一个非常聪明的代换。它把定义在 $[0, \pi]$ 上的 $\theta$ 映射到了 $[-1, 1]$ 区间，使得方程形式更简洁，也更容易应用幂级数解法。记住，$x=\pm1$ 对应于极点（北极和南极），是方程的正则奇点。

---

## A4.1 Legendre多项式

我们首先求解 $m=0$ 的情况，即Legendre方程(3)。

### 方法：幂级数解法

令解的形式为：

$$
y = \sum_{k=0}^{\infty} c_k x^k \quad (4)
$$

将其代入方程(3)，并比较同次幂项的系数，可以得到系数 $c_k$ 的递推关系：

$$
c_{k+2} = \frac{k(k+1) - \lambda}{(k+1)(k+2)} c_k \quad (5)
$$

>这是求解常微分方程的标准方法。递推关系决定了整个级数的结构。注意，这是一个二阶递推，意味着偶数项和奇数项是独立的，分别由 $c_0$ 和 $c_1$ 决定。

因此，通解可以写成两个线性独立解的叠加：

$$
\begin{cases}
y_1(x) = c_0 + c_2 x^2 + c_4 x^4 + \cdots \\
y_2(x) = c_1 x + c_3 x^3 + c_5 x^5 + \cdots
\end{cases} \quad (6)
$$

### 关键：物理边界条件

对于物理上可接受的解，我们需要波函数在全空间有界。这意味着在 $x = \pm1$ （即 $\theta = 0, \pi$）处，解不能发散。

从递推关系(5)可以看出，当 $k \to \infty$ 时，$c_{k+2}/c_k \sim 1 - 2/k$。这与 $\ln(1+x) + \ln(1-x) = \ln(1-x^2)$ 的泰勒展开的相邻项系数之比相同。而我们知道，$\ln(1-x^2)$ 在 $x=\pm1$ 处发散。

因此，为了保证解在 $|x|=1$ 处有界，我们必须让无穷级数在某一项之后截断，即变成一个**多项式**。

### 确定本征值 $\lambda$

为了让级数截断，递推关系(5)的分子必须在某个 $k=l$ 时为零：

$$
l(l+1) - \lambda = 0 \implies \lambda = l(l+1), \quad l = 0, 1, 2, \cdots \quad (7)
$$

> 这就是量子化的来源！$\lambda$ 不能取任意值，只能取离散的 $l(l+1)$。这里的 $l$ 就是我们熟悉的角量子数。这完美地体现了量子力学中“不连续”的特性。

当 $\lambda = l(l+1)$ 时，如果 $l$ 是偶数，则 $y_1(x)$ 是一个 $l$ 次多项式；如果 $l$ 是奇数，则 $y_2(x)$ 是一个 $l$ 次多项式。我们通常规定多项式的最高次项系数为：

$$
c_l = \frac{(2l)!}{2^l (l!)^2} \tag{8}
$$

这样得到的多项式称为 **Legendre多项式**，记作 $P_l(x)$。

其显式表达式为：

$$
P_l(x) = \sum_{k=0}^{\lfloor l/2 \rfloor} \frac{(-1)^k (2l-2k)!}{2^l k! (l-k)! (l-2k)!} x^{l-2k} \quad (9)
$$

> 这个公式看起来复杂，但它是从递推关系一步步推出来的。记住它的前几项即可：$P_0(x)=1$, $P_1(x)=x$, $P_2(x)=(3x^2-1)/2$。它们构成了 $[-1,1]$ 区间上的正交基。

### Legendre多项式的性质

#### 1. 奇偶性

$$
P_l(-x) = (-1)^l P_l(x) \quad (10)
$$

#### 2. Rodrigues公式

$$
P_l(x) = \frac{1}{2^l l!} \frac{d^l}{dx^l} (x^2 - 1)^l \quad (11)
$$

>Rodrigues公式是一个非常优美的表达式，它直接给出了多项式的微分形式。这在计算积分或证明正交性时非常有用。

#### 3. 生成函数

$$
(1 - 2xt + t^2)^{-1/2} = \sum_{l=0}^{\infty} P_l(x) t^l \quad (12)
$$

> 生成函数是特殊函数理论中的一个强大工具。它把一个无穷级数压缩成一个简单的闭合表达式，便于进行各种操作。

#### 4. 正交归一性

$$
\int_{-1}^{1} P_l(x) P_{l'}(x) dx = \frac{2}{2l+1} \delta_{ll'} \quad (13)
$$

#### 5. 递推关系

$$
\begin{aligned}
(l+1)P_{l+1} - (2l+1)xP_l + lP_{l-1} &= 0 \\
xP_l' - P_{l-1}' &= lP_l \\
P_{l+1}' - P_{l-1}' &= (2l+1)P_l \\
(x^2 - 1)P_l' &= lxP_l - lP_{l-1} \\
(2l+1)(x^2 - 1)P_l' &= l(l+1)(P_{l+1} - P_{l-1})
\end{aligned}\tag{14}
$$

---

## A4.2 连带Legendre多项式

现在我们回到一般情况，即 $m \neq 0$ 的连带Legendre方程(2)[[#角动量本征方程]]。

### 解法：利用Legendre多项式

我们猜测解的形式为：

$$
y(x) = (1-x^2)^{|m|/2} v(x) \quad (15)
$$
>好猜兄弟，好猜

代入方程(2)，经过一番微分运算，可以得到关于 $v(x)$ 的新方程：

$$
(1-x^2)v'' - 2(|m|+1)xv' + [\lambda - |m|(|m|+1)]v = 0 \quad (16)
$$

> 这个代换的物理意义在于，因子 $(1-x^2)^{|m|/2}$ 能够“吸收”掉原方程中因 $m \neq 0$ 而产生的奇异性，从而将问题转化为一个类似Legendre方程的形式。

对比方程(16)和Legendre方程(3)，我们发现它们的形式完全一样，只是参数不同。具体来说，方程(16)等价于一个“新的”Legendre方程，其本征值为 $\lambda - |m|(|m|+1)$。

因此，当 $\lambda = l(l+1)$ 时，方程(16)有物理上可接受的解的条件是：

$$
l(l+1) - |m|(|m|+1) = l'(l'+1)
$$

其中 $l'$ 是一个非负整数。显然，只有当 $l' = l - |m|$ 时才可能成立。这意味着 $l \geq |m|$。

于是，$v(x)$ 的解就是 $P_{l-|m|}(x)$，或者更准确地说，是 $P_l(x)$ 的 $|m|$ 阶导数。

最终，我们得到连带Legendre多项式的定义：

$$
P_l^{|m|}(x) = (1-x^2)^{|m|/2} \frac{d^{|m|}}{dx^{|m|}} P_l(x) \quad (17)
$$

> 这个定义非常直观。它表明，连带Legendre多项式是在普通Legendre多项式的基础上，乘上一个“权重”因子 $(1-x^2)^{|m|/2}$ 并求导 $|m|$ 次。这个“权重”因子正是为了满足角动量在 $z$ 方向的投影 $m$ 所带来的额外约束。

对于 $m \geq 0$，我们有：

$$
P_l^m(x) = (1-x^2)^{m/2} \frac{d^m}{dx^m} P_l(x) \quad (18)
$$

利用Rodrigues公式(11)[[#2. Rodrigues公式]]，我们可以写出：

$$
P_l^m(x) = \frac{1}{2^l l!} (1-x^2)^{m/2} \frac{d^{l+m}}{dx^{l+m}} (x^2 - 1)^l \quad (19)
$$

这个公式对 $m < 0$ 也同样适用，只需通过下面的关系式来定义：

$$
P_l^{-m}(x) = (-1)^m \frac{(l-m)!}{(l+m)!} P_l^m(x) \quad (20)
$$

> 负 $m$ 的定义是为了保证球谐函数的正交性和完备性。它本质上是利用了 $e^{im\varphi}$ 和 $e^{-im\varphi}$ 的关系。

#### 正交归一性

连带Legendre多项式也满足正交关系：

$$
\int_{-1}^{1} P_l^m(x) P_{l'}^m(x) dx = \frac{2}{2l+1} \frac{(l+m)!}{(l-m)!} \delta_{ll'} \quad (21)
$$

---

## A4.3 球谐函数

有了角度部分的解，我们现在可以构造完整的角动量本征函数——**球谐函数**。

### 定义

球谐函数 $Y_{lm}(\theta, \varphi)$ 定义为：

$$
Y_{lm}(\theta, \varphi) = (-1)^m \sqrt{\frac{2l+1}{4\pi} \frac{(l-m)!}{(l+m)!}} P_l^m(\cos\theta) e^{im\varphi} \quad (22)
$$

其中：
- $l = 0, 1, 2, \cdots$
- $m = -l, -l+1, \cdots, 0, \cdots, l-1, l$

> 球谐函数是量子力学中最重要的函数之一。它不仅是角动量算符的本征函数，也是氢原子波函数的角度部分。记住这个定义，因为它会贯穿整个量子力学课程。

### 性质

#### 1. 本征值
$$
\begin{aligned}
\hat{L}^2 Y_{lm} &= l(l+1)\hbar^2 Y_{lm} \\
\hat{L}_z Y_{lm} &= m\hbar Y_{lm}
\end{aligned} \quad (23)
$$

#### 2. 正交归一性

$$
\int_0^{2\pi} d\varphi \int_0^{\pi} \sin\theta d\theta Y_{lm}^*(\theta, \varphi) Y_{l'm'}(\theta, \varphi) = \delta_{ll'} \delta_{mm'} \quad (24)
$$

#### 3. 空间反射对称性

考虑空间反射 $\vec{r} \to -\vec{r}$，对应的球坐标变换为 $(r, \theta, \varphi) \to (r, \pi-\theta, \pi+\varphi)$。

可以证明：

$$
Y_{lm}(\pi-\theta, \pi+\varphi) = (-1)^l Y_{lm}(\theta, \varphi) \quad (25)
$$

>它表明，当 $l$ 为偶数时，球谐函数在空间反射下不变（偶宇称）；当 $l$ 为奇数时，球谐函数变号（奇宇称）。这在讨论原子光谱的跃迁选择定则时至关重要。

### 最简单的几个球谐函数

| $lm$  | $Y_{lm}(\theta, \varphi)$                                                 | $r^l Y_{lm}(\theta, \varphi)$                     |
| :---- | :------------------------------------------------------------------------ | :------------------------------------------------ |
| $00$  | $\frac{1}{\sqrt{4\pi}}$                                                   | $\frac{1}{\sqrt{4\pi}}$                           |
| $10$  | $\sqrt{\frac{3}{4\pi}} \cos\theta$                                        | $\sqrt{\frac{3}{4\pi}} z$                         |
| $1±1$ | $\mp \sqrt{\frac{3}{8\pi}} \sin\theta e^{\pm i\varphi}$                   | $\mp \sqrt{\frac{3}{8\pi}} (x \pm iy)$            |
| $20$  | $\sqrt{\frac{5}{16\pi}} (3\cos^2\theta - 1)$                              | $\sqrt{\frac{5}{16\pi}} (2z^2 - x^2 - y^2)$       |
| $2±1$ | $\mp \sqrt{\frac{15}{8\pi}} \cos\theta \cdot \sin\theta e^{\pm i\varphi}$ | $\mp \sqrt{\frac{15}{8\pi}} (x \pm iy)z$          |
| $2±2$ | $\frac{1}{2} \sqrt{\frac{15}{8\pi}} \sin^2\theta e^{\pm 2i\varphi}$       | $\frac{1}{2} \sqrt{\frac{15}{8\pi}} (x \pm iy)^2$ |

> 这些函数对应着s、p、d轨道。例如，$Y_{00}$ 是球对称的s轨道；$Y_{10}, Y_{1\pm1}$ 是p轨道；$Y_{20}, Y_{2\pm1}, Y_{2\pm2}$ 是d轨道。理解它们的几何形状对于化学和材料科学至关重要。

---

## A4.4 几个有用的展开式

### (a) 平面波的球面波展开

$$
e^{ikz} = e^{ikr\cos\theta} = \sum_{l=0}^{\infty} (2l+1)i^l j_l(kr) P_l(\cos\theta) = \sum_{l=0}^{\infty} \sqrt{4\pi(2l+1)} i^l j_l(kr) Y_{l0}(\theta) \quad (26)
$$

其中 $j_l(kr)$ 是球Bessel函数。

> 这个展开式在散射理论中极其重要。它表明，一个沿z轴传播的平面波，可以被看作是无数个具有不同角动量 $l$ 的球面波的叠加。

### (b) 库仑势的多极展开

两个点电荷之间的库仑势可以展开为：

$$
\frac{1}{|\vec{r}_1 - \vec{r}_2|} =
\begin{cases}
\frac{1}{r_2} \sum_{l=0}^{\infty} \left( \frac{r_1}{r_2} \right)^l P_l(\cos\theta), & r_1 < r_2 \\
\frac{1}{r_1} \sum_{l=0}^{\infty} \left( \frac{r_2}{r_1} \right)^l P_l(\cos\theta), & r_2 < r_1
\end{cases} \quad (27)
$$

其中 $\theta$ 是 $\vec{r}_1$ 与 $\vec{r}_2$ 之间的夹角。

> 这是静电学中的经典结果。它说明，远处的电势可以用一系列“多极矩”来描述，最低阶的是单极子（$l=0$），然后是偶极子（$l=1$）、四极子（$l=2$）等等。

### (c) 球谐函数的加法定理

$$
P_l(\cos\theta) = \frac{4\pi}{2l+1} \sum_{m=-l}^{l} Y_{lm}^*(\theta_1, \varphi_1) Y_{lm}(\theta_2, \varphi_2) \quad (28)
$$

其中 $\theta$ 是两个方向 $(\theta_1, \varphi_1)$ 和 $(\theta_2, \varphi_2)$ 之间的夹角。

当两个方向相同时，即 $\theta=0$，利用 $P_l(1)=1$，可得：

$$
\sum_{m=-l}^{l} |Y_{lm}(\theta, \varphi)|^2 = \frac{2l+1}{4\pi} \quad (29)
$$

> 加法定理是球谐函数理论中最深刻的结果之一。它揭示了球谐函数作为旋转群不可约表示的内在联系。这个公式在计算角动量耦合和张量算符矩阵元时非常有用。