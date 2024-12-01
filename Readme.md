![image-20241119210100900](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119210100900.png)

​	J积分的积分区域如上图所示。J积分计算公式：
$$
J=\sum_{m=1}^{M}\left[\left(w^{m}-\sigma_{11}^{m} \frac{\partial u_{1}^{m}}{\partial x_{1}}-\sigma_{12}^{m} \frac{\partial u_{2}^{m}}{\partial x_{1}}\right) \frac{d g^{m}}{d x_{1}}-\left(\sigma_{21}^{m} \frac{\partial u_{1}^{m}}{\partial x_{1}}+\sigma_{22}^{m} \frac{\partial u_{2}^{m}}{\partial x_{1}}\right) \frac{d g^{m}}{d x_{2}}\right] d A^{m}
$$
​	在上式中，$w^m$表示第m个单元的应变能密度，$\sigma_{11},\sigma_{12},\sigma_{22},$分别表示第m个单元在方向1和方向2以及方向12上的单元应力, $x_1,x_2$表示单元在方向1和方向2上的直角坐标，$u_1,u_2$则表示单元在方向1和方向2上的位移。$dA_m$为第m个单元的面积, g为一个辅助函数，其具体的值与点的极坐标r有关,  计算公式为:
$$
g(r) = \frac{r-r_0}{r_1-r_0}
$$
​	首先我们来看一下单元面积$dA_m$的计算方法:
$$
\mathrm{d} A^{m}=\frac{1}{2}\left[(r+d r)^{2}-r^{2}\right] \mathrm{d} \theta
$$
![image-20241119211025817](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119211025817.png)

​	三个参数$r, dr, d\theta$如上图所示.

​	接下来是应变能密度:

![image-20241119211126377](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119211126377.png)

​	上市中, $p^{im}$表示第m个单元中的第i个原子在变形之后的势能，相对应的，$p^{im}_0$表示初始状态下的原子势能。$V_m$为第m个单元的体积，可以用上面计算出来的$dA_{m}$乘以样品的厚度计算获得。

​	接下来是计算单元应力：

![image-20241119211357937](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119211357937.png)

​	看公式很简单，就是把在单元内部的所有原子的应力分量加起来除以单元体积就可以获得。但是这里存在一个疑问，分子动力学获得的原子应力单位不是应力单位，而是应力乘以体积单位，所以这里的应力究竟指的是分子动力学里面的原子应力，还是指实际应力单位的应力呢？通过量纲分析，我倾向于认为这里的应力指的是原子应力，不过我换了一种平均化的方法，因为我提前在lammps的in文件中将原子应力除以了体积，所以获得的应力单位是GPa，就是真正的应力单位，所以我写程序的时候是直接把第m个单元里面的所有原子的应力分量加起来，再除以第m个单元内部的原子数量N，最终结果就代表了第m个单元内部的应力分量数值。

​	接下来我们需要计算每个单元四个角点的在水平和竖直方向上的位移分量:

![image-20241119211900028](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119211900028.png)

​	这里作者采用了空间平均的方法，找到距离每个单元的4个节点3埃米以内的所有原子，然后将每个坐落在对应角点距离内的原子质量$m^{\alpha}$乘以原子位移$u^{\alpha}$, 将它们加起来，再除以这些原子的总质量，就得到了每个单位的4个节点的水平和竖直方向的位移。

​	接下来就是有限元分析中的的形函数:

![image-20241119212417733](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119212417733.png)

​	在二维有限元分析中，我们常常需要计算大量的二重积分，为了便于计算二次积分，我们需要将一个任意形状的四面体单元映射到一个标准的2x2正方形区域中，如下图所示:

![image-20241119213548047](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119213548047.png)

​	在上图中，作图中的坐标称为全局坐标或者直角坐标，用x和y表示，右图中的坐标系称为自然坐标系，用$\xi$和$\eta$表示。它们之间存在着如下的映射函数关系:
$$
x = P(\xi, \eta) = x_1N_1(\xi, \eta) + x_2N_2(\xi, \eta)+x_3N_3(\xi, \eta)+x_4N_4(\xi, \eta)\\
y = Q(\xi, \eta) = y_1N_1(\xi, \eta) + y_2N_2(\xi, \eta)+y_3N_3(\xi, \eta)+y_4N_4(\xi, \eta)
$$
​	其中$N_1(\xi, \eta)$到$N_4(\xi, \eta)$为线性四边形单元的形函数。四个形函数, 及其对应的对自然坐标系下的偏导数为:
$$
N_1(\xi, \eta) = \frac{1}{4}(1-\xi)(1-\eta), \frac{\partial N1}{\partial \xi} = -\frac{1}{4}(1-\eta), \frac{\partial N1}{\partial \eta} = -\frac{1}{4}(1-\xi) \\
N_2(\xi, \eta) = \frac{1}{4}(1+\xi)(1-\eta), \frac{\partial N1}{\partial \xi} = \frac{1}{4}(1-\eta), \frac{\partial N1}{\partial \eta} = -\frac{1}{4}(1+\xi)	\\
N_3(\xi, \eta) = \frac{1}{4}(1+\xi)(1+\eta) , \frac{\partial N1}{\partial \xi} = \frac{1}{4}(1+\eta), \frac{\partial N1}{\partial \eta} = \frac{1}{4}(1+\xi)\\
N_4(\xi, \eta) = \frac{1}{4}(1-\xi)(1+\eta), \frac{\partial N1}{\partial \xi} = -\frac{1}{4}(1+\eta), \frac{\partial N1}{\partial \eta} = \frac{1}{4}(1-\xi)
$$
​	可以看到这四个形函数和文献中描述的公式（8）形式一致。接下来我们需要计算原子的应变和辅助函数g的值：

![image-20241119213118553](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119213118553.png)

​	首先从上图的第一行可以看到，原子在**j**方向上的位移$\varepsilon$等于单元的4个节点在j方向上的位移$u_j^{im}$分别乘以4个形函数$N_i$对**x1**的偏导数之和。其中单元的4个节点在j方向上的位移$u_j^{im}$在公式（7）中已经解决，辅助函数g只需要把直角坐标换成极坐标之后，也只是简单的四则运算而已，难点在于后续的形函数$N(\xi,\eta)$偏导。首先我们需要明确，形函数N是自然坐标$\xi$和$\eta$的函数，但是自然坐标又是全局坐标$x$和$y$的函数, 即:
$$
x = x(\xi, \eta)\\
y = y(\xi, \eta)\\
N(\xi , \eta)  =N(x(\xi, \eta), y(\xi, \eta)) = N(x, y)
$$
​	所以形函数也可以是全局坐标x和y的函数， 即$N(\xi, \eta)=N(x, y)$。那么根据链式求导法则，那么我们可以得到:
$$
\frac{\partial N}{\partial \xi}  = \frac{\partial N}{\partial x}\frac{\partial x}{\partial \xi} +\frac{\partial N}{\partial y}\frac{\partial y}{\partial \xi}\\
\frac{\partial N}{\partial \eta}  = \frac{\partial N}{\partial x}\frac{\partial x}{\partial \eta} +\frac{\partial N}{\partial y}\frac{\partial y}{\partial \eta}
$$
将上式写成矩阵的形式:
$$
\begin{bmatrix}
 \frac{\partial N}{\partial \xi}\\
\frac{\partial N}{\partial \eta}
\end{bmatrix}
=
\begin{bmatrix}
  \frac{\partial y}{\partial \xi}& \frac{\partial x}{\partial \xi} \\
  \frac{\partial x}{\partial \eta}&\frac{\partial y}{\partial \eta}
\end{bmatrix}
\begin{bmatrix}
 \frac{\partial N}{\partial x}\\
\frac{\partial N}{\partial y}
\end{bmatrix}
$$
​	其中矩阵:
$$
\begin{bmatrix}
\frac{\partial x}{\partial \xi}&  \frac{\partial y}{\partial \xi}\\
\frac{\partial x}{\partial \eta}&\frac{\partial y}{\partial \eta}
\end{bmatrix}
$$
​	被称为雅可比矩阵**J**，这个矩阵与文献中的公式(10)完全一致：

![image-20241119220206737](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119220206737.png)而根据上面的公式，我们发现只要知道了自然坐标$\xi$和$\eta$那么就可以很容易的求出来雅可比矩阵**J**。这里以上述矩阵的第一个元素为例, 展示一下具体的计算公式：
$$
\frac{\partial x}{\partial \xi} = x_1\frac{\partial N_1(\xi, \eta)}{\partial \xi} + x_2\frac{\partial N_2(\xi, \eta)}{\partial \xi} + x_3\frac{\partial N_3(\xi, \eta)}{\partial \xi} +x_4\frac{\partial N_4(\xi, \eta)}{\partial \xi}\\
=-\frac{x_1}{4}(1-\eta)+\frac{x_2}{4}(1-\eta)+\frac{x_3}{4}(1+\eta)-\frac{x_4}{4}(1+\eta)
$$


​	同时我们发现，如果已知自然坐标$\xi$和$\eta$，求出来雅可比矩阵的逆矩阵$J^{-1}$之后, 可以得到如下等式:
$$
\begin{bmatrix}
 \frac{\partial N}{\partial x}\\
\frac{\partial N}{\partial y}
\end{bmatrix}
=
J^{-1}
\begin{bmatrix}
 \frac{\partial N}{\partial \xi}\\
\frac{\partial N}{\partial \eta}
\end{bmatrix}
$$
​	而4个形函数$N(\xi, \eta)$关于自然坐标$\xi$和$\eta$的表达式以及对应偏导我们都已经在上文中展示过了，所以利用上述等式我们可以轻松求出来4个形函数$N(\xi, \eta)$对全局坐标x和y的偏导。

​	那么现在让我们回顾一下单元应变和单元辅助函数g：

![image-20241119220117867](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119220117867.png)

​	计算过程就非常清楚了。那么获得单元应变和单元辅助函数g之后，根据公式11，我们就需要求出它们分别对水平坐标x1和竖直坐标x2的偏导了：

![image-20241119220446617](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119220446617.png)

​	上图中的公式推导利用的是一种叫高斯积分的方法：对于一个在区域**[-1, 1]**上的**连续且可导**的函数f(x)来说，有如下等式成立:
$$
\int_{-1}^{1} f(x) = f(-\frac{1}{\sqrt{3} } ) + f(\frac{1}{\sqrt{3} } ) + o(x^4)
$$
​	其中$o(x^4)$是一个较小的误差项，可以忽略，此外，通过线性缩放，我们可以将任意区间**[a, b]**转换为区域**[-1, 1]**, 所以我们有：
$$
\int_{-1}^{1} f(x) = f(-\frac{1}{\sqrt{3} } ) + f(\frac{1}{\sqrt{3} } )
$$
​	而对于二重积分来说，我们可以将二次积分视作一重积分之后再进行一次一重积分：
$$
\int_{a}^{b} \int_{c}^{d} f(x, y)dxdy = \int_{a}^{b} (\int_{c}^{d} f(x, y)dx)dy
$$
​	同样通过线性缩放，我们可以把[a, b], [c, d]区间变化到[-1, 1], [-1, 1]区间，也就是自然坐标系$(\xi, \eta)$的区间，那么我们就会有:

![image-20241119221358467](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119221358467.png)

​	然后我们再回过头来看看文献的公式（11）：

![image-20241119221426685](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119221426685.png)

​	可以看到，这就是一个非常简单的应用了高斯-勒让德积分法的二重积分，对应的4个自然坐标分别为:
$$
(\xi_1, \eta_1) = (\frac{1}{\sqrt 3},\frac{1}{\sqrt 3} )\\
(\xi_1, \eta_2) = (-\frac{1}{\sqrt 3},\frac{1}{\sqrt 3} )\\
(\xi_2, \eta_1) = (\frac{1}{\sqrt 3},-\frac{1}{\sqrt 3} )\\
(\xi_2, \eta_2) = (-\frac{1}{\sqrt 3},-\frac{1}{\sqrt 3} )
$$
​	而根据我们之前的推导，已知自然坐标$(\xi, \eta)$的具体的值，那么我们就可以轻易的计算出雅可比矩阵$J_{jacobi}$，以及单元在1和2方向上的应变$\varepsilon_1,\varepsilon_2 $, 以及辅助函数g(r)的值。那么求出对应值之后，带入公式(3)即可计算出样品的J积分:

![image-20241119221923515](C:\Users\huangyijing\AppData\Roaming\Typora\typora-user-images\image-20241119221923515.png)





