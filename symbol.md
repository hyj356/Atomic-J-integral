$w_m$: 第m个单元的应变能密度

$\sigma_{11}^m$:第m个单元在方向1上的应力分量

$u_1^m$:第m个单元在方向1上的变形量

$x_1$:直角坐标系下方向1的坐标

$x_2$:直角坐标系下方向2的坐标

$\sigma_{12}^m$:第m个单元在方向12上的应力分量, 也就是剪切应力分量

$u_2^m$:第m个单元在方向2上的变形量

$\sigma_{21}^m$:第m个单元在方向21上的应力分量, 其值在分子动力学中完全等于$\sigma_{12}^m$

$\sigma_{22}^m$:第m个单元在方向2上的应力分量

$g^m$:第m个单元的辅助函数$g(r) = \frac {(r-r_0)}{(r_1-r_0)}$ , 其中r为在极坐标系下的坐标, 等于点到原点的距离

$dA^m$:第m个单元的面积



​						对于函数$f(x)$, 若其在闭区间$[-1, 1]$上是连续的，那么有：
$$
\int_{-1}^{1} f(x) \approx f(\frac{1}{\sqrt[]{3} } )+f(- \frac{1}{\sqrt[]{3} } )
$$

$$
\iint_{D}f(x,y) = \int_{-1}^{1} \int_{-1}^{1} f(\xi ,\eta)d \xi d\eta = \int_{-1}^{1} (\int_{-1}^{1} f(\xi ,\eta)d \xi)d\eta\\
=f(\frac{1}{\sqrt[]{3} } , \frac{1}{\sqrt[]{3} } ) + f(-\frac{1}{\sqrt[]{3} } , \frac{1}{\sqrt[]{3} })+f(\frac{1}{\sqrt[]{3} } , -\frac{1}{\sqrt[]{3} }) + f(-\frac{1}{\sqrt[]{3} } , -\frac{1}{\sqrt[]{3} })
$$

$$
u(\xi, \eta)  = u_1N_1(\xi, \eta) + u_2N_2(\xi, \eta)+u_3N_3(\xi, \eta)+u_4N_4(\xi, \eta)\\
v(\xi, \eta)  = v_1N_1(\xi, \eta) + v_2N_2(\xi, \eta)+v_3N_3(\xi, \eta)+v_4N_4(\xi, \eta)
$$

$$
(\frac{1}{\sqrt[]{3} } , \frac{1}{\sqrt[]{3} } ) (-\frac{1}{\sqrt[]{3} } , \frac{1}{\sqrt[]{3} })(\frac{1}{\sqrt[]{3} } , -\frac{1}{\sqrt[]{3} })(-\frac{1}{\sqrt[]{3} } , -\frac{1}{\sqrt[]{3} })
$$

