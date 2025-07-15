# Calculate w with u, v input

Created by Zilong Qin (<zilong.qin@gmail.com>), 2022/3/15, @GBA-MWF, Shenzhen

## From the incompressible equation

$$
\nabla\mathbf{V}=0
$$
$$
\frac {\partial u} {\partial x}+\frac {\partial v} {\partial y}+\frac {\partial w} {\partial z} = 0 \label{ref1}\tag{1}
$$

## The incompressible equation under $\sigma$-coordination

With Eq. \eqref{ref1} and Eq. \eqref{ref5}, we have the incompressible equation under $\sigma$-coordination:

$$
\frac {\partial u} {\partial x'}+\frac {\partial v} {\partial y'}-\frac{z_{x'}}{z_{\sigma}}\frac {\partial u} {\partial \sigma}-\frac{z_{y'}}{z_{\sigma}}\frac {\partial v} {\partial \sigma}+\frac{1}{z_{\sigma}}\frac {\partial w} {\partial \sigma} = 0
$$

From \eqref{ref2}:

$$
\frac{z_{x'}}{z_{\sigma}}=Z_t\frac {Z_t-z} {({Z_t-Z_s)}^2}\frac {\partial Z_s} {\partial x'}=\frac {Z_t-\sigma} {Z_t-Z_s}\frac {\partial Z_s} {\partial x'}=c_1
$$

$$
\frac{z_{y'}}{z_{\sigma}}=Z_t\frac {Z_t-z} {({Z_t-Z_s)}^2}\frac {\partial Z_s} {\partial y'}=\frac {Z_t-\sigma} {Z_t-Z_s}\frac {\partial Z_s} {\partial y'}=c_2
$$

$$
\frac{1}{z_{\sigma}}=\frac {Z_T} {Z_T-Z_s}=c_3
$$

$$
\frac {\partial c_1}{\partial \sigma} = \frac {-1} {Z_t-Z_s}\frac {\partial Z_s} {\partial x'}
$$

$$
\frac {\partial c_2}{\partial \sigma} = \frac {-1} {Z_t-Z_s}\frac {\partial Z_s} {\partial y'}
$$


Finally:

$$
\frac {\partial w} {\partial \sigma} = \frac {1}{c_3}(c_1\frac {\partial u} {\partial \sigma}+c_2\frac {\partial v} {\partial \sigma}-\frac {\partial u} {\partial x'}-\frac {\partial v} {\partial y'}) \tag{6}
$$

## Solve the w with u and v

Make a differential over $\sigma$:

$$
\frac {\partial^2 w} {\partial \sigma^2} = \frac {1}{c_3} (\frac {\partial c_1}{\partial \sigma} \frac {\partial u}{\partial \sigma} +\frac {\partial c_2}{\partial \sigma} \frac {\partial u}{\partial \sigma} +c_1\frac {\partial^2 u} {\partial \sigma^2}+c_2\frac {\partial^2 v} {\partial \sigma^2} - \frac { \partial (\frac {\partial u} {\partial x'}+\frac {\partial v} {\partial y'})} {\partial \sigma})
$$

And through a poisson equation:
$$
w=\nabla_{\sigma}^{-2}(\frac{1}{c_3} \frac {\partial c_1}{\partial \sigma} \frac {\partial u}{\partial \sigma} +\frac{1}{c_3} \frac {\partial c_2}{\partial \sigma} \frac {\partial u}{\partial \sigma} +\frac {c_1}{c_3}\frac {\partial^2 u} {\partial \sigma^2}+\frac {c_2}{c_3}\frac {\partial^2 v} {\partial \sigma^2} - \frac {1}{c_3}\frac { \partial (\frac {\partial u} {\partial x'}+\frac {\partial v} {\partial y'})} {\partial \sigma})
$$

## UnitTest

Assume a function that:

$$
u = zsin(x) \\
v = zsin(y)
$$

So that:
$$
\frac {\partial w}{ \partial z} = -\frac {\partial u}{ \partial x} -\frac {\partial v}{ \partial y} =-zcos(x)-zcos(y)
$$

And :
$$
w = -\frac {1}{2}z^2(cos(x)+cos(y))+c
$$

In the TF coordinate:
$$
z=\sigma+(1-
\frac{\sigma} {Z_{t}})Z_s(x',y')=(1-\frac{Z_s}{Z_t}) \sigma+Z_s
$$

So that:
$$
w= -\frac {1}{2}((1-\frac{Z_s}{Z_t}) \sigma+Z_s))^2(cos(x)+cos(y))+c
$$

And:
$$
\frac {\partial w}{ \partial \sigma} = -(cos(x)+cos(y)) ((1-\frac{Z_s}{Z_t}) \sigma+Z_s)(1-\frac{Z_s}{Z_t})
$$

$$
\frac {\partial^2 w}{ \partial \sigma^2} = -(cos(x)+cos(y))(1-\frac{Z_s}{Z_t})^2
$$

## UnitTest2
$$
\nabla \cdot A=\frac{1}{R^2} \frac{\partial}{\partial R}\left(R^2 A_R\right)+\frac{1}{R \sin \theta} \frac{\partial}{\partial \theta}\left(A_\theta \sin \theta\right)+\frac{1}{R \sin \theta} \frac{\partial A_\phi}{\partial \phi}
$$

If $\nabla \cdot A=0$:
$$
\frac{1}{R} \frac{\partial}{\partial R}\left(R^2 A_R\right)+\frac{1}{ \sin \theta} \frac{\partial}{\partial \theta}\left(A_\theta \sin \theta\right)+\frac{1}{\sin \theta} \frac{\partial A_\phi}{\partial \phi}=0
$$

