# Transformation of terrain-following coordinate

Created by Zilong Qin (<zilong.qin@gmail.com>), 2022/3/30, @GBA-MWF, Shenzhen

## The coordinate transformation

Assume that the coordinate system in $\sigma$-coordinate is $x^{\prime}, y^{\prime}, \sigma$, and original coordinate is $x, y, z$. Since the height of $\sigma$ is not constant. We have the relationship:

$$
\left\{\begin{array}{c}
x=x^{\prime} \\
y=y^{\prime} \\
z=z\left(x^{\prime}, y^{\prime}, \sigma\right)
\end{array}\right.
$$

Here
$$
\sigma=Z_{t}\frac{z-Z_{s}(x,y)}{Z_{t}-Z_{s}(x,y)}
$$

$$
z=\frac{\sigma}{Z_{t}}(Z_{t}-Z_{s}(x',y'))+Z_{s}(x',y') =\sigma+(1-
\frac{\sigma} {Z_{t}})Z_s(x',y')
$$

From the chain rule of derivatives:
$$
\left\{\begin{array}{l}
\frac{\partial}{\partial x'}=\frac{\partial}{\partial x} \frac{\partial x}{\partial x'}+\frac{\partial}{\partial y} \frac{\partial y}{\partial x'}+\frac{\partial}{\partial z} \frac{\partial z}{\partial x'} \\
\frac{\partial}{\partial y'}=\frac{\partial}{\partial x} \frac{\partial x}{\partial y'}+\frac{\partial}{\partial y} \frac{\partial y}{\partial y'}+\frac{\partial}{\partial z} \frac{\partial z}{\partial y'} \\
\frac{\partial}{\partial \sigma}=\frac{\partial}{\partial x} \frac{\partial x}{\partial \sigma}+\frac{\partial}{\partial y} \frac{\partial y}{\partial \sigma}+\frac{\partial}{\partial z} \frac{\partial z}{\partial \sigma}
\end{array}\right.
$$
In matrix:
$$
\left[\begin{array}{c}
\frac{\partial}{\partial x'} \\
\frac{\partial}{\partial y'} \\
\frac{\partial}{\partial \sigma}
\end{array}\right]=\left[\begin{array}{ccc}
\frac{\partial x}{\partial x'} & \frac{\partial y}{\partial x'} & \frac{\partial z}{\partial x'} \\
\frac{\partial x}{\partial y'} & \frac{\partial y}{\partial y'} & \frac{\partial z}{\partial y'} \\
\frac{\partial x}{\partial \sigma} & \frac{\partial y}{\partial \sigma} & \frac{\partial z}{\partial \sigma}
\end{array}\right]\left[\begin{array}{c}
\frac{\partial}{\partial x} \\
\frac{\partial}{\partial y} \\
\frac{\partial}{\partial z}
\end{array}\right]=\mathbf{T}\left[\begin{array}{c}
\frac{\partial}{\partial x} \\
\frac{\partial}{\partial y} \\
\frac{\partial}{\partial z}
\end{array}\right]
$$

$$
\left[\begin{array}{c}
\frac{\partial}{\partial x} \\
\frac{\partial}{\partial y} \\
\frac{\partial}{\partial z}
\end{array}\right]=\mathbf{T}^{-1}\left[\begin{array}{c}
\frac{\partial}{\partial x'} \\
\frac{\partial}{\partial y'} \\
\frac{\partial}{\partial \sigma}
\end{array}\right]
$$

So that we have Jacobi matrix:
$$
T=\left[\begin{array}{ccc}
\frac{\partial x}{\partial x'} & \frac{\partial y}{\partial x'} & \frac{\partial z}{\partial x'} \\
\frac{\partial x}{\partial y'} & \frac{\partial y}{\partial y'} & \frac{\partial z}{\partial y'} \\
\frac{\partial x}{\partial \sigma} & \frac{\partial y}{\partial \sigma} & \frac{\partial z}{\partial \sigma}
\end{array}\right]=\left[\begin{array}{lll}
1 & 0 & z_{x^{\prime}} \\
0 & 1 & z_{y^{\prime}} \\
0 & 0 & z_{\sigma}
\end{array}\right]
$$

and:

$$
T^{-1}=\left[\begin{array}{ccc}
1 & 0 & -\frac{z_{x^{\prime}}}{z_{\sigma}} \\
0 & 1 & -\frac{z_{y^{\prime}}}{z_{\sigma}} \\
0 & 0 & \frac {1} {z_{\sigma}}
\end{array}\right]
$$

Therefore:
$$
\frac{\partial}{\partial x}=\frac{\partial}{\partial x'}\cdot1+\frac{\partial}{\partial y'} \cdot0-\frac{\partial}{\partial \sigma} \frac{\partial \sigma}{\partial x'}= \frac{\partial}{\partial x'}-\frac{z_{x'}}{z_{\sigma}}\frac{\partial}{\partial \sigma} \\
\frac{\partial}{\partial y}=\frac{\partial}{\partial x'} \cdot0+\frac{\partial}{\partial y'} \cdot1-\frac{\partial}{\partial \sigma} \frac{\partial \sigma}{\partial y'}=\frac{\partial}{\partial y'} -\frac{z_{y'}}{z_{\sigma}}\frac{\partial}{\partial \sigma} \\
\frac{\partial}{\partial z}=\frac{\partial}{\partial x'} \cdot0+\frac{\partial}{\partial y'} \cdot0+\frac{\partial}{\partial \sigma} \frac{\partial \sigma}{\partial z}=\frac{1}{z_{\sigma}} \frac{\partial}{\partial \sigma}
$$

From :
$$
\frac{\partial z}{\partial x'}=(1-
\frac{\sigma} {Z_{t}})\frac{\partial Z_s}{\partial x'}
$$

$$
\frac {\partial z} {\partial \sigma}=1-\frac {Z_s} {Z_t}
$$

$$
\frac {z_{x'}} {z_{\sigma}} = (1-\frac {\sigma} {Z_t})\frac {\partial Z_s} {\partial x'} \frac {Z_t} {Z_t-Z_s} = \frac {Z_t-\sigma} {Z_t-Z_s}\frac {\partial Z_s} {\partial x'}
$$

$$
\frac {z_{y'}} {z_{\sigma}} = (1-\frac {\sigma} {Z_t})\frac {\partial Z_s} {\partial x'} \frac {Z_t} {Z_t-Z_s} = \frac {Z_t-\sigma} {Z_t-Z_s}\frac {\partial Z_s} {\partial y'}
$$

$$
\frac {1} {z_{\sigma}}=\frac {Z_t} {Z_t-Z_s}
$$

So that:

$$
\frac{\partial}{\partial x}= \frac{\partial}{\partial x'}-\frac {Z_t-\sigma} {Z_t-Z_s}\frac {\partial Z_s} {\partial x'}\frac{\partial}{\partial \sigma} \\
\frac{\partial}{\partial y}=\frac{\partial}{\partial y'} -\frac {Z_t-\sigma} {Z_t-Z_s}\frac {\partial Z_s} {\partial y'}\frac{\partial}{\partial \sigma} \\
\frac{\partial}{\partial z}=\frac {Z_t} {Z_t-Z_s} \frac{\partial}{\partial \sigma}
$$

## For tansformation of the divergence of a vector field:

$$
\nabla \mathbf{A} = \nabla ' \mathbf{A} - \frac {1} {z_{\sigma}} \nabla ' z \frac {\partial \mathbf{A}} {\partial \sigma}
$$

$$
\nabla \mathbf{A} = \nabla ' \mathbf{A} - (\nabla '(\frac {1} {z_{\sigma}}  z \frac {\partial \mathbf{A}} {\partial \sigma})-z \nabla ' (\frac {1} {z_{\sigma}} \frac {\partial \mathbf{A}} {\partial \sigma}))
$$

To implement it in a integral form:


