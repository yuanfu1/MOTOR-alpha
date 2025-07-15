# Constraint

## GeoBalance Constraint

GeoBalance is the balance inbetween the Coreolis force and the pressure gradient force:

$$
-fv=-\frac{1}{\rho}\frac{\partial p}{\partial x}
$$

$$
fu=-\frac{1}{\rho}\frac{\partial p}{\partial y}
$$

## Hydrostatic balance

While the $p_s$ is the control variable:

From $p=\rho R_dT_v$ , and $\frac {dp} {dz} = -\rho g$:

$$
p=p_s e^ {\int {-\frac{g}{R_dT_v} dz}}
$$

Where, $T_v = T(1+0.608q_v)$, therefore:

$$
p=p_s exp( {\int {-\frac{g}{R_dT(1+0.608q_v)} dz}})
$$

And its derivative:

$$
\frac {\partial p} {\partial p_s}= exp( {\int {-\frac{g}{R_dT(1+0.608q_v)} dz}})
$$

$$
\frac {\partial p} {\partial T}=p_s exp( {\int {-\frac{g}{R_dT(1+0.608q_v)} dz}})\frac{g}{R_dT^2(1+0.608q_v)} 
$$

$$
\frac {\partial p} {\partial q_v}=p_s exp( {\int {-\frac{g}{R_dT(1+0.608q_v)} dz}})\frac{0.608g}{R_dT(1+0.608q_v)^2}
$$

The discrete form of the equation:
$$
P_1=P_s
$$

$$
P_2=P_1*exp(\frac{-g*(h_2-h_1)}{R_d * \frac {T_1+T_2} {2}*(1+0.608*\frac{q_{v1}+q_{v2}}{2})})
$$

$$
P_n=P_{n-1}*exp(\frac{-g*(h_n-h_{n-1})}{R_d * \frac {T_{n-1}+T_n} {2}*(1+0.608*\frac{q_{v{n-1}}+q_{vn}}{2})})
$$

So that, its jacobian is:

$$
\frac {\partial P_n} {\partial P_{n-1}}=exp(\frac{-g*(h_n-h_{n-1})}{R_d * \frac {T_{n-1}+T_n} {2}*(1+0.608*\frac{q_{v{n-1}}+q_{vn}}{2})})
$$

$$
\frac {\partial P_n} {\partial T_{n-1}}=P_{n-1}*exp(\frac{-g*(h_n-h_{n-1})}{R_d * \frac {T_{n-1}+T_n} {2}*(1+0.608*\frac{q_{v{n-1}}+q_{vn}}{2})})*\frac{g*(h_n-h_{n-1})}{R_d * \frac {(T_{n-1}+T_n)^2} {2}*(1+0.608*\frac{q_{v{n-1}}+ q_{vn}}{2})}
$$

$$
\frac {\partial P_n} {\partial T_{n}}=P_{n-1}*exp(\frac{-g*(h_n-h_{n-1})}{R_d * \frac {T_{n-1}+T_n} {2}*(1+0.608*\frac{q_{v{n-1}}+q_{vn}}{2})})*\frac{g*(h_n-h_{n-1})}{R_d * \frac {(T_{n-1}+T_n)^2} {2}*(1+0.608*\frac{q_{v{n-1}}+q_{vn}}{2})}
$$

$$
\frac {\partial P_n} {\partial q_{v{n-1}}}=P_{n-1}*exp(\frac{-g*(h_n-h_{n-1})}{R_d * \frac {T_{n-1}+T_n} {2}*(1+0.608*\frac{q_{v{n-1}}+q_{vn}}{2})})*\frac{0.608g*(h_n-h_{n-1})}{2*R_d * \frac {T_{n-1}+T_n} {2}*(1+0.608*\frac{q_{v{n-1}}+q_{vn}}{2})^2}
$$

$$
\frac {\partial P_n} {\partial q_{vn}}=P_{n-1}*exp(\frac{-g*(h_n-h_{n-1})}{R_d * \frac {T_{n-1}+T_n} {2}*(1+0.608*\frac{q_{v{n-1}}+q_{vn}}{2})})*\frac{0.608g*(h_n-h_{n-1})}{2*R_d * \frac {T_{n-1}+T_n} {2}*(1+0.608*\frac{q_{v{n-1}}+q_{vn}}{2})^2}
$$

$$
\delta P_n = \frac {\partial P_n} {\partial P_{n-1}} \delta P_{n-1} + \frac {\partial P_n} {\partial T_{n-1}} \delta T_{n-1} + \frac {\partial P_n} {\partial T_{n}} \delta T_{n} + \frac {\partial P_n} {\partial q_{v{n-1}}} \delta q_{v{n-1}} + \frac {\partial P_n} {\partial q_{vn}} \delta q_{vn}
$$

## Hydrostatic balance using $lnP$ as control variable

From $p=\rho R_dT_v$ , and $\frac {dp} {dz} = -\rho g$:

$$
\frac {dp}{p} = -\frac {g}{R_dT_v} dz \Rightarrow d(lnp)= -\frac {g}{R_dT_v} dz
$$

So that in discrete form::

$$
lnp_0=lnp_s
$$

$$
lnp_1=lnp_0-\frac {g}{R_d\frac{T_2+T_1}{2}(1+0.608\frac {qv_2+qv_1}{2})} dz
$$

$$
lnp_n=lnp_{n-1}-\frac {g}{R_d\frac{T_n+T_{n-1}}{2}(1+0.608\frac {qv_n+qv_{n-1}}{2})} dz
$$


$$
\frac {\partial lnp_n} {\partial lnp_{n-1}}=1
$$

$$
\frac {\partial lnp_n} {\partial T_{n-1}}=-\frac {g}{R_d\frac{T_n+T_{n-1}}{2}(1+0.608\frac {qv_n+qv_{n-1}}{2})} dz
$$

$$
$$