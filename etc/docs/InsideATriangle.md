To determine if a given point $v$ lies in the interior of a given triangle, consider an individual vertex, denoted $v_0$, and let ![v_1] and $v_2$ be the vectors from $v_0$to the other two vertices. Expressing the vector from $v_0$ to $v$ in terms of $v_1$ and $v_2$ then gives

$v=v_0 + a v_1 + b v_2$ where $a$ and $b$ are constants. Solving for $a$ and $b$ gives

| ![a](https://mathworld.wolfram.com/images/equations/TriangleInterior/Inline14.svg) | ![=](https://mathworld.wolfram.com/images/equations/TriangleInterior/Inline15.svg) | ![(det(v v_2)-det(v_0 v_2))/(det(v_1 v_2))](https://mathworld.wolfram.com/images/equations/TriangleInterior/Inline16.svg) | (2)  |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ---- |
| ![b](https://mathworld.wolfram.com/images/equations/TriangleInterior/Inline17.svg) | ![=](https://mathworld.wolfram.com/images/equations/TriangleInterior/Inline18.svg) | ![-(det(v v_1)-det(v_0 v_1))/(det(v_1 v_2)),](https://mathworld.wolfram.com/images/equations/TriangleInterior/Inline19.svg) | (3)  |

where

| ![ det(u v)=uxv=u_xv_y-u_yv_x ](https://mathworld.wolfram.com/images/equations/TriangleInterior/NumberedEquation2.svg) | (4)  |
| ------------------------------------------------------------ | ---- |
|                                                              |      |

is the [determinant](https://mathworld.wolfram.com/Determinant.html) of the matrix formed from the [column vectors](https://mathworld.wolfram.com/ColumnVector.html) $u$ and $v$. Then the point $v$ lies in the interior of the triangle if $a, b > 0$ and $a+b < 1$.

If the [convex hull](https://mathworld.wolfram.com/ConvexHull.html) of the triangle vertices plus the point $v_0$ is bounded by four points, the point $v_0$ lies outside the triangle. However, if it contains three points, the point $v_0$ may lie either in the interior or in the exterior.