# 当前主版本 (v0.1.4)

## Update 2025-06-04
By Zilong QIN
1. GRAPES 模型数据读取 (IOGrapesModelvar_500m.F90)，修正了单位转换的 bug，处理顺序优化：调整了 Calc_sg_vars 和 update_restrictionOfStatics 的调用位置；
2. 网格高度处理 (SingleGrid.F90)，修正了单垂直层情况下的高度设置问题，地形数据处理：改进了地形高度的初始化；
3. 为泊松方程求解器中的 calResidual 添加伴随；




## Update 2025-05-13
By Jiongming Pang
1. Fixed bugs in compile. When compiling on Ubuntu,the function names in the 'CMakeLists.txt' file are case-sensitive.

## Update 2025-05-12
By Jiongming Pang
1. Fixed bugs in compile of develop brand.
2. Hybrid-EnLoc parallel computing optimization.
3. Optimization of parallel reading modelvar ensemble files by BEC-EnLoc.

## Update 2025-05-07
1. 优化了yaml_get_var的结构，现在可以读取最多 4 层结构的yaml文件。


## Update 2025-04-09
1. 修正了 IOERA5 中的 配置文件读取 bug；
2. 修改了 RField 中雷达和探空风场的权重曲线；
3. 修改了 Application 中地面温度和地面气压的同化层：气压为 N-1 层开始同化，温度 N-2 层开始同化；
4. 增加了 AnalysisSummary.txt输出，以方面查看各层的分析尺度；

## Update 2025-03-21
By Zilong Qin
1. 调整了 RField 中雷达的权重，调整了 ObsRadarVelEach 中雷达的权重；
2. 优化了 UV2Divg 的代码，构造测试了 UV2Divg 的伴随；
3. 优化了 ObsRadarVel的共线搜索性能；
4. 优化了 domainCheck 的性能；
5. 在 Application 加入了ERA5 替换的选项；
6. 修复了 Application 中不能正确输出增量的 bug。

## Update 2025-03-05
By Zilong Qin
1. 将一些僵尸程序以及静态文件转移到了/etc/historicalCode, 如有误删，可从此目录找回；
2. 将 Application 种的 bcgFillin 转移到了 ObsTools 中，变成通用工具模块；

## Update 2025-03-11
By Yali Wu

1. 当迭代步长很小时，FRCG算法中原负值检查部分更新的d0可能会过于大，最终导致程序崩溃。通过控制rng_b>10^-6次方修复此问题。
2. SolverFRCG中增加是否做负值检查（BoundaryCheck: Projected or inner）和是否做线搜索（LineSearch）的判断。
3. Application中控制多重外循环的参数k变得更为灵活，可通过YAML设置。
   

## Update 2025-03-04
By Zilong Qin， Yongjian Huang

1. 在 state 里加入了 fillPresWithHydrostatic， 用于在输出分析场之前使用静力平衡关系更新三维的气压场；
2. 在 Laplace 算子加入系数考虑 lat，lon不相等的权重修正；
3. 为雷达数据加入了数据筛选参数，现在不会同化离地高度 1km 以下的雷达数据，不同化分析时刻 900秒以外的数据，数据跳点 1，以加快处理速度；
4. 降低雷达数据的观测权重，以优化同化效果。
5. 在static 目录下加入标准 yaml 参考文件：
```bash
/static/Application/App_3DVar_TH-XY_CMA-GD_6km.yaml
```
1. era5 数据检验脚本加入到了，大家可用作调试参考
```bash
/etc/VerifyERA5
```
## Update 2025-02-11
By Zilong Qin， Guangqiu Ma

1. 新增了一个BFieldGauss.F90，用于调试高斯相关函数的B 矩阵作用效果，目前效率极低，还无法使用；
2. 将拉普拉斯的权重针对 lat，lon网格不均匀的情况，做了一个权重调整。
3. 优化了 uv2divg 的实现方式，现在以 corner 位置保证连续性守恒，修正跳点的问题；
4. 为 IOWRF 添加了一个读取外部地形输入的功能与接口。

## Update 2025-02-27
### 风廓线 数据同化
1. 风廓线QC
- 参考张宇博士的方法，采用FAST-MCD对风廓线数据进行质控，质控后同化效果有明显提升；
- 在MOTOR-DA B矩阵Laplace的条件下，确定了一组相对较优风廓线的B矩阵的配置参数；
- 在Application中的ObsList数组中加入profiler可开启风廓线同化功能；

### GNSS RO 数据同化
1. GNSS RO 参数设置及说明（根据当前天河业务版本区域及网格配置设置）
- 支持的卫星： COSMIC-2, SPIRE, METOP, 和 天目星座
- 多尺度同化设置: 根据GNSS RO观测的空间代表性设置，GNSS RO观测的空间空间分辨率为200km ~ 300 km，根据当前业务网格G08为6km，从最小网格的6 * (4**3)=384km
```fortran
Application.F90
...
      CASE ('gnssrefr')
        IF (sg%gLevel >= this%mgEnd - 3) THEN
          CALL this%gnssro%ObsThinning(XbRef, Yi(12), mpObs, .TRUE., .FALSE.)
        END IF 
      END SELECT
...
```
- 同化变量： 支持同化 温度，湿度和气压 三个变量， 通过配置YAML中IO.GNSSRO_adj_mode实现对同化变量的控制， 设置为1只同化温度，设置为2同化温度和湿度，设置为其他值同时同化三个变量，默认的情况下设置为2.
```yaml
IO：
   ...
   GNSSRO_adj_mode: 2
   ...
```
- 观测时间设置： 通过设置YAML中IO.set_all_gnss_to_fcst_time为TRUE，可以强制把GNSSRO观测.的时间设置为预报时刻。默认情况设置为TRUE。
 
 ```yaml
IO：
   ...
   set_all_gnss_to_fcst_time: TRUE
   ...
```



## Update 2025-02-06
By Yuanfu Xie

1. Added a new solver of Neuman BC boundary condition for a Poisson solver, GMRES.
2. Extended the edge stencil to 7 from 6 in order to address the mixed Poisson solver need.

## Update 2025-02-06
1. 分开了 HbyInc 方案和其他业务方案的代码，以便于后续的维护和开发，涉及的代码包括：
```
src/Application/Applications_HybInc.F90
src/MOTOR-DA/Minimization/JFunc/JFunc_HybInc.F90
src/MOTOR-DA/Minimization/MiniSolver/SolverFRCG_HybInc.F90
src/Application/App_CMA_GD_HybInc.F90
```
有关 HybInc 方案的相关调试和测试，通过以上代码实现。

## Update 2025-02-05
By Yuanfu Xie (yuanfu_xie@yahoo.com)

1. I deveoped and tested a new minimization interface routine, FRCGSolver.F90, interfacing with a unified variational new cost function and its gradient that intends to provide flexibility of background error covariances and controls. I can use it to reproduce the analysis from App\_CMA\_GD\_Opr.exe.
2. I created a new data type of C2M, designed to unify all control options such as vorticity/divergence, 3DVar or 4DVar. It is still under a 4DVar testing now.

## Update 2025-02-05
By Yongjian Huang
1. 添加了fortran的stdlib，其作用类似于C++的标准库，提供一些常用的基础功能：
- 实用程序（容器、字符串、文件、操作系统/环境集成、单元测试和断言、日志记录等）
- 算法（搜索和排序、合并等）
- 数学（线性代数、稀疏矩阵、特殊函数、快速傅里叶变换、随机数、统计、常微分方程、数值积分、优化等）

具体参考： 中文官方网站 [https://stdlib.fortran-lang.cn/index.html]；GITHUB主页 [https://github.com/fortran-lang/stdlib] 或 [https://opensimi.coding.net/p/MOTOR/d/STDLIB_Fortran/git]

主要修改：（1）在external中添加了fortran_stdlib的submodule；（2）在compileExt.sh中添加了编译的功能；（3）coding上的版本是经修改后的兼容f2008编程规范版本。

## Update 2025-02-04
By Zilong QIN (zilong.qin@gmail.com)

1. Merge the feature_HybridEnLoc and add_motor-dpp branch;
2. Test pass for general case.
One more parameter should be specified in the YAML file:
``` yaml
RunMode:
   ...
  Framework: FullState #FullState #Incremental
```

## Update 2025-02-04 
By Zilong QIN (zilong.qin@gmail.com)

### 功能改进与更新
1. 在代码格式化中加入`git diff`的判断，现在代码格式化近会对本次修改的代码进行格式化，而不是仓库内所有的文件；
2. BFieldLaplace_t调整：
   1. 在BFieldLaplace_t成员函数中加入了是否启用Neumann BC的选项，如果使用为true，则会在边界处使用Neumann边界条件；否则仅使用内点Laplace作为约束条件，调整边界；默认设置为true;
   2. 调整了scaleParaX， scaleParaY， scaleParaZ, scaleParaT随网格层数的调整变化，现在变化倍数与两层网格的分辨率变化成正比，更贴近算子离散格式的变化比例;
   3. 现在调整了scaleParaX， scaleParaY， scaleParaZ, scaleParaT使用归一化配置，加入了RelativeWeightJb2Jo作为扩散球整体大小，scaleParaX， scaleParaY， scaleParaZ, scaleParaT作为相对椭球半径；
   4. 加入了对谈空温度的随网格指数权重调整，以单独控制温度的影响半径等。
   ``` yaml
   BMat:
      ScaleParaX: 49    # Relative radius at longtitude direction
      ScaleParaY: 49    # Relative radius at latitude direction
      ScaleParaZ: 1.9   # Relative radius at vertical direction
      ScaleParaT: 0.1   # Relative radius at time direction
      RelativeWeightJb2Jo: 200.0    # Relative weight between Jb and Jo, radius of diffusion sphere
   ```
3. RField_t调整：
   1. 将所有观测数据，按照网格/观测数量比做了一个权重归一化，以配合B矩阵调整，实现在不同网格之间固定扩散半径的效果；
   2. 现在各变量之间的sigma调整相对独立，调整单个观测的sigma不会影响其他观测数据的收敛；
   3. 加入了一个多重网格相关的指数权重，通过调整权重可以影响观测的扩散情况：
   4. 探空水汽观测从1.0调整为0.25，温度权重从2.0调整为1.0；
4. SingleGrid_t调整：
   中加入了mgParent变量，以从sg中访问mg相关的变量，用体使用如下所示:
   ``` fortran
   SELECT TYPE (mg => sg%mgParent)
    TYPE is (MultiGrid_t)
      ! mg is the object of MultiGrid
    END SELECT
   END SELECT
   ```
5. State_t 调整：
   添加clearHalo函数，用于清除halo区域的数据, 以配合在伴随中使用；
6. C2O_t调整：
   将uv2w的是否启用改为yaml参数启动，不再根据控制变量检测结果启用；
7. oprRadarVel_t调整：
   1. 在执行伴随之前提前clearHalo；
   2. oprRadarVel_t不再依赖是否有wwnd输入变量；
8. 修改GrapesIO_t中的读取背景场格式的选项，对应：
``` yaml
  grapes_model_type: CMA-GD # CMA-GD / CMA-MESO / CMA-GD-V3
# CMA-GD-V3:  对接CMA-GD V3.1及以前版本
# CMA-GD: 对接CMA-GD V3.2及以后版本
# CMA-MESO: 对应CMA-MESO
```
1. 在ObsSound_t中加入了过滤机制，对于18km以上的qvapor观测数据设置为missing，以提高模式积分的稳定度；
2.  在ObsSet_t中加入了rmVar函数，用于各个多重网格上快速挑选观测数据，Application中原修改观测命名的筛选方式变为：
``` fortran
   ! 新写法
   CALL Yi(1)%rmVar('SYNOP_pcpa')
   CALL Yi(1)%rmVar('SYNOP_pcpa5min')

   ! 原写法
   ! DO i = LBOUND(Yi(1)%obsFields, 1), UBOUND(Yi(1)%ObsFields, 1)
   !   IF (Yi(1)%ObsFields(i)%Get_Name() == 'pcpa') THEN
   !      Yi(1)%ObsFields(i)%values = 0.0D0
   !   END IF
   !   IF (Yi(1)%ObsFields(i)%Get_Name() == 'pcpa5min') THEN
   !      Yi(1)%ObsFields(i)%values = 0.0D0
   !   END IF
   ! END DO
```
对应在Application以及在obsSound中使用修改观测命名的方式筛选数据的代码，均已更新修改；
1.  针对风场观测数据，例如雷达径向风、风廓线、地面站、以及探空风等直接或者间接观测，加入了针对uwnd 和wwnd是否开启B矩阵扩散的选项：
``` yaml
BMat:
   ...
   disableBmatForUV: TRUE
```
如果此项为TRUE，则不会对uwnd和wwnd进行B矩阵扩散，此选项一般配合Junc中Use_JcTerm_InCompres一起使用：
``` yaml
Minimization:
  Use_JcTerm_InCompres: TRUE
  Weight_Jc_InCompres: 100
  Type_Jc_InCompres: UV2Divg  
```
例如在Junc中使用了UV2Divg，此时可选择不对uwnd和wwnd进行B矩阵扩散，因为Jc已经包含边界扩散信息；
1.  针对地面水汽和气压的扩散尺度情况，仅在最细的三层网格加入地面温度、地面水汽和气压观测：
```fortran
   IF (sg%gLevel < this%mgEnd - 2) THEN
      CALL Yi(1)%rmVar('SYNOP_temp')
      CALL Yi(1)%rmVar('SYNOP_qvapor')
      CALL Yi(1)%rmVar('SYNOP_pres')
      CALL Yi(1)%rmVar('SYNOP_psl')
   END IF
```

### 框架简化
1. 删除了nbCoarst等相关语句，获取粗一层网格数，可使用mgParent上到上一层的num_icell_global获取；
2. 删除了以下文件（**commit c6f2e99f3097c1b08f27438652f96ecd956cc18b**）：
   ```
    src/MOTOR-DA/BMatrix/BField/BFieldLaplace_Alpha-0.1.3.F90
    src/MOTOR-DA/BMatrix/BField/BFieldLaplace_Dirichlet.F90
    src/MOTOR-DA/BMatrix/BField/BFieldLaplace_HalfLap.F90
    src/MOTOR-DA/BMatrix/BField/BFieldLaplace_HalfLapPlusDirichlet1.F90
    src/MOTOR-DA/BMatrix/BField/BFieldLaplace_NullBdy.F90
    src/MOTOR-DA/BMatrix/BField/BFieldLaplace_OneSideLap.F90
    src/MOTOR-DA/BMatrix/BField/BFieldLaplace_Orig.F90
    setEnv_cjl.sh
   ```


## Update 2025-01-21
By Yali Wu (wuyali@gbamwf.com)
1. 选择需要同化的卫星观测资料种类不再依赖于App_CMA_GD_Opr.F90，而是通过YAML参数来设置。举例如下。同时，Application.F90中卫星资料同化部分改为更加简洁的界面。
   ```
   platform_name: [fy4_1, fy3_4, fy3_4, fy3_5, fy3_5]
   inst_name: [agri, mwts2, mwhs2, mwts3, mwhs3]
   turnOn: [FALSE, TRUE, FALSE, FALSE, FALSE]
   ```
2. 增加src/MOTOR-DP(data pre-processing)工作目录，目前用于对卫星亮温资料的前处理，包括读写IO、质控、偏差订正。注意：同一子时间窗内的观测资料被处理到子时间窗的末端。
3. 新增增量法同化框架（控制变量为增量形式）。测试结果：增量法和全量法可得到接近的分析结果。
4. 新增FY3 MWTS/MWHS两种微波辐射计资料的同化功能（By Yali Wu and Yuan Tang)。
5. 新增RTTOV CLOUD的模拟功能，用于同化IR cloudy radiance。
6. Add dynamic cloud error model for cloudy radiance assimilation.
7. 新增水物质作为控制变量。
8. 新增水物质的scaling廓线用于控制变量转换。
9. 新增读入/写出qcqr.dat作为同化背景场/分析场的模块。
10. 新增FY4B功能。同时，大幅优化了同一卫星不同系列仪器的代码设计。需要在FY4-AGRI模块下指定satellite为fy4_1还是fy4_2.
11. 新增cloudy static BC模块。
12. 增加JFunc各项的归一化功能。归一化系数设置为当前网格的背景场和各种类观测的数目。
13. 水物质的scaling廓线目前在singleGrid中执行，需要在YAML中配合```CV_Transform```使用。
14. 新增小波分析稀释化方法（By Ting Shu, 2024）。
15. Add U.S. standard atmospheric profile for calculating standard weighting function.
16. Add tests for understanding RTTOV nonlinearity.

## Update 2025-01-13 
By Jiongming PANG

1. Add Applied EnLoc as Hybrid.
2. Added localization by using RE refered to GSI.
3. Added incremental analysis scheme, only available in HybInc, which refered to Hybrid method used in GSI.


## Update 2025-01-10
By Yongjian HUANG & Zilong QIN

1. 构造了存储下标的结构体，以及qsort算法，加速了ObsSurface中的数据查重过程，速度提高了几十倍，测试与原过程完全一致；

## Update 2024-11-27 
By Yongjian HUANG

1. Add the halo width in the parallel halo exchange. 
``` yaml
geometry:
  # topoSmooth: [3,2,2,2,2]
  haloWidth: 2
```

## Update (V XXX)
By Guanqiu MA, Zilong QIN (zilong.qin@gmail.com)

### 功能改进与更新
