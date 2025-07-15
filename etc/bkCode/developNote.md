# Develop note

### Tuesday, June 21, 2022

~~`alpha-0.1.1`~~ `develop`  修改yaml文件中, mgStart和mgEnd项目的位置，从Minimization转移到geometry, 所有代码中相关读取的部分都做了调整.
~~`alpha-0.1.1`~~ `develop`  改动geometry，geometry只生成mgStart到mgEnd的网格，不再生成其他层  #633 

### Wednesday, June 22, 2022

~~`alpha-0.1.1`~~ `develop`  修改为多重网格的执行层号由yaml文件中的mgStart与mgEnd来确定，不再直接使用代码中的设定值
~~`alpha-0.1.1`~~ `develop`  在obsSurface与obsSound中加入了带有判断条件的转换，如果控制变量中有`rho` 或者`rho_ctl`，则会将观测中的`pres`转换为`rho`
~~`alpha-0.1.1`~~ `develop`  在Application中加入 submodule dumpAna来输出分析之后的背景场文件

### Monday, June 27, 2022

~~`alpha-0.1.1`~~ `develop`  构建App_CMA_GD_Opr.F90，用于MOTOR-3DVar的业务化测试，此程序对应超算上部署在/public/home/simi/optest/3DVarVerification的同化测试系统
~~`alpha-0.1.1`~~ `develop`  在Application module中加入了雷达径向风同化的代码与配置
~~`alpha-0.1.1`~~ `develop`  修改OprRadarVel module中的bug，在生成径向风计算角度的辅助参数时，加入是否rwnd的检测，只对rwnd进行此计算
~~`alpha-0.1.1`~~ `develop`  优化application.F90 在并行计算时的时序，删除了一些barrier的使用，当对sub进程实施barrier时，其他进程恢复需要较长的时间
~~`alpha-0.1.1`~~ `develop`  在conversion.F90中加入了RH_to_qvapor，将相对湿度转换为qvapor
~~`alpha-0.1.1`~~ `develop`  在obsSurface加入了带有判断条件的转换，如果控制变量中有`qvapor` 则将观测变量中的`RH`，转换为`qvapor`，并在转换之前加入了数据有效性的检测（检测missing数据）
~~`alpha-0.1.1`~~ `develop`  修正gridGenLatLon.F90中chld_indx生成错误的bug，现在MOTOR框架可以使用非对称网格结构了，比如，使用641 X 769（640 X 768），多重网格会按照实际实际格点划分到不可被2整除
~~`alpha-0.1.1`~~ `develop`  在IOWRF中加入了t2m 和 pres的地面数据数据输入
~~`alpha-0.1.1`~~ `develop`  为适应业务化配置，在yaml配置文件中加入了RunMode-Mode条目和RunMode-Task条目，Mode设置为Debug/UnitTest/Alpha/Beta/Release，Task条目设置为SfcAna/3DAna/3DVar/4DVar，当时Mode设置为Debug时，Application中会输出背景场、O-A、ThinningObs和最终的分析场，如果设置为其他配置，则只输出最终的分析场
~~`alpha-0.1.1`~~ `develop`  改进了State2NC的输出逻辑，用mo_netcdf库替换了使用原生库的配置，将long，lat调整为了变量的前两个维度，用于ncview绘图方便
~~`alpha-0.1.1`~~ `develop`  在SolverFRCG中的有约束条件中加入了pcpa变量大于0的约束条件
~~`alpha-0.1.1`~~ `develop`  注释了ObsBase中差值部分的权重系数检测：IF (useBkg) gaussian = gaussian * EXP(-(bkgdAtObs-forward(k,j,i))**2)
~~`alpha-0.1.1`~~ `develop`  改进了B矩阵边界上的导数守恒处理，导数方向全部改完了由内向外方向

### Tuesday, July 5, 2022

`alpha-0.1.1`  调整了Geometry中scatter/gather，scatterv/gatherv等语句的时序，发现当两个以上阻塞进程使用时，mpi容易发生网络错误，可能是openmpi对目前石岩超算的网络构架支持性不好，时序调整后目前正常
`alpha-0.1.1`  改进了Minimization中迭代步长检测的逻辑，目前一维搜索不会一直收敛到步长的局部最优，在获取到下降步长之后立刻返回，同时加入了一个J与J梯度的重算检测，每当有步长搜索与约束条件收缩时，会回算JJgrad
`alpha-0.1.1`  nc输出逻辑不再支持重复变量名的输出，使用Obs2StateBaseTypeName来替换Obs2StateBaseName输出thinning之后的观测
`alpha-0.1.1`  在mpddBase中加入了AllRedLOR函数，用于求或逻辑的并行规约
`alpha-0.1.1`  改进在ObsSurface中的unique check逻辑，由字符串检查优化为时间数值检查，加速约8s
`alpha-0.1.1`  改进求解器的边界检查机制，让边界检查在第一次评估JJGrad之前，确保在有约束情况下减少一次计算迭代的时间

### Wednesday, July 6, 2022

`alpha-0.1.1`  改进yaml_read的运行逻辑，现在每次读取配置，只在第一次读取文件，后面再次使用不再重复读取文件
`alpha-0.1.1`  加入了PostProc-Export2SelDomain模块，可以在配置文件中设定，输出给定区域、给定分辨率的分析场
`alpha-0.1.1`  修改了IOGrapes的pi、th、pip和thp的处理方式，由grapes直接读入的th为全量的theta，pi相同，现在读如th、pi后，计算thBase和piBase，最后拿分析场减去两者，得到最终的pip和thp输出
`alpha-0.1.1`  优化了aggr和scatter并行函数，两者目前不再会频繁的获取分配缓存的大小，同时优化了并行函数的时序
`alpha-0.1.1`  将IOGrapes在并行插值之前的时间插值改为单精度计算，以减小计算量

### Wednesday, July 13, 2022
`alpha-0.1.1`  将IOGrapes中加入了q值的负值检测修正，所有q中的负值读入之后否会被设置为1e-10

### Thursday, July 14, 2022
`alpha` 加入新的functions在Applications.F90实现新的thinning策略。但是在对比中出现segmentation fault。 出现在JFunc.F90的程序里。我尝试了两天没有找出segmentation fault的原因。我暂时推送这个版到coding上请子龙帮助检查一下。在第131行和137行，我分别做了debugging打印，131行输出正常，但是第137行没有输出。最终的结论是第135行的减法出现segmentation。

经过几天的debugging，找到原因：在生成ObsSet数据结构时，没有调用其对应的constructor, ObsSet_t。还没有解释为什么前面实验这个问题并没有出现，在测试新的thinning策略时，发现。解决方案：在调用ObsConcat_s生成新合并数据Y之前，先叫Y=ObsSet_t。经过反复试验，系统没有在出现问题。于7月18日，将修改和清理debugging信息版推送到coding。
### Thursday, July 14, 20
`alpha-0.1.1`  在conversion_m中加入露点温度dew point 到qvapor的转换，并实施到ObsSound观测算子中
`alpha-0.1.1`  在雷达算子子中加入unambigious range变量的读入，同时，在筛选数据的时候加入了对于不模糊距离的筛选检测
`alpha-0.1.1`  注释了将地面气压在输入时插值到最底层的语句，ps变量在分析时并不使用，地面气压目前使用theta和pi的外插来填充
`alpha-0.1.1`  优化了雷达数据读取中的经纬度与高度计算函数，这部分是雷达数据输入时主要的计算热点，改进之后速度比之前快40%
`alpha-0.1.1`  将同一个站点的雷达数据做了数据的合并，按照之前的计算逻辑，每一个单独的雷达数据文件作为一个独立的观测来约束，由于雷达的站点位置不变，同一站点不同时间的数据可以合并为同一项，thinning的径向风数据更加合理，速度更快

### Saturday, July 16, 2022
`alpha-0.1.1`  更新GrapesIO中顶层和底层的插值，现在使用突出边界的外插来填补，之前使用边界的外插
`alpha-0.1.1`  修正了GrapesIO中对于unstaggerred grid，顶层，底层高度场外插的bug

## Wednesday, July 20, 2022
`alpha` 与alpha-0.1.1合并。
`feature_geosbal` 编译出现问题：GeosBal.F90。

### Monday, July 18, 2022
`alpha-0.1.1`  改进了ObsSurface中的查重的效率，在GeoTools中加入elemental函数，比之前速度提高2倍
`alpha-0.1.1`  在IOGrapes和GrapsIO module中加入了对于pres作为控制变量的读取和导出支持
`alpha-0.1.1`  发现地面站数据中大部分的站点高度都是NaN，这里把所有地面站数据的高度都设置为了2e-7，使之强制置于地面上

### Wednesday, July 20, 2022
`alpha-0.1.1`  PI输出时加入了静力平衡计算的pi输出，使pi在输出时只参考Ps的数值，不再参考其他三位场变量
`alpha-0.1.1`  在Laplace中加入了条件判断，对于pres变量，只在水平和时间维度上扩散，不再在垂直维度上扩散
`alpha-0.1.1`  简化了IOGrapes，pi和theta值现在读入MOTOR-DA网格之后再计算T、p和rho，不再在输入之前计算



#Monday, August 1, 2022

`alpha` 增添yaml选项，允许用户选择同化风向风速或风场分量uv。

`alpha` 改进surface Ingest，根据yaml信息只读入相关变量，替换原来hard coded地面变量。

`alpha` 修改State.F90 getVarIdx函数功能，允许风场采用风向和风速作为观测值。
### Wednesday, July 27, 2022
`alpha-0.1.1`  在IOWRF中加入了空文件的检测
`alpha-0.1.1`  在SurfaceObs中加入了国家站的筛选，把读取风场变量改为了10米平均风
`alpha-0.1.1`  合并谢老师的 `alpha` 分支，修复了在Release模式下，多重网格插值不启动的bug，地面站站点检测中的CYCLE改为了EXIT，把从静态目录读取solosite文件改为了从配置文件中获取目录文件
`alpha-0.1.1`  修改BFieldLaplace，现在在垂直方向上在ScalePara小于50时，将使用无分母的Laplace算子来计算垂直方向的算子
`alpha-0.1.1`  修改BFieldLaplace的伴随，在垂直方向上对于地面层的计算增量为0，使地面增量可以向上扩散，上部增量不会对地面扩散

### Monday, August 1, 2022
`alpha-0.1.1`  在IOGrapes中加入了对非均匀网格的插值，可以对任意给定层的的网格分层实施垂直多重网格
`alpha-0.1.1`  修改了MultiGrid中prolongation中插值的系数，针对非均匀网格，把0.5的插值系数修正为根据网格实际分层来计算系数
`alpha-0.1.1`  修改了多重网格的分层数限制，让垂直分层最少不少于5层，能够计算垂直方向上的Laplace算子

### Tuesday, August 2, 2022
`alpha-0.1.1`  发现q水汽廓线在同化时出现异常，修改回了有关边界部分的拉普拉斯算子的改进
`alpha-0.1.1`  为了平衡各算子的量级，现在暂时将IOGrapes，ObsSound，ObsSurface中pres的单位由Pa缩小了100倍，变为hPa，qvapor的单位由kg/kg变为k/kg，增大了1000倍

### Sunday, August 7, 2022
`alpha-0.1.1`  在obsBase中domainCheck之后加入了对于地面观测变量的修正，将气压按照模式面高度差换算到模式面上，温度按照0.65度的温度递减率换算到模式面上
`alpha-0.1.1`  调整application，使雷达径向风只在最细的三层网格上使用
`alpha-0.1.1`  调整RadarRAW，使雷达数据只使用9000以下
`alpha-0.1.1`  调整BFieldLaplace，现在垂直方向上，qvapor的扩散系数比其他变量大1倍
`alpha-0.1.1`  调整RField，使雷达径向风观测的sigma值设置为4，量级比其他变量小4倍
`alpha-0.1.1`  调整了多重网格的分层策略，现在最细网格不少于33层

### Wednesday, August 10, 2022
`alpha-0.1.1`  Laplace系数里水汽的系数从2倍改回为1，使其与其他变量的扩散保持一致
`alpha-0.1.1`  调整RField，使探空观测的sigma值设置为2，量级比其他变量小2倍
`alpha-0.1.1`  温度的边界层系数修正从0.65改为0.35，根据查阅资料，地面边界层的温度下降系数为0.3-0.4
`alpha-0.1.1`  探空现在修改为水平位移为0，所有数据都放置在初始时刻
`alpha-0.1.1`  地面站观测现在不再筛选只使用国家基本站，使用所有基本站点
`alpha-0.1.1`  垂直多重网格对于300层以下网格不再切分，垂直方向上暂时不使用多重网格
`alpha-0.1.1`  IOGrapes的垂直分层曲线现在参考unstaggerred grid

### Tuesday, August 16, 2022
`alpha-0.1.1`  IOGrapes里现在使用增量的方式来插值回到Grape文件，使结果更加合理
`alpha-0.1.1`  现在在输出气压场使用静力平衡插值时，使用虚温Tv来计算静力平衡Tv=T(1+0.608*qv)
`alpha-0.1.1`  探空数据中气压改为pres_d，使之不进入同化，探空数据改为了以高度来定标
`alpha-0.1.1`  添加了地面高度与模式面高度差异超过2000米，则不再修正气压和水汽

### Thursday, August 18, 2022
`alpha-0.1.1`  温度的递减率改会为0.65，气压根据温度的递减率进行修正
`alpha-0.1.1`  探空的气压观测数据不再进入同化分析

### Monday, August 22, 2022
`alpha-0.1.1`  加入相控阵雷达的读取程序，现在程序可以读取广东地区格式的相控阵雷达数据，程序会根据雷达数据类型选择不同的读取程序
`alpha-0.1.1`  在application中，加入了innovation的输出，用于查看同化分析扩散的效果

### Tuesday, August 30, 2022
`alpha-0.1.1`  MOTOR-DA现在已经接入了相控阵雷达数据格式，接入深圳两部雷达测试效果比较微小；
`alpha-0.1.1`  针对Laplace扩散不合理的问题，调整了每一层多重网格的Laplace系数，让每两层的系数增大1.7倍，借此减小粗网格的系数，使之更靠近观测，同时增大细网格的系数，使之更加平滑，使结果更加合理；
`alpha-0.1.1`  对于背景场插值输入时，有时会有一些点插值异常的问题，单独测试了slint的package，把编译选项改为default_real8之后，插值异常的问题解决，更新了slint的库文件；
`alpha-0.1.1`  在SolverFRCG中加入了iter的输入选择，如果输入的迭代次数大于给定的迭代次数，则使用给定的迭代次数。
`alpha-0.1.1`  径向风数据现在只在最细的两层网格同化

`alpha` Merged with alpha-0.1.1 and added a function of getVarIdx to ObsBase.F90. Fixed variable unit changes with checks to make sure the variables exist.

### Tuesday, August 16, 2022

`alpha` 新功能：在观测稀疏地方增加背景场作为观测，弥补观测缺失造成的分析问题。
`alpha` 在计算ObsSet内积时增添debugging信息，以便检查minimization对不同观测数据的收敛情况

### Friday, September 16, 2022

`alpha` passed unit test of Geostrophic balance with terrain following coordinate. As the terrain following coordinate test uses the background fields, this background is turned off in the automatic unit test.
`alpha` Fixed a few unit tests of MOTOR-DA. So now, all unit tests are passed except the satellite related tests.
`alpha` Fixed ObsSurface.F90, ObsSound.F90 and ObsVwpw.F90 hard coded variable indices.

### Saturday, October 8, 2022
`alpha-0.1.1` 筛除雷达径向风100m以下的数据
`alpha-0.1.1` 修正地面站观测中气压修正程序中的错误

### Friday, October 28, 2022
`alpha-0.1.2` 增加在MOTOR代码在天河超算上的适配，规避结构体数组拼接的使用范式

### Sunday, October 30, 20
10月份针对框架的诸多修改，暂时没有同步到更新日志中
`alpha-0.1.2` 增加了 Export2HASCoordInSelDomain 的 module， 用于将MOTOR-DA的分析场在输出时输出到离地高度坐标系的三维分析场上
`alpha-0.1.2` 在输出文件中增加了地形高度的输出
`alpha-0.1.2` 移除modeNO，直接使用任务配置字符串来判断任务属性

### Sunday, October 30, 2022
`alpha-0.1.2` IOGRAPES.F90: 1. printing *more info on saturated qvapor;* 2. *changing gravity acceleration to parameter;* 3. *modified the system call of cp inputfile to output, from -Lvr to -Lv as r requires other options*
`alpha-0.1.2` Applications.F90: cleaned up a little and modified the background fill-in code;
`alpha-0.1.2` parameters.F90: changed the gravitational acceleration parameter to 9.80665D0

### Wednesday, November 9, 2022
`alpha-0.1.2` 将ObsSurface融合谢老师修改的根据实际变量读取地面观测的版本；
`alpha-0.1.2` 在IOWRF中加入10分钟平均风速和2分钟平均风速风向的读入，取值为使用瞬时风速风向值；
`alpha-0.1.2` 在求解器、网格选取等位置加入对于pcpa5min的处理
`alpha-0.1.2` 在UV2W中加入网格层数筛选，只在网格大于3层以上的网格上实施
`alpha-0.1.2` 在OprRadarVel和RadarRAW中讲使用的地球几何半径设置为等效地球半径 8500km
`alpha-0.1.2` ObsRadarRef中设置为读取reflectivity
`alpha-0.1.2` 在YAMLRead中加入了针对文件名的校验，当文件名变化时，重新读取YAML文件

### Thursday, December 1, 2022
`alpha` 将张华老师的云导风模块merged到alpha版本，debugging，转换观测时间格式，增添分析时间窗口检验。

### Sunday, December 4, 2022
`alpha` 统一了MOTOR-QC里missing和invalid，都采用了在param.F90的常数。

### Monday, December 19, 2022
`alpha` Merged with alpha-0.1.2. Added features of cloud drift wind, unified missing and invalid_value in obsBase.F90 and formatted some print statements for more clear information for debugging.

### Monday, March 6, 2023
`alpha` Merged with alpha-0.1.3. Starting porting Z-grid model codes to alpha branch. Adding Z-grid model arrays to singleGrid.F90.