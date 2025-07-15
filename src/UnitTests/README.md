# MOTOR 单元测试目录

本目录包含 MOTOR 项目的所有单元测试文件，按模块组织。

## 目录结构

```text
UnitTests/
├── Application/        # 应用程序相关测试 (7 个文件)
├── MOTOR-DA/          # 数据同化模块测试 (51 个文件)
├── MOTOR-DP/          # 数据处理模块测试 (6 个文件)
├── MOTOR-PS/          # 物理求解模块测试 (43 个文件)
├── MOTOR-QC/          # 质量控制模块测试 (23 个文件)
├── MOTOR-Repository/  # 基础仓库模块测试 (6 个文件)
├── Template/          # 模板模块测试 (2 个文件)
├── Utilities/         # 工具函数测试 (28 个文件)
├── CMakeLists.txt     # CMake 构建配置
├── test_petsc.f90     # PETSc 测试文件
├── Test_Poisson.F90   # 泊松求解器测试
├── TestMixedSolver.F90 # 混合求解器测试
└── TestShallowWater.F90 # 浅水方程测试
```

## 测试文件统计

- **总计**: 166 个测试文件
- **MOTOR-DA**: 51 个文件 (数据同化相关)
- **MOTOR-PS**: 43 个文件 (物理求解相关)
- **Utilities**: 28 个文件 (工具函数相关)
- **MOTOR-QC**: 23 个文件 (质量控制相关)
- **Application**: 7 个文件 (应用程序相关)
- **MOTOR-DP**: 6 个文件 (数据处理相关)
- **MOTOR-Repository**: 6 个文件 (基础仓库相关)
- **Template**: 2 个文件 (模板相关)

## 使用说明

1. 所有测试文件已从原始模块目录迁移到此统一目录
2. 测试文件按功能模块分类组织，便于管理和维护
3. 每个子目录包含对应模块的所有单元测试
4. 使用 CMake 构建系统进行编译和运行

## 注意事项

- 测试文件的依赖关系可能需要调整 CMakeLists.txt 配置
- 运行测试前请确保相关模块已正确编译
- 某些测试可能需要特定的输入文件或环境配置

## 历史记录

- 2025-06-30: 将所有分散的单元测试文件迁移到统一的 UnitTests 目录
