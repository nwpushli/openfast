# 在 MATLAB 中使用 `FASTLinearizationFile`

仓库的 `docs/FASTLinearizationFile.m` 提供了 MATLAB 版本的 `FASTLinearizationFile` 类，用于读取 OpenFAST 线性化结果 `.lin` 文件。
下面给出在 MATLAB 中使用该类的典型步骤。

## 1. 将工具脚本加入路径
```matlab
rootPath = 'C:/path/to/openfast';  % 修改为本地仓库路径
addpath(fullfile(rootPath, 'docs'));
```

## 2. 读取 `.lin` 文件
```matlab
lin = FASTLinearizationFile(fullfile(rootPath, 'build', 'MyTurbine.1.lin'));
```
读取完成后，线性化矩阵和操作点信息以属性形式存储在对象中：

- `lin.A`, `lin.B`, `lin.C`, `lin.D`
- `lin.u`, `lin.y`, `lin.x`, `lin.xdot`
- `lin.udescr()`, `lin.ydescr()` 用于获取简写后的通道描述。

## 3. 在 MATLAB 中继续分析
结合仓库中的 `compute_hase_transfer.m` 可直接建立传递函数或计算频率响应：

```matlab
[Hase, G, sys, omega] = compute_hase_transfer(lin.A, lin.B, lin.C, lin.D, ...
    'Inputs', {'HWindSpeed'}, 'Outputs', {'TwrBsFys'}, ...
    'InputDescriptions', lin.udescr(), 'OutputDescriptions', lin.ydescr(), ...
    'FrequencySpan', [0.1 5]);
```

如需查看原始描述，可访问 `lin.u_info.Description`、`lin.y_info.Description` 等字段；
若 `.lin` 文件包含 `dUdu`、`dUdy` 等雅可比矩阵，类也会自动读取并保存在对应属性中。

通过上述步骤，即可在 MATLAB 中无缝读取和使用 OpenFAST 的线性化结果，完成塔底疲劳载荷估算、模态分析或控制设计等任务。
