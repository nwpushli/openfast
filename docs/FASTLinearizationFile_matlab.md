# 在 MATLAB 中调用 `FASTLinearizationFile`

`FASTLinearizationFile` 是 OpenFAST 仓库随附的 Python 类，位于 `reg_tests/lib/fast_linearization_file.py`，用于读取线性化生成的 `.lin` 文件。MATLAB 可以通过其内置的 Python 接口直接调用该类，从而在 MATLAB 脚本中获得 `A/B/C/D` 矩阵以及输入输出描述。

## 1. 准备 Python 环境
1. 确认 MATLAB 能够调用 Python (`pyenv` 查看版本)。
2. 将 OpenFAST 仓库根目录加入 Python 的搜索路径：
   ```matlab
   rootPath = 'C:/path/to/openfast';   % 修改为本地仓库路径
   if count(py.sys.path, rootPath) == 0
       insert(py.sys.path, int32(0), rootPath);
   end
   ```

## 2. 读取 `.lin` 文件
```matlab
mod = py.importlib.import_module('reg_tests.lib.fast_linearization_file');
lin = mod.FASTLinearizationFile('MyTurbine.1.lin');
```
`lin` 是一个 Python 字典对象，键名包括 `A`、`B`、`C`、`D`、`x`、`u`、`y` 等。

## 3. 将矩阵转换为 MATLAB 数组
MATLAB 通过 `double()` 即可把 `numpy.ndarray` 转换为 `double` 类型：
```matlab
A = double(lin{'A'});
B = double(lin{'B'});
C = double(lin{'C'});
D = double(lin{'D'});
```

## 4. 获取输入/输出描述
```matlab
uDesc = cellstr(lin.udescr());
yDesc = cellstr(lin.ydescr());
```

## 5. 与 `compute_hase_transfer` 联用示例
```matlab
[Hase, G, sys, omega] = compute_hase_transfer(A, B, C, D, ...
    'Inputs', {'HWindSpeed'}, 'Outputs', {'TwrBsFys'}, ...
    'InputDescriptions', uDesc, 'OutputDescriptions', yDesc, ...
    'FrequencySpan', [0.1 5]);
```

以上步骤即可在 MATLAB 中直接使用 `FASTLinearizationFile` 读取线性化文件，并继续完成传递函数、频率响应等分析。
