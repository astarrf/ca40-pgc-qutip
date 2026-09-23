# ⁴⁰Ca⁺ 偏振梯度冷却：QuTiP 仿真与阶段性结论

## 397 nm Joshi 四能级模型

[397 nm Joshi 模型说明](simulation/JOSHI_397NM.md) · [QuTiP 脚本](simulation/joshi_397nm.py)

这一脚本使用 S₁/₂↔P₁/₂ 四能级和完整的运动位置指数算符，默认参数为 690 kHz 轴向频率、初态 n̄=15、+30 MHz 蓝失谐；并提供 +210 MHz 对照模式。它与下文的 393 nm Li 有效模型是两个独立模型，尚未完成默认大截断的收敛计算。


基于 **Li et al., “Robust polarization gradient cooling of trapped ions”, New Journal of Physics 24, 043028 (2022)** 的 393 nm 路线，建立单离子、单运动模式的有效主方程模型。

> **研究状态：机制演示与参数探索。尚未定量复现 Li 2022，也没有预测实际装置的最终温度。**
> 数值收敛与物理模型完整性是两项独立要求；本仓库明确保留尚未解决的截断误差。

## 主要结论

1. **本次有限基底计算出现了净降温。** 在 285 kHz 模式、−180 MHz 失谐、50 kHz 双束频差、单束 Ω/2π=20 MHz 下，N=48 的实际初始平均声子数 5.971，在 150 μs 后降至 5.093。驱动强度与初态是演示选择，并非论文全部参数的复现。
2. **静态梯度的相位很重要。** 在相同名义热初态 n̄₀=3、N=32 下，50 μs 后相位 0 给出 n̄≈2.856，相位 π/2 给出 n̄≈3.888。相位 0 在 N=48 下复算为 2.890，短时间降温趋势仍存在。不能将这一结果推广为“静态梯度普遍更好”。
3. **不能认为 PGC 总会继续降温。** 某些初态和驱动组合出现升温，因此需要实际扫描初态、光强、失谐及相位。
4. **目前主要数值限制是 Fock 空间截断。** 较热初态的 N=40 与 N=48 轨迹最大差为 0.22431 个声子；弱驱动长时间轨迹也尚未充分收敛。上述终点不能当作高精度稳态冷却极限。

![固定初态的参数与相位对照](simulation/results/parameter_comparison.png)

完整数据、参数、收敛比较和解释见 **[结果报告](simulation/RESULTS.md)**；公式、近似和运行方法见 **[模型说明](simulation/README.md)**。

## 模型包含与省略的物理

包含两个 S₁/₂ Zeeman 基态、一个量子谐振子、偏振梯度光移、光抽运、自发辐射反冲、移动梯度和常数环境加热。P₃/₂ 激发态按远失谐、弱激发条件消去，其跃迁系数保留在有效算符中。

尚未显式包含 D 态和 854/866 nm 再泵浦、激发态 Zeeman 分辨的失谐、ASE、微运动以及三维模式耦合。理想闭合 S–P 循环并不等同于真实再泵浦过程。代码不适用于近共振 Doppler/PGC crossover 或几千声子的热捕获。

## 快速运行

本仓库不包含虚拟环境或缓存。推荐 Python 3.12；已运行版本为 QuTiP 5.3.1，依赖版本保存在 `simulation/requirements.txt`。

Windows PowerShell，在仓库根目录：

```powershell
python -m venv .venv
& '.\.venv\Scripts\python.exe' -m pip install -r simulation/requirements.txt
& '.\.venv\Scripts\python.exe' simulation/test_model.py
& '.\.venv\Scripts\python.exe' simulation/pgc_qutip.py --label my_run
```

macOS / Linux：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r simulation/requirements.txt
.venv/bin/python simulation/test_model.py
.venv/bin/python simulation/pgc_qutip.py --label my_run
```

重新执行完整比较集合，可运行 `simulation/run_suite.py`；顺序执行可能需要几十分钟。已有完整结果随仓库保存，不需要重跑才能阅读。相同输出名称会覆盖同名数据，请为新实验选择新的 `--label`。

## 文件导航

| 路径 | 内容 |
|---|---|
| [simulation/RESULTS.md](simulation/RESULTS.md) | 本次实际结果及误差限制 |
| [simulation/README.md](simulation/README.md) | 主方程、参数定义、适用边界 |
| [simulation/pgc_qutip.py](simulation/pgc_qutip.py) | QuTiP 主程序 |
| [simulation/test_model.py](simulation/test_model.py) | 7 项物理与数值测试 |
| [simulation/run_suite.py](simulation/run_suite.py) | 完整重现入口 |
| [simulation/results/](simulation/results/) | CSV、参数 JSON、最终声子分布及图表 |

初步 `pilot` 数据保留用于追溯，不作为主要结论。每组正式 JSON 保存实际初始声子数、初始截断尾部、最高四态布居、迹误差、最终最小本征值及软件版本。

## 验证状态

- 7 项小型物理/数值测试通过。
- 50 μs 基准算例中，反冲积分从 3 到 5 节点，n̄(t) 最大差约 1.28×10⁻⁵。
- 同一基准算例中，求解容差收紧 100 倍，最大差约 2.2×10⁻⁸。
- 所有检查仅支持对应工作点；不能消除更大声子截断或不完整物理模型的误差。

## 参考与来源

- [Li et al. (2022), DOI: 10.1088/1367-2630/ac6233](https://doi.org/10.1088/1367-2630/ac6233)
- [Reiter & Sørensen (2012), Effective operator formalism, PRA 85, 032111](https://doi.org/10.1103/PhysRevA.85.032111)
- [QuTiP master-equation documentation](https://qutip.readthedocs.io/en/stable/guide/dynamics/dynamics-master.html)

该模型与报告由 AI 辅助编写，并实际执行了所附计算；不是 Li 作者提供的软件。结果记录于 2026-09-22。本仓库不包含原论文 PDF、私人凭据或本地运行环境。
