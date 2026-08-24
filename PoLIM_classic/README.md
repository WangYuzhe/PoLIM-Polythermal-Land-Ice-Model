# PoLIM — Polythermal Land Ice Model

PoLIM is an open-source, two-dimensional higher-order flowband model for the
dynamics and thermodynamics of mountain glaciers, with particular support for
polythermal glaciers. The momentum balance is based on the Blatter–Pattyn
approximation, while the energy balance is formulated in terms of enthalpy.

PoLIM was developed by **Yuzhe Wang and Tong Zhang**, with the scientific work
describing the model contributed by Yuzhe Wang, Tong Zhang, Cunde Xiao and Jiawen Ren.

## Main capabilities

- Two-dimensional flowline and flowband simulations
- Higher-order ice-flow dynamics based on the Blatter–Pattyn approximation
- Isothermal and polythermal configurations
- Standard and modified enthalpy-gradient thermodynamic schemes (SEGM and
  MEGM)
- Temperature- and water-content-dependent ice rheology
- Several basal boundary and sliding-law options
- Ice-thickness evolution with adaptive CFL substepping
- Gradient, temperature-index, and positive-degree-day surface mass-balance
  options
- Optional coupling to a cavity-sheet-type subglacial hydrology model
- Idealized and benchmark glacier geometries

Earlier versions of PoLIM were evaluated against standard experiments,
including ISMIP-HOM, enthalpy benchmark, and SHMIP experiments. Users should
perform validation appropriate to their own configuration, modifications, and
scientific application.

## Getting started

PoLIM is written in MATLAB. Clone or download the repository, start MATLAB in
the repository root, and run one of the two main scripts:

```matlab
main_isoT_forward   % isothermal example
main_polyT          % polythermal thermomechanical example
```

The distributed examples use the Arolla geometry in
`geo_inputs/geo_arolla.mat`. Model constants and most numerical and physical
options are defined in `params.m`; experiment-specific choices, time settings,
geometry, and surface mass-balance forcing are set near the beginning of each
main script.

Before using PoLIM for scientific production, users should inspect these
settings, provide suitable geometry and forcing, test numerical convergence,
and validate the configuration for the intended application.

## Repository structure

| Path | Purpose |
| --- | --- |
| `main_isoT_forward.m` | Isothermal model driver |
| `main_polyT.m` | Polythermal thermomechanical model driver |
| `params.m` | Physical constants and model options |
| `solver_u_core.m` | Higher-order velocity-system assembly |
| `solver_enthalpy_*.m` | Enthalpy solvers |
| `solver_subHydro.m` | Subglacial hydrology solver (beta) |
| `run_smb.m`, `calc_smb_*.m` | Surface mass-balance schemes and driver |
| `geo_inputs/` | Example and benchmark geometries |
| `mini_tools/` | Plotting and geometry utilities |

## Citation and responsible use

PoLIM is shared openly to support reproducible research and further scientific
development. Users are welcome and encouraged to use, test, modify, and extend
the model. In return, proper academic attribution is expected.

If PoLIM, its source code, or a substantial part of its methods or numerical
implementation contributes to research, software, a thesis, a report, a
presentation, or another scholarly output, please:

1. **Cite the PoLIM model paper** listed below.
2. **Acknowledge PoLIM and its developers** in the text or acknowledgements,
   and identify the version or commit used when possible.
3. **Describe material modifications clearly.** If PoLIM is adapted or
   extended, state that the resulting software is based on, derived from, or
   extends PoLIM, and distinguish the original components from the new work.
4. **Do not present a modified, renamed, or extended version of PoLIM as a
   wholly independently developed tool.** A new interface, application,
   calibration, module, or collection of modifications does not remove the
   obligation to credit the model and code on which it is substantially based.
5. **Preserve relevant copyright, authorship, citation, and license notices**
   when redistributing source code or derivative software, and comply with the
   license accompanying the version used.

Citation is part of responsible scientific practice and is distinct from any
permissions granted by a software license. If you are unsure whether your use
requires attribution, please err on the side of transparency and cite PoLIM.

A suitable acknowledgement is:

> This work used PoLIM (Polythermal Land Ice Model), developed by Yuzhe Wang
> and Tong Zhang. We used version [version/commit] and [briefly describe any
> modifications, if applicable].

### Recommended citation

Wang, Y., Zhang, T., Xiao, C., Ren, J., & Wang, Y. (2020). A two-dimensional,
higher-order, enthalpy-based thermomechanical ice flow model for mountain
glaciers and its benchmark experiments. *Computers & Geosciences, 141*,
104526. <https://doi.org/10.1016/j.cageo.2020.104526>

```bibtex
@article{Wang2020PoLIM,
  author  = {Wang, Yuzhe and Zhang, Tong and Xiao, Cunde and Ren, Jiawen and Wang, Yanfen},
  title   = {A two-dimensional, higher-order, enthalpy-based thermomechanical ice flow model for mountain glaciers and its benchmark experiments},
  journal = {Computers \& Geosciences},
  year    = {2020},
  volume  = {141},
  pages   = {104526},
  doi     = {10.1016/j.cageo.2020.104526}
}
```

## 中文说明：引用、致谢与负责任使用

PoLIM 是面向冰川动力学与热力学研究的开源科学模型。我们欢迎并鼓励同行使用、
测试、修改和进一步开发 PoLIM。开放源代码旨在促进科学交流、结果复现和模型共同
发展；与此同时，规范引用、如实说明代码来源并尊重原开发者的工作，是使用开源
科研软件应遵守的基本学术规范。

如果您的论文、学位论文、研究报告、会议报告、软件或其他成果使用了 PoLIM、
PoLIM 源代码，或实质性采用了其方法、数值方案和程序实现，请务必：

1. **引用上述 PoLIM 模型论文**；
2. **在正文或致谢中明确说明使用了 PoLIM，并致谢模型开发者**，同时尽可能注明
   所使用的版本号或 Git commit；
3. **如实、清楚地说明所作修改或扩展**，明确区分 PoLIM 原有内容与您新增的工作；
4. **不得将基于 PoLIM 修改、重命名或扩展而成的程序表述为完全独立开发的工具**。
   增加界面、改变应用对象、重新率定参数、添加模块或进行其他修改，均不能取代对原模型、源代码及开发者的恰当引用与说明；
5. 在传播源代码或衍生软件时，**保留相关的版权、作者、引用和许可信息**，并遵守所使用版本附带的软件许可条款。

软件许可所授予的使用或修改权限，并不等同于免除学术引用义务。若不确定是否需要引用，建议遵循公开透明的原则，主动引用 PoLIM 并说明其对相关工作的贡献。尊重开源科研软件的知识贡献，是维护学术诚信和推动科学共同体健康发展的重要基础。

中文致谢示例：

> 本研究使用了由王玉哲（Yuzhe Wang）和张通（Tong Zhang）开发的多温型陆地冰模型 PoLIM（Polythermal Land Ice Model）。本研究使用的版本为[版本号/commit]，并在此基础上进行了[简要说明所作修改；如无修改可删除此项]。

## Contributing

Bug reports, benchmark results, documentation improvements, and scientifically
well-documented model extensions are welcome. When proposing a change, please
explain its scientific motivation, assumptions, and validation, and include a
minimal reproducible example when practical.

## License

PoLIM_classic is made available under the
[PolyForm Noncommercial License 1.0.0](LICENSE). Noncommercial research, use,
modification, and redistribution are permitted subject to the license terms.
Commercial use is not licensed and requires separate written permission from
the applicable copyright holder or holders. Redistributions must also preserve
the required attribution notice in [NOTICE](NOTICE).

## Disclaimer

PoLIM is research software. Its users are responsible for evaluating numerical
stability, parameter choices, input data, uncertainty, and fitness for their
intended application. Results should be interpreted with appropriate scientific
judgement.
