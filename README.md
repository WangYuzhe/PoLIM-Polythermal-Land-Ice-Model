# PoLIM (Polythermal Land Ice Model)

## Authors
Wang Yuzhe (State Key Laboratory of Cryospheric Sciences, Chinese Academy of Sciences, China; College of Resources and Environment, University of Chinese Academy of Sciences, China. wangyuzhe@ucas.ac.cn)

Zhang Tong (Fluid Dynamics and Solid Mechanics Group, Los Alamos National Laboratory, USA. zhgtong@gmail.com)

## Descriptions
PoLIM is a 2D first-order flowband thermomechanical ice flow model. It is designed for modeling the dynamics and thermodynamics of mountain glaciers particularly for polythermal glaciers.

PoLIM has been verified by standard benchmark problems, including the ISMIP-HOM experiments, the enthalpy benchmark experiments, and the SHMIP experiments. PoLIM shows good performances and agrees well with these benchmark results, indicating its robust capability of simulating the thermomechanical behaviors of glaciers.

## Features
* PoLIM simplifies the momentum balance using the Blatter–Pattyn approximation.
* PoLIM implements an enthalpy formulation of the energy balance equation.
* PoLIM includes a scheme for gravity-driven drainage of water in temperate ice.
* PoLIM couples a cavity-sheet type subglacial hydrology model to ice dynamics.


## Citation
@Article{WangYuzhe2020,
  author   = {Yuzhe Wang and Tong Zhang and Cunde Xiao and Jiawen Ren and Yanfen Wang},
  journal  = {Computers & Geosciences},
  title    = {A two-dimensional, higher-order, enthalpy-based thermomechanical ice flow model for mountain glaciers and its benchmark experiments},
  year     = {2020},
  pages    = {104526},
  volume   = {141},
  doi      = {https://doi.org/10.1016/j.cageo.2020.104526},
  issn     = {0098-3004},
  keywords = {Polythermal glacier, Ice flow model, Thermomechanical coupling, Mountain glacier dynamics},
  url      = {http://www.sciencedirect.com/science/article/pii/S0098300419311458},
}
