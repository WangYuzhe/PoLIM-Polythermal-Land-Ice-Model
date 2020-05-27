# PoLIM (Polythermal Land Ice Model)

## Authors
Wang Yuzhe (State Key Laboratory of Cryospheric Sciences, Chinese Academy of Sciences, China; College of Resources and Environment, University of Chinese Academy of Sciences, China. wangyuzhe@ucas.ac.cn)

Zhang Tong (Fluid Dynamics and Solid Mechanics Group, Los Alamos National Laboratory, USA. zhgtong@gmail.com)

## Descriptions
PoLIM is a 2D flowband thermomechanical ice flow model. It uses Blatter-Pattyn higher-order approximations and describes the energy transportation using the enthalpy method. It is designed for modeling the dynamics of mountain glaciers, and can be used to model the polythermal structures of moutain glaciers.

PoLIM is validated by the ISMIP-HOM benchmark experiments and Kleiner's enthalpy benchmark experiments. PoLIM also implements the Schoof-Hewitt water transportation scheme in temperate ice and a cavity-sheet subglacial hydrology model.

## Features
* PoLIM implements an enthalpy-based thermal model which is particularly convenient to simulate the dynamics of polythermal glaciers.
* PoLIM includes a drainage model to simulate the water transport in the temperate ice layer driven by gravity.
* PoLIM includes a subglacial hydrology model to simulate the subglacial water pressure for the coupling with the basal sliding law.


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
