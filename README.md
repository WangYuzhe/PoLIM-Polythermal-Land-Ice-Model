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
A paper introducing PoLIM has been submitted on Computers & Geosciences.

If you used PoLIM for your work, you could cite the following papers:

@Article{ZhangTong2013,
  author  = {Zhang, Tong and Xiao, Cunde and Colgan, William and Qin, Xiang and Du, Wentao and Sun, Weijun and Liu, Yushuo and Ding, Minghu},
  title   = {Observed and modelled ice temperature and velocity along the main flowline of {E}ast {R}ongbuk {G}lacier, {Q}omolangma ({M}ount {E}verest), {H}imalaya},
  journal = {Journal of Glaciology},
  year    = {2013},
  volume  = {59},
  number  = {215},
  pages   = {438-448},
  doi     = {10.3189/2013JoG12J202}
}

@Article{WangYuzhe2018_tc,
  author  = {Wang, Yuzhe. and Zhang, Tong. and Ren, Jiawen. and Qin, Xiang. and Liu, Yushuo. and Sun, Weijun. and Chen, Jizu. and Ding, Minghu. and Du, Wentao. and Qin, Dahe.},
  title   = {An investigation of the thermomechanical features of Laohugou Glacier No. 12 on Qilian Shan, western China, using a two-dimensional first-order flow-band ice flow model},
  journal = {The Cryosphere},
  year    = {2018},
  volume  = {12},
  number  = {3},
  pages   = {851--866},
  doi     = {10.5194/tc-12-851-2018}
}
