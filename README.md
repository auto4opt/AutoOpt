## AutoOptLib (Automated Design of Metaheuristic Optimization Algorithms)
[![GNU General Public License v3.0](https://img.shields.io/badge/license-GNU%20GPL--v3.0-green.svg)](https://github.com/qz89/AutoOpt/blob/main/LICENSE)
![](https://img.shields.io/badge/Matlab-%3E%3D%202018a%20-blue.svg)
![](https://img.shields.io/badge/Windows-Pass-brightgreen.svg)
![](https://img.shields.io/badge/MacOS-Pass-brightgreen.svg)

## What does it provide?
AutoOptLib is a MATLAB library for automatically designing metaheuristic optimization algorithms. It could relieve manual algorithm design effort and gain enhanced performance beyond human-made algorithms. It provides:

*  **Rich library of design choices** - Over 40 representative metaheuristic components for designing algorithms for continuous, discrete, and permutation problems with/without constraints and uncertainties.
*  **Flexibility to designing diverse algorithms** - Design algorithms with diverse structures in a single run, enables great possibility to find novel and efficient algorithms.
*  **Fair benchmark of various design objectives and techniques** - Various design objectives, e.g., solution quality, runtime, and anytime performance. Different design techniques, e.g., racing, intensification, and surrogate.
*  **Good accessibility** - GUI for users to input problems, manage the algorithm design process, make experimental comparisons, and visualize results with simple one-clicks. 
*  **Easy extensibility** - Easily add new algorithmic components, objectives, and techniques by a uniform interface. 

Real applications of AutoOptLib are reported [here](https://arxiv.org/abs/2303.06536), which illustrates its efficiency and practicability.  

## How to use it?
Follow the steps below to use AutoOptLib:

1. Download source code and add it to MATLAB path.

2. Implement the targeted optimization problem.

3. Define the space for designing algorithms.

4. Run AutoOptLib by command or GUI. 

Refer to the [Documentation](https://github.com/qz89/AutoOpt/blob/main/Documentation.pdf) for user guide.

Refer to the [Introduction](https://arxiv.org/abs/2303.06536) for AutoOptLib's key ideas and features.  

## Copyright
AutoOptLib is developed and maintained by the Swarm Intelligence Lab at the Department of Computer Science and Engineering, Southern University of Science and Technology.

Copyright (C) <2023>  \<Swarm Intelligence Lab\>

Users can use, redistribute, and modify AutoOptLib under the terms of the GNU General Public version 3 or any later version. AutoOptLib is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose.

Users should reference the following papers if using AutoOptLib in your publication:
```
@article{zhao2022autoopt,
  title={AutoOpt: A General Framework for Automatically Designing Metaheuristic Optimization Algorithms with Diverse Structures},
  author={Zhao, Qi and Yan, Bai and Chen, Xianglong and Hu, Taiwei and Cheng, Shi and Shi, Yuhui},
  journal={arXiv preprint arXiv:2204.00998},
  year={2022}
}

@article{zhao2023autooptlib,
  title={AutoOptLib: A Matlab Library for Automatically Designing Metaheuristic Optimization Algorithms},
  author={Zhao, Qi and Yan, Bai and Hu, Taiwei and Chen, Xianglong and Yang, Jian and Shi, Yuhui},
  journal={arXiv preprint arXiv:2303.06536},
  year={2023}
}
```

## Contact
For any question, comment, or suggestion, please contact Dr. Qi Zhao, zhaoq@sustech.edu.cn.
