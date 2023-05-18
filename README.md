# AutoOptLib: MATLAB library for automatically designing metaheuristic optimization algorithms

[![GNU General Public License v3.0](https://img.shields.io/badge/license-GNU%20GPL--v3.0-green.svg)](https://github.com/qz89/AutoOpt/blob/main/LICENSE)


AutoOptLib is a MATLAB library for automatically designing metaheuristic optimization algorithms. It provides:

1. A rich library of design choices covering over 40 representative metaheuristic components ranging from neighborhood search, evolutionary search, to swarm optimizer, as well as components for particular problems, 
e.g., niching for large-scale problems, constraint violation for constrained problems, and multiple sampling for uncertain problems. The components can be used to design algorithms for various types of problems, 
including continuous, discrete, and permutation problems with/without constraints and uncertainties.

2. Flexibility to designing diverse algorithms. AutoOptLib can design algorithms with diverse structures in a single run. The flexibility greatly boosts the possibility of finding novel and efficient algorithms.

3. Fair benchmark of different design objectives and techniques. AutoOptLib contains different types of design objectives, e.g., solution quality, runtime, and anytime performance, for various problem-solving scenarios. It offers many representative
techniques for searching, evaluating, and identifying algorithms, e.g., iterative racing, intensification, and surrogate. These options ensure correct benchmark comparison and convenient use for different demands and interests.

4. Good accessibility to both researchers and practitioners. AutoOptLib is fully written in MATLAB, such that it can directly interface with abundant MATLAB simulation environments. Furthermore, AutoOptLib has a concise GUI. Users can input their problem, 
manage the algorithm design process, make experimental comparisons, and visualize results with simple one-clicks on the GUI. AutoOptLib also supports automatically tuning hyperparameters of a given algorithm, using and comparing with classic algorithms.

5. Easy extensibility. AutoOptLib ollows the open-closed principle. Users are convenient to implement their own design choices, objectives, and techniques based on the current sources, and add the implementations to the library by a uniform interface. 
Such that the library is easy to extend to serve various problem-solving scenarios and update to stay state-of-the-art.

Please read the [Introduction](https://arxiv.org/abs/2303.06536) for understanding key features of AutoOptLib.  Please read the [Documentation](https://github.com/qz89/AutoOpt/blob/main/Documentation.pdf) for user guide.

AutoOptLib is a free software. You can use, redistribute, and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version. 

AutoOptLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public 
License along with this library.

Please reference the paper below if using AutoOptLib in your publication:
@article{zhao2023autooptlib,
  title={AutoOptLib: A Library of Automatically Designing Metaheuristic 
         Optimization Algorithms in Matlab},
  author={Zhao, Qi and Yan, Bai and Hu, Taiwei and Chen, Xianglong and 
          Yang, Jian and Shi, Yuhui},
  journal={arXiv preprint 	arXiv:2303.06536},
  year={2023}
}

AutoOptLib is developed and actively maintained by the Swarm Intelligence Lab at the Department of Computer Science and Engineering, Southern University of Science and Technology.

Copyright (C) <2023>  <Swarm Intelligence Lab>

For any question, comment or suggestion, please contact Dr. Qi Zhao at <zhaoq@sustech.edu.cn>.

