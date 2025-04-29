# Interface of AutoOptLib with Python problem files

#----------------------------Copyright-------------------------------------
# Copyright (C) <2025>  <Swarm Intelligence Lab>

# AutoOptLib is a free software. You can use, redistribute, and/or modify
# it under the terms of the GNU General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or any later 
# version. 
#--------------------------------------------------------------------------
def get_type():
    # get problem type
    # return type = ['continuous'/'discrete'/'permutation', 'static'/'sequential', 'certain'/'uncertain']
    # e.g.,
    type = ['continuous', 'static', 'certain']  # a static, continuous problem without uncertainty
    return type


def get_bound():
    # get solution space boundary
    # shape: [1, D], where D is the dimensionality of solution space, type: 'list'
    # e.g.,
    lower = [0, 0, 0, 0, 0]
    upper = [1, 1, 1, 1, 1]  # a 5D solution space
    return lower, upper


def evaluate(Decs, instanceInd):
    # evaluate solutions
    # 'Decs' is the solutions fetched from Matlab, shape: [N, D], where N and D are the number of solutions and the
    # dimensionality of a solution, respectively.
    # 'instanceInd' is the index of problem instance

    obj, con, acc = your_evaluate_method(Decs, instanceInd)
    # your_evaluate_method contains your code for evaluating solutions on the current problem instance
    # 'obj': solutions' objective values, shape: [N, 1] for single-objective optimization, type: 'list'
    # 'con': solutions' constraint violation values, shape: [N, 1], type: 'list'
    # 'acc': accessory data, shape: [N, P], one row for one solution' accessory data, type: 'list'
    # replace 'con' and 'acc' with '_' if not applicable
    return obj, con, acc
