ANN (Artificial Neural Network) function for plumed
====================

This is plumed ANNB function (annfunc) module.  It implements `ANNB` class, which is a subclass of `Function` class.  `ANNB` class takes multi-dimensional arrays as inputs for a fully-connected feedforward neural network with specified neural network weights and generates corresponding outputs.  The `ANNB` outputs can be used as collective variables, inputs for other collective variables, or inputs for data analysis tools. Compared to the `ANN` class, the `ANNB` class supports batch normalization and `tanh` activation function.

## Installation

Enable compilation by adding the `--enable-modules=annbfunc` to the configure command.

## Usage

It is used in a similar way to the `annfunc` and [other plumed functions](https://www.plumed.org/doc-v2.5/user-doc/html/_function.html).  To define an `ANNB` function object, we need to define following keywords:

- `ARG` (string array): input variable names for the fully-connected feedforward neural network

- `NUM_LAYERS` (int): number of layers for the neural network

- `NUM_NODES` (int array): number of nodes in all layers of the neural network

- `ACTIVATIONS` (string array): types of activation functions of layers, currently we have implemented "Linear", "Tanh", "Circular" layers, it should be straightforward to add other types as well

- `WEIGHTS` (numbered keyword, double array): this is a numbered keyword, `WEIGHTS0` represents flattened weight array connecting layer 0 and layer 1, `WEIGHTS1` represents flattened weight array connecting layer 1 and layer 2, ...  An example is given in the next section.

- `BIASES` (numbered keyword, double array): this is a numbered keyword, BIASES0 represents bias array for layer 1, BIASES1 represents bias array for layer 2, ...

Assuming we have an `ANNB` function object named `ann`, we use `ann.node-0, ann.node-1, ...` to access component 0, 1, ... of its outputs (used as collective variables, inputs for other collective variables, or data analysis tools).

## Examples

Assume we have an ANNB with numbers of nodes being [2, 3, 1], and weights connecting layer 0 and 1 are

```
[[1,2],
[3,4],
[5,6]]
```

weights connecting layer 1 and 2 are

```
[[7,8,9]]
```

Bias for layer 1 and 2 are

```
[10, 11, 12]
```

and 

```
[13]
```

respectively.

All activation functions are `Tanh` with batch normalization applied to each hidden layer.

Then if input variables are `l_0_out_0, l_0_out_1`, the corresponding `ANN` function object can be defined using following plumed script: 

```
ann: ANNB ARG=l_0_out_0,l_0_out_1 NUM_LAYERS=3 NUM_NODES=2,3,1 ACTIVATIONS=BNTanh,BNTanh  WEIGHTS0=1,2,3,4,5,6 WEIGHTS1=7,8,9  BIASES0=10,11,12 BIASES1=13 GAMMAS0=,  BETAS0=, EXPECTATIONS0=, VARIANCES0=
```

This plumed script can be generated with function `Plumed_helper.get_ANNB_expression()` in [this](https://github.com/weiHelloWorld/plumed_helper/blob/master/plumed_helper.py) repository.  Following is the Python code using this function to generate the script above:

```Python
from plumed_helper import Plumed_helper
ANN_weights = [np.array([1,2,3,4,5,6]), np.array([7,8,9])]
ANN_bias = [np.array([10, 11, 12]), np.array([13])]
Plumed_helper.get_ANN_expression('ANN', node_num=[2, 3, 1], 
                                 ANN_weights=ANN_weights, ANN_bias=ANN_bias,
                                 activation_list=['BNTanh', 'BNTanh'])
```
Further detailed templates can be found in [examples](../../examples).

## Authors

1. Version 1: Wei Chen (UIUC, weichen9@illinois.edu) and Andrew Ferguson (University of Chicago, andrewferguson@uchicago.edu)
2. Version 2: Nicholas Herringher, Siva Dasetty, and Andrew Ferguson (U.Chicago)

## Copyright

See ./COPYRIGHT
