/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MIT License

Copyright (c) 2019 Wei Chen and Andrew Ferguson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "function/Function.h"
#include "function/ActionRegister.h"
#include "cassert"

#include <string>
#include <cmath>
#include <iostream>
// #include <stdio.h>

using namespace std;

// #define DEBUG
// #define DEBUG_LAYEROUTPUT
// #define DEBUG_FINALDERIVATIVES
// #define DEBUG_LAYERDERIVATIVES

namespace PLMD {
namespace function {
namespace annfunc {

//+PLUMEDOC ANNMOD_Function ANN
/*
Calculates the ANN-function.

This module implements ANN class, which is a subclass of Function class.
ANN class takes multi-dimensional arrays as inputs for a fully-connected feedforward neural network
with specified neural network weights and generates corresponding outputs.
The ANN outputs can be used as collective variables, inputs for other collective variables,
or inputs for data analysis tools.

\par Examples

Assume we have an ANN with numbers of nodes being [2, 3, 1], and weights connecting layer 0 and 1 are

[[1,2], [3,4], [5,6]]

weights connecting layer 1 and 2 are

[[7,8,9]]

Bias for layer 1 and 2 are [10, 11, 12] and [13], respectively.

All activation functions are Tanh.

Then if input variables are l_0_out_0, l_0_out_1, the corresponding ANN function object can be defined using
following plumed script:

\plumedfile
ANN ...
LABEL=ann
ARG=l_0_out_0,l_0_out_1
NUM_LAYERS=3
NUM_NODES=2,3,1
ACTIVATIONS=Tanh,Tanh
WEIGHTS0=1,2,3,4,5,6
WEIGHTS1=7,8,9
BIASES0=10,11,12
BIASES1=13
... ANN
\endplumedfile

To access its components, we use "ann.node-0", "ann.node-1", ..., which represents the components of neural network outputs.


*/
//+ENDPLUMEDOC

class ANNB : public Function
{
private:
  int num_layers;
  vector<int> num_nodes;
  vector<string> activations;   // activation functions
  vector<vector<double> > weights;  // flattened weight arrays
  vector<vector<double> > biases;
  vector<vector<double> > output_of_each_layer;
  vector<vector<double> > input_of_each_layer;
  vector<double** > coeff;  // weight matrix arrays, reshaped from "weights"
  vector<vector<double> > expectations_of_batchnorm_layer; // expectaton in batch normalized activation function -- SD
  vector<vector<double> > variances_of_batchnorm_layer; // variance in batch normalized activation function -- SD
  vector<vector<double> > gammas_of_batchnorm_layer; // gamma in batch normalized activation function -- SD
  vector<vector<double> > betas_of_batchnorm_layer; // beta in batch normalized activation function -- SD
  bool apply_batch_norm;

public:
  static void registerKeywords( Keywords& keys );
  explicit ANNB(const ActionOptions&);
  virtual void calculate();
  void calculate_output_of_each_layer(const vector<double>& input);
  void back_prop(vector<vector<double> >& derivatives_of_each_layer, int index_of_output_component);
};

PLUMED_REGISTER_ACTION(ANNB,"ANNB")

void ANNB::registerKeywords( Keywords& keys ) {
  Function::registerKeywords(keys);
  keys.use("ARG"); keys.use("PERIODIC");
  keys.add("compulsory", "NUM_LAYERS", "number of layers of the neural network");
  keys.add("compulsory", "NUM_NODES", "numbers of nodes in each layer of the neural network");
  keys.add("compulsory", "ACTIVATIONS", "activation functions for the neural network");
  keys.add("numbered", "WEIGHTS", "flattened weight arrays connecting adjacent layers, "
           "WEIGHTS0 represents flattened weight array connecting layer 0 and layer 1, "
           "WEIGHTS1 represents flattened weight array connecting layer 1 and layer 2, ...");
  keys.add("numbered", "BIASES", "bias array for each layer of the neural network, "
           "BIASES0 represents bias array for layer 1, BIASES1 represents bias array for layer 2, ...");
  keys.addFlag("APPLY_BATCH_NORM", false, "Set to TRUE if you want to apply batch normalization and then pass "
               " EXPECTATIONS, VARIANCES, GAMMAS and BETA ");
  keys.add("numbered", "EXPECTATIONS", "expectations array of batch norm for each layer of the neural "
           "network, EXPECTATIONS0 represent expectations array for layer 1, EXPECTATIONS1 representation "
           "expectations array for layer 2, ...");
  keys.add("numbered", "VARIANCES", "variances array of batch norm for each layer of the neural "
          "network, VARIANCES0 represent variances array for layer 1, VARIANCES1 representation "
          "variances array for layer 2, ...");
  keys.add("numbered", "BETAS", "betas array of batch norm for each layer of the neural "
          "network, BETAS0 represent betas array for layer 1, BETAS1 representation "
          "betas array for layer 2, ...");
  keys.add("numbered", "GAMMAS", "gammas array of batch norm for each layer of the neural "
          "network, GAMMAS0 represent gammas array for layer 1, GAMMAS1 representation "
          "gammas array for layer 2, ...");
  // since v2.2 plumed requires all components be registered
  keys.addOutputComponent("node", "default", "components of ANNB outputs");
}

ANNB::ANNB(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  parse("NUM_LAYERS", num_layers);
  num_nodes = vector<int>(num_layers);
  activations = vector<string>(num_layers - 1);
  output_of_each_layer = vector<vector<double> >(num_layers);
  input_of_each_layer = vector<vector<double> >(num_layers);
  coeff = vector<double** >(num_layers - 1);
  parseVector("NUM_NODES", num_nodes);
  parseVector("ACTIVATIONS", activations);
  parseFlag("APPLY_BATCH_NORM", apply_batch_norm);
  log.printf("\nactivations = ");
  for (auto ss: activations) {
    log.printf("%s, ", ss.c_str());
  }
  log.printf("\nnum_nodes = ");
  for (auto ss: num_nodes) {
    log.printf("%d, ", ss);
  }
  vector<double> temp_single_coeff, temp_single_bias;
  vector<double> temp_single_batchnorm_expectations, temp_single_batchnorm_variances; // -- for batch norm. SD
  vector<double> temp_single_batchnorm_gammas, temp_single_batchnorm_betas; // -- for batch norm. SD
  int count_batch_norm = 0;
  for (int ii = 0; ; ii ++) {
    // parse coeff
    if( !parseNumberedVector("WEIGHTS", ii, temp_single_coeff) ) {
      temp_single_coeff=weights[ii-1];
      break;
    }
    weights.push_back(temp_single_coeff);
    log.printf("size of temp_single_coeff = %lu\n", temp_single_coeff.size());
    log.printf("size of weights = %lu\n", weights.size());
    // parse bias
    if( !parseNumberedVector("BIASES", ii, temp_single_bias) ) {
      temp_single_bias=biases[ii-1];
    }
    biases.push_back(temp_single_bias);
    log.printf("size of temp_single_bias = %lu\n", temp_single_bias.size());
    log.printf("size of biases = %lu\n", biases.size());
    if (apply_batch_norm and activations[ii] == string("BNTanh") ) {
      // parse expectations of batch norm layer -- SD
      if( !parseNumberedVector("EXPECTATIONS", ii, temp_single_batchnorm_expectations) ) {
        temp_single_batchnorm_expectations=expectations_of_batchnorm_layer[count_batch_norm-1];
      }
      expectations_of_batchnorm_layer.push_back(temp_single_batchnorm_expectations);
      log.printf("size of temp_single_batchnorm_expectations = %lu\n", temp_single_batchnorm_expectations.size());
      log.printf("size of expectations_of_batchnorm_layer = %lu\n", expectations_of_batchnorm_layer.size());
      // parse variances of batch norm layer -- SD
      if( !parseNumberedVector("VARIANCES", ii, temp_single_batchnorm_variances) ) {
        temp_single_batchnorm_variances=variances_of_batchnorm_layer[count_batch_norm-1];
      }
      variances_of_batchnorm_layer.push_back(temp_single_batchnorm_variances);
      log.printf("size of temp_single_batchnorm_variances = %lu\n", temp_single_batchnorm_variances.size());
      log.printf("size of variances_of_batchnorm_layer = %lu\n", variances_of_batchnorm_layer.size());
      // parse gammas of batch norm layer -- SD
      if( !parseNumberedVector("GAMMAS", ii, temp_single_batchnorm_gammas) ) {
        temp_single_batchnorm_gammas=gammas_of_batchnorm_layer[count_batch_norm-1];
      }
      gammas_of_batchnorm_layer.push_back(temp_single_batchnorm_gammas);
      log.printf("size of temp_single_batchnorm_gammas = %lu\n", temp_single_batchnorm_gammas.size());
      log.printf("size of gammas_of_batchnorm_layer = %lu\n", gammas_of_batchnorm_layer.size());
      // parse betas of batch norm layer -- SD
      if( !parseNumberedVector("BETAS", ii, temp_single_batchnorm_betas) ) {
        temp_single_batchnorm_betas=betas_of_batchnorm_layer[count_batch_norm-1];
      }
      betas_of_batchnorm_layer.push_back(temp_single_batchnorm_betas);
      log.printf("size of temp_single_batchnorm_betas = %lu\n", temp_single_batchnorm_betas.size());
      log.printf("size of betas_of_batchnorm_layer = %lu\n", betas_of_batchnorm_layer.size());
      count_batch_norm += 1;
    }
  }

  if(getNumberOfArguments() != num_nodes[0]) {
    error("Number of arguments is wrong");
  }

  auto temp_coeff = weights;
  for (int ii = 0; ii < num_layers - 1; ii ++) {
    int num_of_rows, num_of_cols; // num of rows/cols for the coeff matrix of this connection
    num_of_rows = num_nodes[ii + 1];
    num_of_cols = num_nodes[ii];
    assert (num_of_rows * num_of_cols == temp_coeff[ii].size()); // check whether the size matches
    // create a 2d array to hold coefficients
    coeff[ii] = new double*[num_of_rows];
    for (int kk = 0; kk < num_of_rows; kk ++) {
      coeff[ii][kk] = new double[num_of_cols];
    }
    for (int jj = 0; jj < temp_coeff[ii].size(); jj ++) {
      coeff[ii][jj / num_of_cols][jj % num_of_cols] = temp_coeff[ii][jj];
    }
  }
  // check coeff
  for (int ii = 0; ii < num_layers - 1; ii ++) {
    log.printf("coeff %d = \n", ii);
    for (int jj = 0; jj < num_nodes[ii + 1]; jj ++) {
      for (int kk = 0; kk < num_nodes[ii]; kk ++) {
        log.printf("%f ", coeff[ii][jj][kk]);
      }
      log.printf("\n");
    }
  }
  // check bias
  for (int ii = 0; ii < num_layers - 1; ii ++) {
    log.printf("bias %d = \n", ii);
    for (int jj = 0; jj < num_nodes[ii + 1]; jj ++) {
      log.printf("%f ", biases[ii][jj]);
    }
    log.printf("\n");
  }
  if(apply_batch_norm) {
    // check expectations of batch norm layer -- SD
    int count_batch_norm = 0;
    for (int ii = 0; ii < num_layers - 1; ii ++) {
      if (activations[ii] == string("BNTanh") ){
        log.printf("expectations %d = \n", ii);
        for (int jj = 0; jj < num_nodes[ii + 1]; jj ++) {
          log.printf("%f ", expectations_of_batchnorm_layer[count_batch_norm][jj]);
        }
        count_batch_norm += 1;
        log.printf("\n");
      }
    }
    // check variances of batch norm layer -- SD
    count_batch_norm = 0;
    for (int ii = 0; ii < num_layers - 1; ii ++) {
      if (activations[ii] == string("BNTanh") ){
        log.printf("variances %d = \n", ii);
        for (int jj = 0; jj < num_nodes[ii + 1]; jj ++) {
          log.printf("%f ", variances_of_batchnorm_layer[count_batch_norm][jj]);
        }
      count_batch_norm += 1;
      log.printf("\n");
      }
    } 
    // check gammas of batch norm layer -- SD
    count_batch_norm = 0;
    for (int ii = 0; ii < num_layers - 1; ii ++) {
      if (activations[ii] == string("BNTanh") ){
        log.printf("gammas %d = \n", ii);
        for (int jj = 0; jj < num_nodes[ii + 1]; jj ++) {
          log.printf("%f ", gammas_of_batchnorm_layer[count_batch_norm][jj]);
        }
      count_batch_norm += 1;
      log.printf("\n");
      }
    } 
    // check betas of batch norm layer -- SD
    count_batch_norm = 0;
    for (int ii = 0; ii < num_layers - 1; ii ++) {
      if (activations[ii] == string("BNTanh") ) {
        log.printf("betas %d = \n", ii);
        for (int jj = 0; jj < num_nodes[ii + 1]; jj ++) {
          log.printf("%f ", betas_of_batchnorm_layer[count_batch_norm][jj]);
        }
      count_batch_norm += 1;
      log.printf("\n");
      }
    }
  }

  log.printf("initialization ended\n");
  // create components
  for (int ii = 0; ii < num_nodes[num_layers - 1]; ii ++) {
    string name_of_this_component = "node-" + to_string(ii);
    addComponentWithDerivatives(name_of_this_component);
    componentIsNotPeriodic(name_of_this_component);
  }
  checkRead();
}

void ANNB::calculate_output_of_each_layer(const vector<double>& input) {
  // first layer
  output_of_each_layer[0] = input;
  int count_batch_norm = 0;
  // following layers
  for(int ii = 1; ii < num_nodes.size(); ii ++) {
    output_of_each_layer[ii].resize(num_nodes[ii]);
    input_of_each_layer[ii].resize(num_nodes[ii]);
    // first calculate input
    for (int jj = 0; jj < num_nodes[ii]; jj ++) {
      input_of_each_layer[ii][jj] = biases[ii - 1][jj];  // add bias term
      for (int kk = 0; kk < num_nodes[ii - 1]; kk ++) {
        input_of_each_layer[ii][jj] += coeff[ii - 1][jj][kk] * output_of_each_layer[ii - 1][kk];
      }
    }
    // then get output
    if (activations[ii - 1] == string("Linear")) {
      for(int jj = 0; jj < num_nodes[ii]; jj ++) {
        output_of_each_layer[ii][jj] = input_of_each_layer[ii][jj];
      }
    }
    else if (activations[ii - 1] == string("Tanh")) {
      for(int jj = 0; jj < num_nodes[ii]; jj ++) {
        output_of_each_layer[ii][jj] = tanh(input_of_each_layer[ii][jj]);
      }
    }
    else if (activations[ii - 1] == string("BNTanh")) { // Added batch normalized Tanh -- SD
      for(int jj = 0; jj < num_nodes[ii]; jj ++) {
        output_of_each_layer[ii][jj] = (tanh(input_of_each_layer[ii][jj]) - expectations_of_batchnorm_layer[count_batch_norm][jj]) * \
                                       (gammas_of_batchnorm_layer[count_batch_norm][jj] / variances_of_batchnorm_layer[count_batch_norm][jj]);
        output_of_each_layer[ii][jj] += betas_of_batchnorm_layer[count_batch_norm][jj];
      }
      count_batch_norm += 1;
    }
    else if (activations[ii-1] == string("ReLU")) { // Added ReLU -- SD      
      for(int jj = 0; jj < num_nodes[ii]; ii ++) {
        output_of_each_layer[ii][jj] = 0;
        if (input_of_each_layer[ii][jj] > 0) {
          output_of_each_layer[ii][jj] = input_of_each_layer[ii][jj];
        }
      }
    }
    else if (activations[ii - 1] == string("Circular")) {
      assert (num_nodes[ii] % 2 == 0);
      for(int jj = 0; jj < num_nodes[ii] / 2; jj ++) {
        double radius = sqrt(input_of_each_layer[ii][2 * jj] * input_of_each_layer[ii][2 * jj]
                             +input_of_each_layer[ii][2 * jj + 1] * input_of_each_layer[ii][2 * jj + 1]);
        output_of_each_layer[ii][2 * jj] = input_of_each_layer[ii][2 * jj] / radius;
        output_of_each_layer[ii][2 * jj + 1] = input_of_each_layer[ii][2 * jj + 1] / radius;

      }
    }
    else {
      printf("layer type not found!\n\n");
      return;
    }
  }
#ifdef DEBUG_LAYEROUTPUT
  // print out the result for debugging
  printf("output_of_each_layer = \n");
  // for (int ii = num_layers - 1; ii < num_layers; ii ++) {
  for (int ii = 0; ii < num_layers; ii ++) {
    printf("layer[%d]: ", ii);
    if (ii != 0) {
      cout << activations[ii - 1] << "\t";
    }
    else {
      cout << "input \t" ;
    }
    for (int jj = 0; jj < num_nodes[ii]; jj ++) {
      printf("%lf\t", output_of_each_layer[ii][jj]);
    }
    printf("\n");
  }
  printf("\n");
#endif
  return;
}

void ANNB::back_prop(vector<vector<double> >& derivatives_of_each_layer, int index_of_output_component) {
  derivatives_of_each_layer = output_of_each_layer;  // the data structure and size should be the same, so I simply deep copy it
  // first calculate derivatives for bottleneck layer
  for (int ii = 0; ii < num_nodes[num_nodes.size() - 1]; ii ++ ) {
    if (ii == index_of_output_component) {
      derivatives_of_each_layer[num_nodes.size() - 1][ii] = 1; 
    }
    else {
      derivatives_of_each_layer[num_nodes.size() - 1][ii] = 0;
    }
  }
  // the use back propagation to calculate derivatives for previous layers
  int count_batch_norm = expectations_of_batchnorm_layer.size() - 1;
  float tanh_output_batch_norm;
  for (int jj = num_nodes.size() - 2; jj >= 0; jj --) {
    if (activations[jj] == string("Circular")) {
      vector<double> temp_derivative_of_input_for_this_layer;
      temp_derivative_of_input_for_this_layer.resize(num_nodes[jj + 1]);
#ifdef DEBUG
      assert (num_nodes[jj + 1] % 2 == 0);
#endif
      // first calculate the derivative of input from derivative of output of this circular layer
      for(int ii = 0; ii < num_nodes[jj + 1] / 2; ii ++) {
        // printf("size of input_of_each_layer[%d] = %d\n",jj,  input_of_each_layer[jj].size());
        double x_p = input_of_each_layer[jj + 1][2 * ii];
        double x_q = input_of_each_layer[jj + 1][2 * ii + 1];
        double radius = sqrt(x_p * x_p + x_q * x_q);
        temp_derivative_of_input_for_this_layer[2 * ii] = x_q / (radius * radius * radius)
            * (x_q * derivatives_of_each_layer[jj + 1][2 * ii]
               - x_p * derivatives_of_each_layer[jj + 1][2 * ii + 1]);
        temp_derivative_of_input_for_this_layer[2 * ii + 1] = x_p / (radius * radius * radius)
            * (x_p * derivatives_of_each_layer[jj + 1][2 * ii + 1]
               - x_q * derivatives_of_each_layer[jj + 1][2 * ii]);
      }
#ifdef DEBUG
      for (int mm = 0; mm < num_nodes[jj + 1]; mm ++) {
        printf("temp_derivative_of_input_for_this_layer[%d] = %lf\n", mm, temp_derivative_of_input_for_this_layer[mm]);
        printf("derivatives_of_each_layer[%d + 1][%d] = %lf\n", jj, mm, derivatives_of_each_layer[jj + 1][mm]);
      }
#endif
      // the calculate the derivative of output of layer jj, from derivative of input of layer (jj + 1)
      for (int mm = 0; mm < num_nodes[jj]; mm ++) {
        derivatives_of_each_layer[jj][mm] = 0;
        for (int kk = 0; kk < num_nodes[jj + 1]; kk ++) {
          derivatives_of_each_layer[jj][mm] += coeff[jj][kk][mm] \
                                               * temp_derivative_of_input_for_this_layer[kk];
#ifdef DEBUG
          printf("derivatives_of_each_layer[%d][%d] = %lf\n", jj, mm, derivatives_of_each_layer[jj][mm]);
          printf("coeff[%d][%d][%d] = %lf\n", jj, kk, mm, coeff[jj][kk][mm]);
#endif
        }
      }
      // TODO: should be fine, pass all tests, although there seems to be some problems here previously
    }
    else {
      for (int mm = 0; mm < num_nodes[jj]; mm ++) {
        derivatives_of_each_layer[jj][mm] = 0;
        for (int kk = 0; kk < num_nodes[jj + 1]; kk ++) {
          if (activations[jj] == string("Tanh")) {
            // printf("tanh\n");
            derivatives_of_each_layer[jj][mm] += derivatives_of_each_layer[jj + 1][kk] \
                                                 * coeff[jj][kk][mm] \
                                                 * (1 - output_of_each_layer[jj + 1][kk] * output_of_each_layer[jj + 1][kk]);
          }
          else if (activations[jj] == string("BNTanh")) { // support for batch normalized Tanh -- SD
            // printf("bntanh\n");
            tanh_output_batch_norm = ( (output_of_each_layer[jj + 1][kk] - betas_of_batchnorm_layer[count_batch_norm][kk] ) \
                                     * variances_of_batchnorm_layer[count_batch_norm][kk] \
                                     * (1/gammas_of_batchnorm_layer[count_batch_norm][kk]) ) \
                                     + expectations_of_batchnorm_layer[count_batch_norm][kk];
            
            derivatives_of_each_layer[jj][mm] += derivatives_of_each_layer[jj + 1][kk] \
                                                 * coeff[jj][kk][mm] \
                                                 * (1 - tanh_output_batch_norm * tanh_output_batch_norm) \
                                                 * (gammas_of_batchnorm_layer[count_batch_norm][kk]) \
                                                 * (1/variances_of_batchnorm_layer[count_batch_norm][kk]);
          }
          else if (activations[jj] == string("Linear")) {
            // printf("linear\n");
            derivatives_of_each_layer[jj][mm] += derivatives_of_each_layer[jj + 1][kk] \
                                                 * coeff[jj][kk][mm] \
                                                 * 1;
          }
          else if (activations[jj] == string("ReLU")) {
              // printf("ReLU\n") // Added ReLU -- SD
            if (output_of_each_layer[jj + 1][kk] <= 0) {
              derivatives_of_each_layer[jj][mm] += 0;
            }
            else {
              derivatives_of_each_layer[jj + 1][kk] += derivatives_of_each_layer[jj + 1][kk] \                                
                                                    * coeff[jj][kk][mm] \
                                                    * 1;
            }
          }
          else {
            printf("layer type not found!\n\n");
            return;
          }
        }
      }
    }
    if (activations[jj] == string("BNTanh")){
      count_batch_norm -= 1;
    }
  }
#ifdef DEBUG_LAYERDERIVATIVES
  // print out the result for debugging
  printf("derivatives_of_each_layer = \n");
  for (int ii = 0; ii < num_layers; ii ++) {
    printf("layer[%d]: ", ii);
    for (int jj = 0; jj < num_nodes[ii]; jj ++) {
      printf("%lf\t", derivatives_of_each_layer[ii][jj]);
    }
    printf("\n");
  }
  printf("\n");
#endif
  return;
}

void ANNB::calculate() {

  vector<double> input_layer_data(num_nodes[0]);
  for (int ii = 0; ii < num_nodes[0]; ii ++) {
    input_layer_data[ii] = getArgument(ii);
  }
  
  calculate_output_of_each_layer(input_layer_data);
  vector<vector<double> > derivatives_of_each_layer;

  for (int ii = 0; ii < num_nodes[num_layers - 1]; ii ++) {
    back_prop(derivatives_of_each_layer, ii);
    string name_of_this_component = "node-" + to_string(ii);
    Value* value_new=getPntrToComponent(name_of_this_component);
    value_new -> set(output_of_each_layer[num_layers - 1][ii]);
    for (int jj = 0; jj < num_nodes[0]; jj ++) {
      value_new -> setDerivative(jj, derivatives_of_each_layer[0][jj]); // TODO: setDerivative or addDerivative?
    }
#ifdef DEBUG_FINALDERIVATIVES
    printf("derivatives = ");
    for (int jj = 0; jj < num_nodes[0]; jj ++) {
      printf("%f ", value_new -> getDerivative(jj));
    }
    printf("\n");
#endif
  }

}

}
}
}
