# Nonparametric Inference on Network Effects

## Introduction

The code is for reproducing the numerical experiments and data example results in 'Optimal Nonparametric Inference on Network Effects with Dependent Edges'.

## Guidelines for Result Replication

### Results in the Main File
  
1) **Table 3:** Please go to folder `./simulation_main/empirical_typeIerror/` and run `typeIerror_and_qqplot.R` in R to obtain the results. The implementation code for experiments is located in `./configurations/`. To save time, one can use the intermediate results saved in the folder `./outputs/` to calculate the empirical sizes.

2) **Table 4:** Please go to folder `./simulation_main/empirical_power/` and run `power.R` in R to obtain the results. The implementation code for experiments is located in `./configurations/`. To save time, one can use the intermediate results saved in the folder `./outputs/` to calculate the empirical powers.

### Results in Appendix

1) **Figure S1:** Please go to folder `./simulation_main/empirical_typeIerror/` and run `typeIerror_and_qqplot.R` in R to obtain the results. The implementation code for experiments is located in `./configurations/`. To save time, one can use the intermediate results saved in the folder `./outputs/` to calculate the empirical sizes.
   
2) **Table S2 and Figure S2:** Please run `./simulation_main/empirical_typeIerror/typeIerror_and_qqplot.R` to obtain the results.

3) **Figure S3:** Please go to folder `./simulation_appendix/` and run `lambda_plot.R` in R to obtain the results. The implementation code for experiments is located in `./configurations/`. To save time, one can use the intermediate results saved in the folder `./outputs/` to calculate the empirical sizes.

2) **Table S4 and Figure S5:** Please go to folder `./simulation_appendix/C-sensitivity/outputs and plots` and run `plot-mis-diag-rate.R` and `plot-size-and-power.R` in R to obtain the results. The implementation code for experiments is located in `./simulation_appendix/C-sensitivity/configurations`. To save time, one can use the intermediate results saved in the folder `./simulation_appendix/C-sensitivity/outputs and plots` to calculate the empirical sizes.

1) **Table S4:** Please go to folder `./application_faculty_hiring_networks/` and run `*.R` to obtain the results, where `*` is `Business`, `CS`, and `History`. The data matrices are in folder `./ready_to_use_data_matrices/`, and the original data files are given in paper <a href="https://www.science.org/doi/10.1126/sciadv.1400005">Systematic inequality and hierarchy in faculty hiring networks</a>.

4) **Table S5 and Figure S8:** First download the data in paper <a href="https://scholar.harvard.edu/melitz/publications/estimating-trade-flows-trading-partners-and-trading-volumes">Estimating trade flows: Trading partners and trading volumes</a>. Then please run `./application_international_trade_network/trade_network.R` to obtain the results.


## Packed Algorithm Code for Practitioners

`networkEffectsTest.R:` Wrapped functions and examples for proposed network effects tests in paper 'Optimal nonparametric inference on network effects with dependent edges'.


## License

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
