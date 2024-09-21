# 50-BOA (IC50-based optimal approach)
This is a MATLAB and R code for precise and efficient estimation of reversible enzyme inhibition constants from initial velocity data.
## System requirements
1. The MATLAB package is based on MATLAB R2023a.
2. The R package is based on R version 4.3.1.
3. The test Excel file is based on version 16.89.
4. The package is tested on macOS M1.
5. The package does not need any non-standard hardware.
## Installation guide
1. The codes can be run directly without installation.
2. No install time is needed.
# Code description
The computational package 50-BOA (IC_50-Based Optimal Approach) includes three MATLAB and R codes: two auxiliary functions, 'BOA_Condition' and 'CV_Inhibition', and a main function, 'Error_Landscape'.

The 'BOA_Condition' function checks whether the given data is sufficient for precise and efficient estimation.

The 'CV_Inhibition function' performs cross-validation on the given data to select an appropriate regularization constant, which is utilized for estimations. 

Utilizing the two auxiliary functions, the main function, 'Error_Landscape', executes the 50-BOA for given data, estimating inhibition constants and generating the error landscape.

## 1. Generate the enzyme inhibition data file formatted by Excel
To use 50-BOA package, users need to organize their initial velocity data obtained with an inhibitor in a defined format. The data should follow the format below:
|Column 1|Column 2|Column 3|Column 4|
|---|---|---|---|
|$$V_{max}$$|$$K_{M}$$|$$IC_{50}$$|$$S_T$$ for estimating $$IC_{50}$$|
|$$S_{T,1}$$|$$I_{T,1}$$|$$V_{0,1}$$|0|
|...|...|...|...|
|$$S_{T,n}$$|$$I_{T,n}$$|$$V_{0,n}$$|0|

The first row should contain values for the maximal formation rate ($$V_{max}$$) and Michalis-Menten constant ($$K_M$$) of substrate, half maximal inhibitory concentration ($$IC_{50}$$), and substrate concentration ($$S_T$$) used for estimating $$IC_{50}$$, in that order.

The subsequent rows should list the experimental setups and initial velocities obatined with following those setups. The first and second columns contain substrate ($$S_{T,1}$$, ..., $$S_{T,n}$$) and inhibitor ($$I_{T,1}$$, ..., $$I_{T,n}$$) concentrations, respectively. The third column lists the resulting initial velocitiy data ($$V_{0,1}$$, ..., $$V_{0,n}$$). The fourth column is a placeholder to match the size of the first row.
## 2. Run the 'Error_Landscape.m' or 'Error_Landscape.R' file.
Once the data is formatted correctly, users can run the package by executing the following line:

    Error_Landscape(dir)
Here, 'dir' is the file directory of the Excel-formatted data. Then function 'BOA_Condition' in package automatically checks the conditions for precise and efficient estimation.

If there is no $$I_T$$ value greater than or equal to $$IC_{50}$$, or if there are not two distinct $$S_T$$ ranged within $$(0.2K_M, 5K_M)$$, 'BOA_Condition' will indicate that the experimental data is insufficient and explain why.

Otherwise, function 'CV_Inhibition' will perform regularization constant selection, and Error_Landscape will estimate inhibition constants.

The output includes the estimated inhibition constants with their 95% confidence intervals, the regularization constant used, and a heatmap of the error landscape.


# Test examples
We provide two test example files, named as 'test_1.xlsx' and 'test_2.xlsx'.

test_1 allows for precise estimation, while test_2 does not.

For test_1, the output is the following:

    >> Error_Landscape("test_1.xlsx")
    The regularization constant is 756.46.
    Kic: 0.0398, (0.0358, 0.0460)
    Kiu: 0.0403, (0.0337, 0.0482)

For test_2, the output is the following:

    >> Error_Landscape("test_2.xlsx")
    Estimation may be insufficient for precise results:
    inhibitor concentration < IC50
    Substrate concentration should vary
    The regularization constant is 0.01.
    Kic: 0.0367, (0.0244, 1592478110290.1743)
    Kiu: 0.0766, (0.0208, 59985261653.5099)
