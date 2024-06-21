The analysis and simulations were carried out using MATLAB R2021b but might work with previous and latest versions of MATLAB. 

To run these files, download and unzip the files in any folder of your choice. The instructions to run the executable files are in the header comments of those files. Each of these files takes less than a minute to execute, and in most cases, only a few seconds.

Folders with executable codes in them.
1. simulate OETO NETO BEITO - contains executable files to simulate exponentially growing cells as well as cells following a linear/bilinear (BEITO, NETO, OETO) growth pattern.
2. HADA movies - To generate elongation speed vs age and growth rate vs age plots for acidic or neutral medium.
3. NET0_v_BETO - To generate details about polar growth such as when a cell pole starts growing and its elongation speed and growth amount.

We have an executable file here. It is an example of how to generate growth rate vs age and elongation speed vs age plots. The data it uses (stored in data folder) is cell length vs time for cells growing in an unbuffered medium (Fig. 4A and 4D in Chung et al., 2023).
call.m - Calls analysis_06_24.m, analysis_09_03.m, and analysis_09_03_2.m to generate elongation speed vs age and growth rate vs age plots for unbuffered media experimental data.

Additional files
binning.m, binning_with_error.m and, binning_with_error_1.m - Used to do binning in the plots. binning_with_error.m and binning_with_error_1.m also return error estimates for the binned data.
birth_to_div_gr_v_age_sims.m- Simulates linear growth and returns elongation speed, growth rate at particular ages to be plotted in call.m function.

Excel file
params.xlsx - Contains parameters which are passed to birth_to_div_gr_v_age_sims.m
