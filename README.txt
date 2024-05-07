Folders with executable codes in them

1. simulate OETO NETO BEITO - contains executable files to simulate exponentially growing cells as well as cells following a linear/bilinear (BEITO, NETO, OETO) growth pattern.
2. HADA movies - To generate elongation speed vs age and growth rate vs age plots for acidic or neutral medium.
3. NET_v_BETO - To generate details about polar growth such as when does a cell pole start growing and its elongation speed and growth amount.

Folders:
data - Data is stored in this folder. processed data contains length vs time data for cells growing in unbuffered medium.

Files executable 

call.m - Calls analysis_06_24.m, analysis_09_03.m, and analysis_09_03_2.m to generate elongation speed vs age and growth rate vs age plots for unbuffered media experimental data.

Additional files
binning.m, binning_with_error.m and, binning_with_error_1.m - Used to do binning in the plots. binning_with_error.m and binning_with_error_1.m also return error estimates for the binned data.

birth_to_div_gr_v_age_sims.m- Simulates linear growth and returns elongation speed, growth rate at particular ages to be plotted in call.m function.

Excel file
params.xlsx - Contains parameters which are passed to birth_to_div_gr_v_age_sims.m