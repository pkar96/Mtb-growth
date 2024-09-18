Executable file
analyse_fn.m - generate elongation speed vs age and growth rate vs age plots for acidic or neutral medium.

Excel file
pH_7 - Contains HADA label data for multiple cells growing in neutral medium 
pH_59 - Contains HADA label data for multiple cells growing in acidic medium
The cells 
Column 1 is 0 (new pole), Column 4 - Length of the cell

Additional files
binning.m, binning_with_error.m and, binning_with_error_1.m - Used to do binning in the plots. binning_with_error.m and binning_with_error_1.m also return error estimates for the binned data.

birth_to_div_gr_v_age_sims.m- Simulates linear growth and returns elongation speed, growth rate at particular ages to be plotted in call.m function.
birth_to_div_gr_v_age_exp_sims.m- Simulates exponential growth and returns elongation speed, growth rate at particular ages to be plotted in call.m function.

params.xlsx - Contains parameters which are passed to birth_to_div_gr_v_age_sims.m
