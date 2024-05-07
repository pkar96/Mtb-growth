Executable files
analyse.m - Generates histograms of growth asymmetry and elongation speed of each pole in acidic and neutral medium using data from HADA_anl_comb.xlsx

analyse_3_pt.m - Used to generate excel file HADA_anl_3pt_pH_#.xlsx. It takes HADA label data and analyzes only those cells which go from no-HADA label region at one pole to both. It finds the elongation speed and growth amount of each pole, and the time at which they start growing.

analyse_4_pt.m - Used to generate excel file HADA_anl_4pt_pH_#.xlsx. It takes HADA label data and analyzes only those cells which have no-HADA label region at both poles from start. It finds the elongation speed and growth amount of each pole, and the time at which they start growing.

Excel files
1. HADA_anl_3pt_pH_#.xlsx - Stores the elongation speed and growth amount of each pole, and the time at which they start growing for either acidic or neutral media. It takes HADA label data and analyzes only those cells which go from no-HADA label region at one pole to both.
2. HADA_anl_4pt_pH_#.xlsx - Stores the elongation speed and growth amount of each pole, and the time at which they start growing for either acidic or neutral media. It takes HADA label data and analyzes only those cells which have no-HADA label region at both poles from start.
3. HADA_anl_comb.xlsx - Stores the elongation speed and growth amount of each pole, and the time at which they start growing for both acidic and neutral media. Combination of HADA_anl_3pt_pH_#.xlsx and HADA_anl_4pt_pH_#.xlsx with each sheet storing a particular pH information.

Other .m files
1. find_es_time_auto.m - finds the time of start of growth, elongation speed and growth amount of a particular pole of a cell.
2. fit_bilinear.m - does bilinear fit and returns predicted y values, residuals and obtained fit parameters.
3. fit_poly.m - does linear fit and returns predicted y values, residuals and obtained fit parameters.
4. piecewise.m - contains a description of piecewise function.
5. predict_v.m - returns the amount of growth of a particular pole in a cell as a function of time. Takes parameters such as time at which growth starts and elongation speed.
