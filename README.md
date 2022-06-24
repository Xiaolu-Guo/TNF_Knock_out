# TNF_Knock_out
Simulation for non-oscillation NFkB signal in TNF-/- in repsonse to TNF.
(The codes are tested in Matlab R2020a)

[Step 1]: Download the whole folder to local.

[Step 2]: In matlab, CHANGE the current folder to 'TNF_Knock_our'(or the folder where 'run_me.m' is in). 

Step 2.1 (Optional) : One can set different values of the variables run_, draw_ to control which task should be done. For example:
    run_sample_TNFo = 0;

Step 2.2 (Optional) : One can set different simulation size (10 cells is simulated by default) through changing the value of 'Num_sample'. For example:
    Num_sample = 100;

Step 2.3 (Optional) : change the data/figure saving location by changing the value of 'data_save_file_path' and 'fig_save_path'. For example:
    data_save_file_path = '../raw_data/';

[Step 3]: Run 'run_me.m' to get all the simulation and results.

Output:
By default, it will create ./data ./figure under the current folder in Matlab for saving data and figures.

( ignore update.sh. this is for updating the codes in lib which comes from other projects)
