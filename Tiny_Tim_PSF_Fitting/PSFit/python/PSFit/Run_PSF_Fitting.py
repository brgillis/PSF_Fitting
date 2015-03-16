#!/usr/bin/env python

import sys
import subprocess as sbp

# Magic values
template_job_name = "test_job"
template_job_file = "/disk4/brg/Program_Files/PSF_Fitting/test_job.pbs"
template_field_name = "/disk4/brg/Data/PSF_Fitting/jb1606cuq_sci2_cor"
template_field_name_rx = template_field_name.replace("\\","\\\\").replace("/","\/").replace(".","\.")
template_job_parameters = "0.99 22 8 4 True True"
template_job_parameters_rx = template_job_parameters.replace("\\","\\\\").replace("/","\/").replace(".","\.")

def main(argv):
    # Check that the name of a file listing fields was passed at the command line
    if(len(argv)) <= 1:
        raise Exception("Name of fields list must be passed at command-line.\n" + \
                        "eg. python Run_PSF_Fitting.py fields_list.txt")
        
    with open(argv[1], 'r') as fi:
        # Read in the file, except for comment lines
        lines = []
        for line in fi:
            line.strip()
            if((line[0][0] != '#') and (len(line) > 0)):
                lines.append(line.split())
                
    if(len(lines) <=0):
        print("No fields found to process.")
        return
    
    if(len(argv) >= 3):
        # The next value in the arguments list should be the name for the summary file
        summary_file_name = argv[2]
        
        # Set up header for the summary file
        with open(summary_file_name, 'w') as fo:
            fo.write("# field_name\tBest_focus\tBest_Chi2\tBest_Chi2_core\tBest_Chi2_wings\t" +
                     "Best_Chi2_size_shape_core\tBest_Chi2_size_shape_wings\t" +
                     "Star_Chi2_mean\tStar_Chi2_stddev\tStar_Chi2_stderr\t" +
                     "Star_Chi2_outlier_frac\tChip\t" +
                     "Core_dp_x_Chi2\tCore_dp_y_Chi2\tCore_qp_xx_Chi2\tCore_qp_xy_Chi2\t" +
                     "Wings_dp_x_Chi2\tWings_dp_y_Chi2\tWings_qp_xx_Chi2\tWings_qp_xy_Chi2\t" +
                     "Core_size_diff_Chi2\tCore_e1_diff_Chi2\tCore_e2_diff_Chi2\t" +
                     "Wings_size_diff_Chi2\tWings_e1_diff_Chi2\tWings_e2_diff_Chi2\t" +
                     "Core_dp_x\tCore_dp_y\tCore_qp_xx\tCore_qp_xy\t" +
                     "Wings_dp_x\tWings_dp_y\tWings_qp_xx\tWings_qp_xy\t" +
                     "Core_size_diff\tCore_e1_diff\tCore_e2_diff\t" +
                     "Wings_size_diff\tWings_e1_diff\tWings_e2_diff\t" +
                     "\n")
        
    else:
        summary_file_name = ""
                
    for line in lines:
        field_name = line[0]
        job_name = field_name.replace(".fits","")
        pbs_file = job_name + ".pbs"
        if(summary_file_name==""):
            job_parameters = template_job_parameters
        else:
            job_parameters = template_job_parameters + " " + summary_file_name
        
        # Set up the run script, using the template as a base
        
        # We'll do three substitute commands with awk, the first to change the name of the job:
        job_name_change_script = "awk '{sub(/" + template_job_name + "/,\"" + job_name + \
            "\")}; 1' "
        
        # The second to change the field it acts on:
        field_change_script = "awk '{sub(/" + template_field_name_rx + "/,\"" + field_name + \
            "\")}; 1' "
        
        # And the third to change the parameters to the job
        job_change_script = "awk '{sub(/" + template_job_parameters_rx + "/,\"" + job_parameters + \
            "\")}; 1' "
        
        # And put the commands together and call:
        cmd = job_name_change_script + " " + template_job_file + " | " + \
              field_change_script + " | " + job_change_script + " > " + pbs_file
        sbp.call(cmd,shell=True)
        
        # Now submit the script
        cmd = "qsub " + pbs_file
        sbp.call(cmd,shell=True) # Comment out if just testing generation
    

if __name__ == "__main__":
    main(sys.argv)
