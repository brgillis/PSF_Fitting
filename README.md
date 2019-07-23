# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

In the file call_tinytim.py, you'll see a line at the top that says: `tinytim_path = "/home/brg/Program_Files/tinytim-7.5"` You'll have to change this to the local path to your installation of Tiny Tim.

Once you've done that, it should be able to be run easily. You can invoke `python PSF_Testing/run_psf_testing.py --help` to see the command-line options. Some example commands are:

Test a single focus:

python run_psf_testing.py --image_filename /disk2/brg/Data/HST_Fields/jb5d07hoq_sci1_cor.fits --min_lowest_separation 1.0 --min_class_star 0.95 --min_star_mag 22 --max_star_mag 25 --focus 2.1875 --norm_errors --focus_sample_x_points 32 --focus_sample_y_points 16 --min_star_snr 50 --subsampling_factor 8 --tinytim_data_path /home/brg/Data/HST_Fields/PSF_models

Fit the focus:

python run_psf_testing.py --image_filename /disk2/brg/Data/HST_Fields/jb5d07hoq_sci1_cor.fits --min_lowest_separation 1.0 --min_class_star 0.95 --min_star_mag 22 --max_star_mag 25 --focus_samples 7 --norm_errors --focus_sample_x_points 32 --focus_sample_y_points 16 --min_star_snr 50 --subsampling_factor 8 --tinytim_data_path /home/brg/Data/HST_Fields/PSF_models

Fit the focus using the updated set of optical parameters I determined in my paper (this should provide a slightly better quality of fit in most cases):

python run_psf_testing.py --image_filename /disk2/brg/Data/HST_Fields/jb5d07hoq_sci1_cor.fits --min_lowest_separation 1.0 --min_class_star 0.95 --min_star_mag 22 --max_star_mag 25 --focus_samples 7 --norm_errors --focus_sample_x_points 32 --focus_sample_y_points 16 --min_star_snr 50 --subsampling_factor 8 --tinytim_data_path /home/brg/Data/HST_Fields/PSF_models --z2 -0.0042 --z3 0.0046 --astigmatism_0 0.0241 --astigmatism_45 0.0300 --coma_x 0.0159 --coma_y -0.0000 --clover_x 0.0074 --clover_y 0.0163 --spherical_3rd -0.0217 --z12 0.0037 --z13 0.0001 --z14 0.0043 --z15 0.0059 --z16 -0.0061 --z17 0.0059 --z18 0.0039 --z19 0.0020 --z20 -0.0008 --z21 0.0072 --spherical_5th 0.0101 --kernel_adjustment 0.9978

At present, the script only works with the WFC detector. For others, you'll have to modify the if/elsif ladder within the get_image_info() function in main.py. Run tiny1 to find out the number for the detector you're using, and enter the name of it in the ladder using the name in the DETECTOR entry of the header file, in all lowercase.

### Who do I talk to? ###

* Talk to Bryan Gillis (brgillis)