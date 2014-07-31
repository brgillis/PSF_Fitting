# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

In the file call_tinytim.py, you'll see a line at the top that says: `tinytim_path = "/home/brg/Program_Files/tinytim-7.5"` You'll have to change this to the local path to your installation of Tiny Tim.

Once you've done that, it should be able to be run easily. You can invoke main.py at the command-line with the name of the fits file you wish to fit the PSF for, with optional other commands as detailed within the main.py file. For a quick test, see if it works with the following command:

`python main.py ???.fits 0.99 14 1 1`

This will generate a single PSF for the image (using the centre of the image as the position and the default focus of 0) and test subtracting it from the brightest stars (those with magnitude < 14), then report the quality of the fit. If that works, you can then ask it to attempt to fit
the best focus, which will take a bit longer. For a full fit, I'd recommend you use the command:

python main.py ???.fits 0.99 22 8 4 True

This will take a bit of time to run, but it should get a much better fit in the end.

At present, the script only works with the WFC detector. For others, you'll have to modify the if/elsif ladder within the get_image_info() function in main.py. Run tiny1 to find out the number for the detector you're using, and enter the name of it in the ladder using the name in the DETECTOR entry of the header file, in all lowercase.

### Who do I talk to? ###

* Talk to Bryan Gillis (brgillis)