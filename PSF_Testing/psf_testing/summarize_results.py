""" @file summarize_results.py

    Created 11 Nov 2015

    @TODO: File docstring

    ---------------------------------------------------------------------

    Copyright (C) 2015 Bryan R. Gillis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from astropy.io import fits
import numpy as np

from utility import smart_logging
from psf_testing import magic_values as mv

param_colnames = (("Z2", "Z 2"),
                  ("Z3", "Z 3"),
                  ("ASTIG_0", "0 degree astigmatism"),
                  ("ASTIG_45", "45 degree astigmatism"),
                  ("COMA_X", "X coma"),
                  ("COMA_Y", "Y coma"),
                  ("CLOVER_X", "X clover"),
                  ("CLOVER_Y", "Y clover"),
                  ("SPHERE_3", "Spherical 3rd"),
                  ("Z12", "Z 12"),
                  ("Z13", "Z 13"),
                  ("Z14", "Z 14"),
                  ("Z15", "Z 15"),
                  ("Z16", "Z 16"),
                  ("Z17", "Z 17"),
                  ("Z18", "Z 18"),
                  ("Z19", "Z 19"),
                  ("Z20", "Z 20"),
                  ("Z21", "Z 21"),
                  ("SPHERE_5", "Spherical 5th"),
                  ("KERN_ADJ", "Kernel adjustment")
                  )

def make_results_summary(results_filename_roots,
                         summary_filename):
    
    image_filenames = []
    chips = []
    obs_times = []
    exp_times = []
    
    ras = []
    decs = []
    
    focii = []
    
    X_squareds = []
    chi_squareds = []
    
    num_stars = []
    dofs = []
    
    m0_diff_diffs = []
    m0_Zs = []

    m0_noisy_diff_diffs = []
    m0_noisy_Zs = []
    
    Qs = {}
    param_vals = {}
    
    logger = smart_logging.get_default_logger()

    for Q_label in ("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"):
        Qs[Q_label + "_DIF"] = []
        Qs[Q_label + "_Z2"] = []
        Qs[Q_label + "NDIF"] = []
        Qs[Q_label + "NZ2"] = []

    for Q_label in ("QXC", "QYC", "QPC", "QCC", "QSC", "QXW", "QYW", "QPW", "QCW", "QSW"):
        Qs[Q_label + "_DIF"] = []
        
    for param in param_colnames:
        param_vals[param[0]] = []
        
    for weight_label in ("core_", "wings_"):
        for label in ("star_", "model_"):
            for Q_label in ("Qx", "Qy", "Qplus", "Qcross", "Qsize"):
                Qs[weight_label+label+Q_label + "_mean"] = []
                
    
    for results_filename_root in results_filename_roots:
        results_filename = results_filename_root + mv.results_tail
        try:
            results_file = fits.open(results_filename)
            _ = results_file[1].header["QXC_DIF"]
        except IOError as _e:
            logger.warn("File " + results_filename + " cannot be opened and will be skipped.")
            continue
        except KeyError as _e:
            logger.warn("File " + results_filename + " is out of date and will be skipped.")
            continue
            
        header = results_file[1].header
        
        image_filenames.append(results_filename.split('/')[-1].replace(mv.results_tail,
                                                                       mv.image_extension))
        
        chips.append(header['CCDCHIP'])
        obs_times.append(header['OBS_TIME'])
        exp_times.append(header['EXP_TIME'])
        
        ras.append(header['RA_TARG'])
        decs.append(header['DEC_TARG'])

        focii.append(header["FOCUS"])
        
        X_squareds.append(header["X_SQR"])
        chi_squareds.append(header["CHI_SQR"])
        
        try:
            num_stars.append(header["NSTAR"])
            dofs.append(header["DOF"])
        except KeyError as _e:
            dof = header["CDOF"]+1
            num_stars.append((header["XDOF"]+1)/dof)
            dofs.append(dof)
            
        
        m0_diff_diffs.append(header["M0D_DIF"])
        m0_Zs.append(header["M0D_Z2"])
    
        m0_noisy_diff_diffs.append(header["M0DNDIF"])
        m0_noisy_Zs.append(header["M0DNZ2"])
    
        for Q_label in ("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"):
            Qs[Q_label + "_DIF"].append(header[Q_label + "_DIF"])
            Qs[Q_label + "_Z2"].append(header[Q_label + "_Z2"])
            Qs[Q_label + "NDIF"].append(header[Q_label + "NDIF"])
            Qs[Q_label + "NZ2"].append(header[Q_label + "NZ2"])
            
        for Q_label in ("QXC", "QYC", "QPC", "QCC", "QSC", "QXW", "QYW", "QPW", "QCW", "QSW"):
            Qs[Q_label + "_DIF"].append(header[Q_label + "_DIF"])
            
        for param in param_colnames:
            param_vals[param[0]].append(header[param[0]])
            
        data = results_file[1].data
        good_stars = data["is_not_outlier"]
        mask = np.logical_not(good_stars)

        for weight_label in ("core_", "wings_"):
            for label in ("star_", "model_"):
                for Q_label in ("Qx", "Qy", "Qplus", "Qcross", "Qsize"):
                    good_values = np.ma.masked_array(data[weight_label+label+Q_label],mask)
                    mean = np.mean(good_values)
                    Qs[weight_label+label+Q_label + "_mean"].append(mean)
                    
                    
    
    columns = [fits.Column(name="filename", format='30A', array=image_filenames),
               fits.Column(name="chip", format='B', array=chips),
               fits.Column(name="obs_time", format='E', array=obs_times),
               fits.Column(name="exp_time", format='E', array=exp_times),
               fits.Column(name="ra", format='E', array=ras),
               fits.Column(name="dec", format='E', array=decs),
               fits.Column(name="focus", format='E', array=focii),
               fits.Column(name="X_squared", format='E', array=X_squareds),
               fits.Column(name="chi_squared", format='E', array=chi_squareds),
               fits.Column(name="num_stars", format='E', array=num_stars),
               fits.Column(name="dofs", format='E', array=dofs),
               fits.Column(name="m0_diff_diff", format='E', array=m0_diff_diffs),
               fits.Column(name="m0_Z2", format='E', array=m0_Zs)]
    
    for Q_label, colname in zip(("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"), 
                                ("Qx_diff", "Qy_diff",
                                 "Qplus_sum", "Qcross_sum", "Qsize_sum",
                                 "Qplus_diff", "Qcross_diff", "Qsize_diff")):
        columns.append(fits.Column(name=colname + "_diff", format='E', array=Qs[Q_label + "_DIF"]))
        columns.append(fits.Column(name=colname + "_Z2", format='E', array=Qs[Q_label + "_Z2"]))
        
    for Q_label, colname in zip(("QXC", "QYC", "QPC", "QCC", "QSC", "QXW", "QYW", "QPW", "QCW", "QSW",), 
                                ("Qx_core", "Qy_core",
                                 "Qplus_core", "Qcross_core", "Qsize_core",
                                 "Qx_wings", "Qy_wings",
                                 "Qplus_wings", "Qcross_wings", "Qsize_wings",)):
        columns.append(fits.Column(name=colname + "_diff_mean", format='E', array=Qs[Q_label + "_DIF"]))
        
    columns += [fits.Column(name="m0_noisy_diff_diff", format='E', array=m0_noisy_diff_diffs),
                fits.Column(name="m0_noisy_Z2", format='E', array=m0_noisy_Zs)]
    
    for Q_label, colname in zip(("QXD", "QYD", "QPS", "QCS", "QSS", "QPD", "QCD", "QSD"), 
                                ("Qx_diff", "Qy_diff",
                                 "Qplus_sum", "Qcross_sum", "Qsize_sum",
                                 "Qplus_diff", "Qcross_diff", "Qsize_diff")):
        columns.append(fits.Column(name=colname + "_noisy_diff", format='E', array=Qs[Q_label + "NDIF"]))
        columns.append(fits.Column(name=colname + "_noisy_Z2", format='E', array=Qs[Q_label + "NZ2"]))
        
    
    for weight_label in ("core_", "wings_"):
        for label in ("star_", "model_"):
            for Q_label in ("Qx", "Qy", "Qplus", "Qcross", "Qsize"):
                name = weight_label+label+Q_label + "_mean"
                columns.append(fits.Column(name=name, format='E', array=Qs[name]))
                
    for param in param_colnames:
        columns.append(fits.Column(name=param[1],format='E',array=param_vals[param[0]]))

    tbhdu = fits.BinTableHDU.from_columns(columns)

    tbhdu.writeto(summary_filename,clobber=True)
    
    logger.info("Results summary output to " + summary_filename + ".")
    
    return
        

def make_stack_stacks(results_filename_roots,
                      stack_stack_filename_root):
    
    logger = smart_logging.get_default_logger()
    
    stack_types = ["star", "model", "noisy_model", "residual"]
    
    stacks = {}
    
    for stack_type in stack_types:
            
        tot_weight = 0
        
        # Combine all stacks of this type together
        for results_filename_root in results_filename_roots:
            
            # Get the filename
            filename = results_filename_root + "_" + stack_type + "_stack" + mv.image_extension
            try:
                hdu = fits.open(filename)[0]
                image = hdu.data
            except IOError as _e:
                logger.warn("File " + filename + " cannot be opened and will be skipped.")
                continue
            
            if "TOTAL_W" in hdu.header:
                weight = hdu.header["TOTAL_W"]
            else:
                weight = 1
            
            if stack_type in stacks:
                # Presently assuming all will be the same size
                stacks[stack_type] += image*weight
            else:
                stacks[stack_type] = image*weight
                
            tot_weight += weight
                
        # Normalize
        stacks[stack_type] /= tot_weight
        
        # Save it
        output_filename = stack_stack_filename_root + "_" + stack_type + "_full_stack" + mv.image_extension
        
        fits.writeto(output_filename, stacks[stack_type], clobber=True)
        
        logger.info("Stack of " + stack_type + " stacks saved to " + output_filename + ".")
        
    return
        
