""" @file report_results.py

    Created 5 Nov 2015

    Reports the results of a PSF test, printing basic results in the
    console and detailed results in a fits table.

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

from psf_testing import magic_values as mv

def report_results(test_results, filename_root):

    # Set up the columns first
    columns = []

    for label, i, err_i in zip(("star_", "model_", "noisy_model_"), (4, 5, 6), (7, None, 7)):

        if err_i is None:
            get_errs = np.zeros_like
            err_i = i
        else:
            get_errs = lambda x : x

        columns.append(fits.Column(name=label + "m0", format='E', array=test_results[i][0]))
        columns.append(fits.Column(name=label + "m0" + "_err", format='E', array=get_errs(test_results[err_i][0])))

        for Q_label, j in zip(("Qx", "Qy", "Qplus", "Qcross", "Qsize"), range(5)):

            columns.append(fits.Column(name=label + Q_label, format='E', array=test_results[i][1][j]))
            columns.append(fits.Column(name=label + Q_label + "_err", format='E',
                                       array=get_errs(test_results[err_i][1][j])))

    columns.append(fits.Column(name="is_outlier", format='L', array=test_results[8][1]))

    tbhdu = fits.BinTableHDU.from_columns(columns)

    # Add metadata to the header

    tbhdu.header["FOCUS"] = test_results[0]
    tbhdu.header["CHI_SQR"] = test_results[1][0]
    tbhdu.header["DOF"] = test_results[1][1]

    tbhdu.header["M0_DIFF"] = test_results[2][0][0]
    tbhdu.header["M0_Z"] = test_results[3][0][0]

    for Q_label, j in zip(("QX", "QY", "QP", "QC", "QS"), range(5)):
        tbhdu.header[Q_label + "_DIFF"] = test_results[2][0][1][j]
        tbhdu.header[Q_label + "_Z"] = test_results[3][0][1][j]

    tbhdu.header["M0_NDIFF"] = test_results[2][1][0]
    tbhdu.header["M0_NZ"] = test_results[3][1][0]

    for Q_label, j in zip(("QX", "QY", "QP", "QC", "QS"), range(5)):
        tbhdu.header[Q_label + "_NDIFF"] = test_results[2][1][1][j]
        tbhdu.header[Q_label + "_NZ"] = test_results[3][1][1][j]


    results_filename = filename_root + "_results" + mv.table_extension

    tbhdu.writeto(results_filename)

    return
