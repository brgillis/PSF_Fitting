""" @file smart_logging.py

    Created 7 Sep 2015

    Logging functions which will use elements if it's available and the
    default Python logger if not.

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

try:
    import ElementsKernel.Logging as log
    
    # Log with Elements
    
except ImportError:
    import logging as log
    
    # Log with Python logger
    
def getLogger(name=None):
    return log.getLogger(name)

def get_logger(name=None):
    return log.getLogger(name)

def set_up_default_logger(name=None):
    logger = getLogger(name)
    logger.setLevel(log.INFO)
    
    handler = log.StreamHandler()
    handler.setLevel(log.INFO)

    formatter = log.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    handler.setFormatter(formatter)
    
    logger.addHandler(handler)
    
    return logger
    
    