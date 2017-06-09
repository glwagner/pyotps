"""This module interfaces with the shell-operated OSU Tidal Prediction Software (OTPS).
"""
from __future__ import division

import os, sys, shutil, subprocess
import numpy as np


class tidaldriver(object):
    def __init__(self, 
            otpspath  = '../../OTPS2', 
            otpstype= 'v1'
            ):
        """Initialize the tidaldriver object and test the 
        interface with OTPS.

            Args:
                otpspath (str): Path to the directory with a compiled
                    version of OTPS which is correctly linked to 
                    downloaded data.

                otpstype (str): Name of the model from which data is 
                    to be extracted. The two options available now are
                    'netcdf' (for the ncdf version of TPXO8-atlas-compact)
                    and 'v1' for the 'v1' version of TPXO8-atlas-compact 
                    available from http://volkov.oce.orst.edu/tides/tpxo8_atlas.html
                    as of June 2, 2017.

            Returns:
                The tidaldriver object.
        """

        self.otpspath = otpspath
        self.otpstype = otpstype

        if otpstype is 'netcdf': 
            self.controlfilename = 'Model_load'
            self.modellist = ['load_file.nc']
        elif otpstype is 'v1':
            self.controlfilename = 'Model_atlas_v1'
            self.modellist = [
                'hf.tpxo8_atlas_30_v1', 
                'uv.tpxo8_atlas_30_v1',
                'grid_tpxo8atlas_30_v1']
        else:
            raise ValueError("otpstype must be 'netcdf' or 'v1'.")

        # Tests and file management
        self.datapath = "{}/DATA".format(self.otpspath)
        self.controlfilepath = "{}/{}".format(self.datapath, self.controlfilename)
        self.origcontrolfilepath = "{}_orig".format(self.controlfilepath)

        if not os.path.isdir(self.datapath):
            raise RuntimeError("The DATA directory {} \n".format(self.datapath) +
                "does not exist. Check that otpspath is correctly specified \n"
                "and that the OTPS data has been downloaded.")
        elif not os.path.isfile(self.controlfilepath):
            raise RuntimeError("The control file {} \n".format(self.controlfilepath) +
                "does not exist.")
                
        # Copy control file and write new one
        shutil.copy2(self.controlfilepath, format(self.origcontrolfilepath))    

        with open(self.controlfilepath, 'w') as controlfile:
            for model in self.modellist:
                controlfile.write("{}/{}\n".format(self.datapath, model))

        # TODO:
        #   1. Check that the result patches Joern's OBCs for the North Atlantic.


    def testdrive(self, withmsg=True):
        """Test drive the tidal driver."""

        # Test lats and lons
        lats = [    9.795,  -29.676,    2.021 ]
        lons = [  103.471,   59.550,   72.294 ]

        (constits, latsout, lonsout, 
            amps, phases) = self.extract_amp_phase(lats, lons)

        if withmsg:
            (flatlats, flatlons) = (latsout.flatten(), lonsout.flatten())
            (flatamps, flatphases) = ({}, {})
            
            labelline = "{:>7s}, {:>7s}:".format('lat', 'lon')
            for c in constits:
                labelline = "{}  {:>7s}".format(labelline, c)
                flatamps[c] = amps[c].flatten()
                flatphases[c] = phases[c].flatten()

            datalines = ''
            for i in range(latsout.size):
                datalines = "{}{:7.3f}, {:7.3f}:".format(
                    datalines, flatlats[i], flatlons[i])
                for c in constits:
                    datalines = "{}  {:7.3f}".format(datalines, flatamps[c][i])
                datalines = "{}\n".format(datalines)

            msg = (
                "Takin' it for a test drive! Okay!\n" +
                "Here's the tidal amplitudes that were extracted:\n" + 
                "{}\n{}".format(labelline, datalines)    
            )
            print(msg)


    def extract_amp_phase(self, lats, lons, constituents='all', 
        var='z', ocegeo='oce', outname='pyotps_amp_ph', inoutpath=None):
        """Run the OTPS to extract tidal amplitudes and phase at specified
            list of latitudes and longitudes.

            Args:
                lats (array-like): Latitudes of points to extract amp and phase.
                lons (array-like): Longitudes of points to extract amp and phase.
                constituents: Tidal constituents to extract. Either a list of 
                    (lowercase) strings of tidal constituents, or the single
                    string 'all' to extract all available constituents.
                var (str): Tidal variable to extract. Either 'z' (elevation), 
                    'U' (east/west transport), 'V' (north/south transport), 
                    'u' (east/west velocity), 'v' (north/south velocity).
                ocegeo (str): String indicating whether to extract the 
                    'oceancentric' ('oce') or the 'geocentric' ('geo') variables.
                outname (str): Name of the the output text file.

            Returns:
                (lats, lons, amps, phases) with numpy arrays of the input latitude, 
                longitude, amplitude, and phase in the same shape as the input.
                
        """

        # Parameters
        execpath = "{}/extract_HC".format(self.otpspath)
        latlonname = 'latlonfile'
        setupname = 'pyotps_setup'
        validconstits = ['m2', 's2', 'n2', 'k2', 'k1', 'o1', 'p1', 'q1', 'mf', 'mm']
        validvars = ['z', 'U', 'V', 'u', 'v']
        (lats, lons) = (np.array(lats), np.array(lons))

        # Check input
        if lats.shape != lons.shape:
            raise ValueError("Input longitude and latitude are not the same shape!")
        elif constituents is not 'all' and not set(constituents).issubset(validconstits):
            raise ValueError("The parameter 'constituents' must either be 'all' or "
                "a subset of [{}].".format(' '.join(validconstits)))
        elif var not in validvars:
            raise ValueError("The parameter var must be one of {}.".format(' '.join(validvars)))
        elif ocegeo not in ['oce', 'geo']:
            raise ValueError("The parameter 'ocegeo' must either be 'oce' or 'geo'.")

        # Save original dimensions and reshape input lats and lons into 1D list
        dims = lats.shape
        npts = np.array(lats).size
        lats = np.reshape(np.array(lats), (npts,))
        lons = np.reshape(np.array(lons), (npts,))

        # Manage path to store input and output
        if inoutpath is None:
            inoutpath = "{}/pyotps_inout".format(self.otpspath) # Default

        outpath = "{}/{}".format(inoutpath, outname)
        latlonpath = "{}/{}".format(inoutpath, latlonname)
        setuppath = "{}/{}".format(inoutpath, setupname)

        # Generate the list of constituents to be read
        constitstr = ''
        if constituents is not 'all':
            for constit in constituents:
                constitstr = constitstr + "{},".format(constit)

        # Generate a string with a list of the lats and lons to be extracted
        latlonlist = ''
        for i in range(len(lats)):
            latlonlist = "{}{:16.6f} {:16.6f}\n".format(latlonlist, lats[i], lons[i])

        setup = makesetup(self.controlfilepath, latlonpath, var, constitstr,
            'AP', ocegeo, '1', outpath)

        # Manage inout directory and delete output file if it exists
        if not os.path.exists(inoutpath): 
            os.makedirs(inoutpath)

        with open(setuppath, 'w') as setupfile:
            setupfile.write(setup)

        with open(latlonpath, 'w') as latlonfile:
            latlonfile.write(latlonlist)

        # Run_otps
        if os.path.exists(outpath):   os.remove(outpath)

        self.run_tidal_driver(execpath, setuppath)
    
        # Read output, reshape, and return
        constits, lats, lons, amps, phases = read_extract_output(
            outpath, otpstype=self.otpstype)

        for c in amps.keys():
            amps[c] = np.reshape(amps[c], dims)
            phases[c] = np.reshape(phases[c], dims)

        return (constits, 
            np.reshape(lats, dims), np.reshape(lons, dims), amps, phases)



    def run_tidal_driver(self, execpath, setuppath):
        """Call the OTPS executable with a shell command."""

        # Run the model
        with open(setuppath) as setupfile:
            data = subprocess.Popen([execpath], 
                cwd=self.otpspath, stdout=subprocess.PIPE, stdin=setupfile).communicate()



  
    
def makesetup(controlfilepath, latlonfile, var, constitstr, 
        ocegeo, outtype, correctminors, outpath):
    """Returns a string with an OTPS setup.
    
        Args:
            controlfilepath (str): Path to the controlfile.

            latlonfile (str): Path to the file specifying a list of 
                lats and lons to extract data at.

            var (str): The solution variable to extract. Must be one of 
                'z', 'U', 'V', 'u', and 'v'.

            constitstr (str): A string specifying which constituents 
                to extract. A blank string defaults to all available
                constituents. The constituent string must be all 
                lower-case and comma-separated.

            ocegeo (Str): Specify whether to extract ocean ('oce') or
                geocentric ('geo') surface height.

            outtype (str): Type of output which only matters for 
                'extract_HC'. Must be either 'AP' for amplitude/phase 
                or 'RI' for real/imaginary.

            correctminors (str): Set to 1 or 0 to specify whether or 
                not to 'correct for minor constituents'.

            outpath (str): Path to the output file.
                
    """
       
    # Generate the setup
    setup = (
        "{}\n".format(controlfilepath) +
        "{}\n".format(latlonfile) + 
        "{}\n".format(var) + 
        "{}\n".format(constitstr) + 
        "{}\n".format(ocegeo) +
        "{}\n".format(outtype) +
        "{}\n".format(correctminors) +
        "{}\n".format(outpath)
    )

    return setup



def read_extract_output(outpath, otpstype='v1'):
    """"Parse the output of OTPS's extract_HC function and
        return the lats, lons, and constit data in a dict."""

    if otpstype is 'v1':
        lineskip = 0
        headskip = 2
    elif otpstype is 'netcdf':
        lineskip = 0
        headskip = 3
    else:
        raise ValueError("otpstype must be 'v1' or 'netcdf'.")

    # Count lines in output file for preallocation purposes.
    nlines = sum(1 for line in open(outpath) if 'HC extracted' not in line)
    
    if nlines < 3:
        raise RuntimeError("The output file {}".format(outpath) +
                "is invalid or does not have data.")

    npts = nlines-3
    lats = np.empty((npts,), dtype=np.float64)
    lons = np.empty((npts,), dtype=np.float64)

    # Preprocess
    with open(outpath, 'r') as outfile:
        linenum = 0
        irec = 0
        for line in outfile:
            linenum += 1
            if linenum == 3:
                # Parse header and initialize amp and phase arrays
                header = line.split()[headskip:]
                constitdata = np.empty((len(header), npts), dtype=np.float64)

                constits = []
                for i in range(0, len(header), 2):
                    constits.append(header[i][:2])

            elif linenum > 3 and 'HC extracted' not in line:
                # Extract data
                linedata = np.array([float(i) for i in line.split()])

                lats[irec] = linedata[0]
                lons[irec] = linedata[1]
                constitdata[:, irec] = linedata[2:]

                irec += 1

    # Reshape data 
    (amps, phases, i) = ({}, {}, 0)
    for c in constits:
        amps[c]   = constitdata[2*i,   :]
        phases[c] = constitdata[2*i+1, :]
        i += 1

    return constits, lats, lons, amps, phases


