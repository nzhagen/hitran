from __future__ import print_function, division
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import glob

__all__ = ['lorentzian_profile', 'read_hitran2012_parfile', 'translate_molecule_identifier',
           'get_molecule_identifier', 'calculate_hitran_xsec', 'downsample_spectrum', 'draw_block_spectrum']

## ======================================================
def lorentzian_profile(kappa, S, gamma, kappa0):
    '''
    Calculate a Lorentzian absorption profile.

    Parameters
    ----------
    kappa : ndarray
        The array of wavenumbers at which to sample the profile.
    S : float
        The absorption line "strength" (the integral of the entire line is equal to S).
    gamma : float
        The linewidth parameter of the profile.
    kappa0 : float
        The center position of the profile (in wavenumbers).

    Returns
    -------
    L : ndarray
        The sampled absorption profile.
    '''
    L = (S / pi) * gamma / ((kappa - kappa0)**2 + gamma**2)
    return(L)

## ======================================================
def read_hitran2012_parfile(filename):
    '''
    Given a HITRAN2012-format text file, read in the parameters of the molecular absorption features.

    Parameters
    ----------
    filename : str
        The filename to read in.

    Return
    ------
    data : dict
        The dictionary of HITRAN data for the molecule.
    '''

    if filename.endswith('.zip'):
        import zipfile
        zip = zipfile.ZipFile(filename, 'r')
        (object_name, ext) = os.path.splitext(os.path.basename(filename))
        print(object_name, ext)
        filehandle = zip.read(object_name).splitlines()
    else:
        filehandle = open(filename, 'r')

    data = {'M':[],               ## molecule identification number
            'I':[],               ## isotope number
            'linecenter':[],      ## line center wavenumber (in cm^{-1})
            'S':[],               ## line strength, in cm^{-1} / (molecule m^{-2})
            'Acoeff':[],          ## Einstein A coefficient (in s^{-1})
            'gamma-air':[],       ## line HWHM for air-broadening
            'gamma-self':[],      ## line HWHM for self-emission-broadening
            'Epp':[],             ## energy of lower transition level (in cm^{-1})
            'N':[],               ## temperature-dependent exponent for "gamma-air"
            'delta':[],           ## air-pressure shift, in cm^{-1} / atm
            'Vp':[],              ## upper-state "global" quanta index
            'Vpp':[],             ## lower-state "global" quanta index
            'Qp':[],              ## upper-state "local" quanta index
            'Qpp':[],             ## lower-state "local" quanta index
            'Ierr':[],            ## uncertainty indices
            'Iref':[],            ## reference indices
            'flag':[],            ## flag
            'gp':[],              ## statistical weight of the upper state
            'gpp':[]}             ## statistical weight of the lower state

    print('Reading "' + filename + '" ...')

    for line in filehandle:
        if (len(line) < 160):
            raise ImportError, 'The imported file ("' + filename + '") does not appear to be a HITRAN2012-format data file.'

        data['M'].append(uint(line[0:2]))
        data['I'].append(uint(line[2]))
        data['linecenter'].append(float64(line[3:15]))
        data['S'].append(float64(line[15:25]))
        data['Acoeff'].append(float64(line[25:35]))
        data['gamma-air'].append(float64(line[35:40]))
        data['gamma-self'].append(float64(line[40:45]))
        data['Epp'].append(float64(line[45:55]))
        data['N'].append(float64(line[55:59]))
        data['delta'].append(float64(line[59:67]))
        data['Vp'].append(line[67:82])
        data['Vpp'].append(line[82:97])
        data['Qp'].append(line[97:112])
        data['Qpp'].append(line[112:127])
        data['Ierr'].append(line[127:133])
        data['Iref'].append(line[133:145])
        data['flag'].append(line[145])
        data['gp'].append(line[146:153])
        data['gpp'].append(line[153:160])

    if filename.endswith('.zip'):
        zip.close()
    else:
        filehandle.close()

    for key in data:
        data[key] = array(data[key])

    return(data)

## ======================================================
def translate_molecule_identifier(M):
    '''
    For a given input molecule identifier number, return the corresponding molecular formula.

    Parameters
    ----------
    M : int
        The HITRAN molecule identifier number.

    Returns
    -------
    molecular_formula : str
        The string describing the molecule.
    '''

    trans = { '1':'H2O',    '2':'CO2',   '3':'O3',      '4':'N2O',   '5':'CO',    '6':'CH4',   '7':'O2',     '8':'NO',
              '9':'SO2',   '10':'NO2',  '11':'NH3',    '12':'HNO3', '13':'OH',   '14':'HF',   '15':'HCl',   '16':'HBr',
             '17':'HI',    '18':'ClO',  '19':'OCS',    '20':'H2CO', '21':'HOCl', '22':'N2',   '23':'HCN',   '24':'CH3Cl',
             '25':'H2O2',  '26':'C2H2', '27':'C2H6',   '28':'PH3',  '29':'COF2', '30':'SF6',  '31':'H2S',   '32':'HCOOH',
             '33':'HO2',   '34':'O',    '35':'ClONO2', '36':'NO+',  '37':'HOBr', '38':'C2H4', '39':'CH3OH', '40':'CH3Br',
             '41':'CH3CN', '42':'CF4',  '43':'C4H2',   '44':'HC3N', '45':'H2',   '46':'CS',   '47':'SO3'}
    return(trans[str(M)])

## ======================================================
def get_molecule_identifier(molecule_name):
    '''
    For a given input molecular formula, return the corresponding HITRAN molecule identifier number.

    Parameters
    ----------
    molecular_formula : str
        The string describing the molecule.

    Returns
    -------
    M : int
        The HITRAN molecular identified number.
    '''

    trans = { '1':'H2O',    '2':'CO2',   '3':'O3',      '4':'N2O',   '5':'CO',    '6':'CH4',   '7':'O2',     '8':'NO',
              '9':'SO2',   '10':'NO2',  '11':'NH3',    '12':'HNO3', '13':'OH',   '14':'HF',   '15':'HCl',   '16':'HBr',
             '17':'HI',    '18':'ClO',  '19':'OCS',    '20':'H2CO', '21':'HOCl', '22':'N2',   '23':'HCN',   '24':'CH3Cl',
             '25':'H2O2',  '26':'C2H2', '27':'C2H6',   '28':'PH3',  '29':'COF2', '30':'SF6',  '31':'H2S',   '32':'HCOOH',
             '33':'HO2',   '34':'O',    '35':'ClONO2', '36':'NO+',  '37':'HOBr', '38':'C2H4', '39':'CH3OH', '40':'CH3Br',
             '41':'CH3CN', '42':'CF4',  '43':'C4H2',   '44':'HC3N', '45':'H2',   '46':'CS',   '47':'SO3'}
    ## Invert the dictionary.
    trans = {v:k for k,v in trans.items()}
    return(int(trans[molecule_name]))

## ======================================================
def calculate_hitran_xsec(data, wavemin=None, wavemax=None, npts=20001, units='m^2'):
    '''
    Given the HITRAN data (line centers and line strengths) for a molecule, digitize the result into a spectrum of
    absorption cross-section in units of cm^2.

    Parameters
    ----------
    data : dict of ndarrays
        The HITRAN data corresponding to a given molecule.
    wavemin : float, optional
        The minimum wavelength os the spectral region of interest.
    wavemax : float, optional
        The maximum wavelength os the spectral region of interest.
    units : str, optional
        A string describing in what units of the output cross-section should be given in. Choices available are:
        {'cm^2/mole', 'cm^2.ppm', 'm^2/mole', 'm^2.ppm', 'm^2', cm^2}.

    Returns
    -------
    waves : ndarray
        The wavelengths at which the cross-section data is evaluated.
    xsec : array_like
        The mean absorption cross-section (in cm^2) per molecule, evaluated at the wavelengths given by input `waves`.
    '''

    if (wavemin == None):
        wavemin = amin(10000.0 / data['linecenter']) - 0.1
    if (wavemax == None):
        wavemax = amax(10000.0 / data['linecenter']) + 0.1

    ## First step: remove any data points that do not correspond to the primary isotope. (If we want to use isotopes,
    ## then we need to figure out mixing ratios.) For most HITRAN gases, the primary isotope is about 99% of the total
    ## atmospheric composition.
    okay = (data['I'] == 1)
    linecenters = array(data['linecenter'][okay])       ## line centers in wavenumbers
    linestrengths = array(data['S'][okay])
    linewidths = array(data['gamma-air'][okay])
    nlines = alen(linecenters)

    ## Convert the wavelengths (um) to wavenumbers (cm^{-1}). Create a spectrum linearly sampled in wavenumber (and
    ## thus nonuniformly sampled in wavelength).
    wavenumbers = linspace(10000.0/wavemax, 10000.0/wavemin, npts)
    waves = 10000.0 / wavenumbers
    xsec = zeros_like(wavenumbers)

    ## Define the list of channel boundary wavelengths.
    dk = wavenumbers[1] - wavenumbers[0]

    for i in arange(nlines):
        linecenter = linecenters[i]
        linestrength = linestrengths[i]
        linewidth = linewidths[i]

        ## If the spectral line is well outside our region of interest, then ignore it.
        if (linecenter < amin(wavenumbers-0.5)):
            continue
        elif (linecenter > amax(wavenumbers+0.5)):
            continue

        ## Note: the quantity sum(L * dk) should sum to "S"!
        L = lorentzian_profile(wavenumbers, linestrength, linewidth, linecenter)
        xsec += L

    if units.endswith('/mole'):
        xsec = xsec * 6.022E23
    elif units.endswith('.ppm'):
        xsec = xsec * 2.686E19

    if units.startswith('cm^2'):
        pass
    elif units.startswith('m^2'):
        xsec = xsec / 10000.0

    return(waves, xsec)

## ======================================================
def downsample_spectrum(waves, spectrum, downsampled_waves=None, downsampled_channel_boundaries=None):
    '''
    (Right now, we assume uniformly sampled wavelengths/wavenumbers.)
    '''

    nwaves = alen(waves)

    ## If it is not already defined, make the list of channel boundary wavelengths.
    if (downsampled_waves != None) and (downsampled_channel_boundaries == None):
        dw = downsampled_waves[1] - downsampled_waves[0]
        downsampled_channel_boundaries = append(amin(downsampled_waves)-(dw/2.0), downsampled_waves+(dw/2.0))
    elif (downsampled_waves == None) and (downsampled_channel_boundaries == None):
        raise ValueError, 'Either "downsampled_waves" or "downsampled_channel_boundaries" is required as an input.'

    ## Generate the channel basis functions used to represent the low-resolution spectral channels in terms
    ## of the high-resolution data.
    nchannels = alen(downsampled_channel_boundaries) - 1
    #print('downwaves=', downwaves)
    #print('downsampled_channel_boundaries=', downsampled_channel_boundaries)
    downspectrum = zeros(nchannels)

    ## From the list of channel boundary wavelengths, generate the channel basis functions.
    for n in arange(nchannels):
        okay = (waves > downsampled_channel_boundaries[n]) * (waves <= downsampled_channel_boundaries[n+1])
        downspectrum[n] = nanmean(spectrum[okay])

    return(downsampled_channel_boundaries, downspectrum)

## =================================================================================================
def draw_block_spectrum(channel_boundaries, spectrum, newfigure=True, title=None, **kwargs):
    '''
    Draw a plot where the spectral channels are nonuniform in width and shaped by histogram-like rectangles.

    Parameters
    ----------
    channel_boundaries : array_like
        A vector of length Nw+1 giving the wavelengths defining the boundaries of the Nw spectral channels.
    spectrum : array_like
        A Nw length vector.
    newfigure : bool
        Whether or not to call matplotlib's `figure()` function.
    title : str
        The plot title.
    **kwargs : any
        Any keyword arguments that can be used by matplotlib's `plot()` function.
    '''

    channel_boundaries = asarray(channel_boundaries)
    spectrum = asarray(spectrum)
    assert (alen(channel_boundaries) == 1 + alen(spectrum)), 'Input "channel_boundaries" must have length 1 more than input "spectrum".'

    cb = channel_boundaries
    nchannels = alen(cb) - 1

    x = []
    y = []
    x.append(cb[0])
    y.append(0.0)

    for n in arange(nchannels):
        x.append(cb[n])
        x.append(cb[n+1])
        y.append(spectrum[n])
        y.append(spectrum[n])

    x.append(cb[-1])
    y.append(0.0)

    if newfigure:
        fig = plt.figure()
        if (title != None): fig.canvas.set_window_title(title)
        xmin = amin(x)
        xmax = amax(x)
        xptp = xmax - xmin
        xmean = 0.5 * (xmin + xmax)
        xlo = xmean - 0.55 * xptp
        xhi = xmean + 0.55 * xptp

        ymin = amin(y)
        ymax = amax(y)
        yptp = ymax - ymin
        ymean = 0.5 * (ymin + ymax)
        ylo = ymean - 0.55 * yptp
        yhi = ymean + 0.55 * yptp
    else:
        ## First grab the existing plot limits. If these were previously determined by draw_block_spectrum(),
        ## then we need only rescale the plot range by (0.5/0.55) to get to the original data limits.
        ## Appending these to the current data vector, we can update the plot limits using both old and new
        ## data. First grab the existing limits.
        (x0,x1,y0,y1) = plt.axis()

        x0mean = 0.5 * (x0 + x1)
        x0ptp = (0.5 / 0.55) * (x1 - x0)
        x0min = x0mean - 0.5 * x0ptp
        x0max = x0mean + 0.5 * x0ptp

        y0mean = 0.5 * (y0 + y1)
        y0ptp = (0.5 / 0.55) * (y1 - y0)
        y0min = y0mean - 0.5 * y0ptp
        y0max = y0mean + 0.5 * y0ptp

        ## Next determine the new plot range using the old and new data limits.
        xmin = amin(append(array(x), x0min))
        xmax = amax(append(array(x), x0max))
        xptp = xmax - xmin
        xmean = 0.5 * (xmin + xmax)
        xlo = xmean - 0.55 * xptp
        xhi = xmean + 0.55 * xptp

        ymin = amin(append(array(y), y0min))
        ymax = amax(append(array(y), y0max))
        yptp = ymax - ymin
        ymean = 0.5 * (ymin + ymax)
        ylo = ymean - 0.55 * yptp
        yhi = ymean + 0.55 * yptp

    plt.plot(x, y, **kwargs)
    if (title != None): plt.title(title)
    if newfigure:
        plt.axis([xlo,xhi,ylo,yhi])

    return


## ==================================================================================================
## ==================================================================================================

if (__name__ == "__main__"):
    molecule = 'H2O'       ## water
    #molecule = 'CO2'       ## carbon dioxide
    #molecule = 'NH3'       ## ammonia
    #molecule = 'SO2'       ## sulfur dioxide
    #molecule = 'CH4'       ## methane
    #molecule = 'H2S'       ## hydrogen sulfide
    #molecule = 'O3'        ## ozone
    #molecule = 'C2H6'      ## ethane

    #units = 'm^2/mole'
    #units = 'm^2.ppm'
    #units = 'cm^2/mole'
    #units = 'cm^2.ppm'
    units = 'm^2'
    #units = 'cm^2'

    wavemin = 1.0
    wavemax = 14.0
    #wavemin = 1.0
    #wavemax = 1.7

    show_downsampled_spectrum = (wavemax - wavemin) > 0.8

    M = get_molecule_identifier(molecule)
    filenames = glob.glob('./par/%02i_hit12.*' % M)
    filename = filenames[0]

    filename = os.path.normpath(os.path.abspath(filename))
    data = read_hitran2012_parfile(filename)

    ## Next set up for calculating the absorption cross-section, given the HITRAN database's values for the
    ## line centers and line strengths.
    nlines = len(data['S'])
    print('Found %i lines' % nlines)
    print('Calculating the absorption cross-section spectrum ...')
    (waves, xsec) = calculate_hitran_xsec(data, wavemin, wavemax, units=units)

    fig = plt.figure()
    fig.canvas.set_window_title(molecule)
    plt.semilogy(waves, xsec, 'k-')
    plt.title(molecule)
    plt.ylabel('Cross-section (' + units + ')')
    plt.xlabel('wavelength (um)')

    if show_downsampled_spectrum:
        nchannels = int((amax(waves) - amin(waves)) / 0.2)
        downwaves = linspace(amin(waves),amax(waves),nchannels)
        (downsampled_channel_boundaries, downspectrum) = downsample_spectrum(waves, xsec, downwaves)
        draw_block_spectrum(downsampled_channel_boundaries, downspectrum, linewidth=3.0, color='red',
                            newfigure=False, zorder=0)

    xmin = amin(waves)
    xmax = amax(waves)
    xptp = xmax - xmin
    xmean = 0.5 * (xmin + xmax)
    xlo = xmean - 0.55 * xptp
    xhi = xmean + 0.55 * xptp

    ymin = amin(xsec)
    ymax = amax(xsec)
    yptp = ymax - ymin
    ymean = 0.5 * (ymin + ymax)
    ylo = ymean - 0.55 * yptp
    yhi = ymean + 0.55 * yptp

    plt.xlim([xlo,xhi])

    plt.show()
