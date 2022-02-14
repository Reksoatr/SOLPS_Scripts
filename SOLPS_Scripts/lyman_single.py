'''
Obtain C-Mod neutral density profiles for a single shot and time interval. 
sciortino, Aug 2020
'''
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import os, copy
#import fit_2D
import pickle as pkl
from scipy.interpolate import interp1d, RectBivariateSpline
import lyman_data
import aurora
from IPython import embed
import mtanh_fitting
from scipy import stats

def get_lya_nn_prof(shot,tmin,tmax, roa_kp, ne, ne_std, Te, Te_std,
                p_ne, p_Te, geqdsk=None, lya_shift=0,
                tomo_inversion=False, zero_pos=0.93, tomo_err=5,
                SOL_exp_decay=True, decay_from_LCFS=False, emiss_min=5e-3): #W/cm^3

    ''' Process Lyman-alpha data for a single shot/time interval. 
        Performs calculation by fitting kinetic profiles and mapping them onto Ly-a data.



    Parameters
    ----------
    shot : int
        CMOD shot number
    tmin : float
        Lower bound of time window
    tmax : float
        Upper bound of time window
    roa_kp : 1D array
        r/a grid
    ne : 1D array
        Electron density in units of :math:`10^{20} m^{-3}`
    ne_std : 1D array
        Uncertainties on electron density in units of :math:`10^{20} m^{-3}`
    Te : 1D array
        Electron temperature in units of :math:`keV`
    Te_std : 1D array
        Uncertainties on electron temperature in units of :math:`keV`
    p_ne : profiletools object
        Electron density object containing experimental data from all loaded diagnostics.
    p_Te : profiletools object
        Electron temperature object containing experimental data from all loaded diagnostics.
    geqdsk : dict
        Dictionary containing processed geqdsk file
    lya_shift : float
        Allows for manual setting of relative shift between Ly-alpha data and kinetic profiles.
        Comes from uncertainty in EFIT position of separatrix.
    tomo_inversion : bool
        If True, use Tomas' tomographic inversion (else use data from tree - Matt's inversion)
    zero_pos : float
        Position at which 0 for brightness data is placed if using Tomas' tomographic inversion
    tomo_err : float
        Systematic error applied to brightness data for Tomas' tomographic inversion (in units of %)
    SOL_exp_decay : bool
        If True, apply an exponential decay (fixed to 1cm length scale) in the region outside of the 
        last radius covered by experimental data.
    decay_from_LCFS : bool
        If True, force exponential decay everywhere outside of the LCFS, ignoring any potential
        data in that region.
    emiss_min : float
        Minimum emissivity to be applied, in W/cm^3.

    Returns
    -------
    res: dict 
        Contains following keys: 

        R : 1D array
            Major radius grid
        roa : 1D array
            r/a grid
        rhop : 1D array
            rhop grid
        ne : 1D array
            Interpolated electron density in units of :math:`cm^{-3}`
        ne_unc : 1D array
            Interpolated electron density uncertaities in units of :math:`cm^{-3}`
        Te : 1D array
            Interpolated electron temperature in units of :math:`eV`
        Te_unc : 1D array
            Interpolated electron temperature uncertainties in units of :math:`eV`
        emiss : 1D array
            Emissivity in units of :math:`W/cm^3`
        emiss_unc : 1D array
            Emissivity uncertainties in units of :math:`W/cm^3`
        nn : 1D array
            Neutral density in units of :math:`cm^{-3}`
        nn_unc : 1D array
            Neutral density uncertainties in units of :math:`cm^{-3}`
        S_ion : 1D array
            Ionization rate in units of :math: `cm^{-3}s^{-1}`
        S_ion_unc : 1D array
            Ionization rate uncertainties in units of :math: `cm^{-3}s^{-1}`

    '''
    time = (tmax+tmin)/2.
    
    ne_decay_len_SOL = 0.01 # m
    Te_decay_len_SOL = 0.01 # m

    # transform coordinates:
    if geqdsk is None:
        geqdsk = lyman_data.get_geqdsk_cmod(
            shot, time*1e3, gfiles_loc = '/home/sciortino/EFIT/lya_gfiles/')

    R_kp = aurora.rad_coord_transform(roa_kp, 'r/a', 'Rmid', geqdsk)
    Rsep = aurora.rad_coord_transform(1.0, 'r/a', 'Rmid', geqdsk)
    rminor = Rsep - geqdsk['RMAXIS'] # minor radius at the midplane
    
    # exponential decays of kinetic profs from last point of experimental data:
    # ne and Te profiles can miss last point depending on filtering?
    max_roa_expt = np.maximum(np.max(p_ne.X[:,0]), np.max(p_Te.X[:,0]))  
    max_rhop_expt = aurora.get_rhop_RZ(geqdsk['RMAXIS']+max_roa_expt*rminor, 0., geqdsk)
    print('Experimental TS data extending to r/a={:.4},rhop={:.4}'.format(max_roa_expt,max_rhop_expt))

    indLCFS = np.argmin(np.abs(R_kp - Rsep))
    ind_max = np.argmin(np.abs(roa_kp - max_roa_expt))
    if SOL_exp_decay and decay_from_LCFS: # exponential decay from the LCFS
        ind_max = indLCFS
        
    ne_std_av = copy.deepcopy(ne_std)
    Te_std_av = copy.deepcopy(Te_std)
    
    if SOL_exp_decay:
        # apply exponential decay outside of region covered by data
        ne_sol = ne[ind_max-1]*np.exp(-(R_kp[ind_max:] - R_kp[ind_max-1])/ne_decay_len_SOL)
        ne_av = np.concatenate((ne[:ind_max], ne_sol))
        Te_sol = Te[ind_max-1]*np.exp(-(R_kp[ind_max:] - R_kp[ind_max-1])/Te_decay_len_SOL)
        Te_av = np.concatenate((Te[:ind_max], Te_sol))

        # set all ne/Te std outside of experimental range to mean of outer-most values
        edge_unc = np.mean(ne_std[ind_max-3:ind_max])
        edge_unc = np.mean(Te_std[ind_max-3:ind_max])
        ne_std_av[ind_max:] = edge_unc if edge_unc<5e13 else 5e13
        Te_std_av[ind_max:] = edge_unc if edge_unc<30e-3 else 30e-3
    else:
        ne_av = copy.deepcopy(ne)
        Te_av = copy.deepcopy(Te)

    # no uncertainties larger than 30 eV outside of LCFS
    Te_std_av[np.logical_and(R_kp>Rsep, Te_std_av>30e-3)] = 30e-3
    #Te_std_av_lcfs = Te_std_av[indLCFS:]
    #Te_std_av_lcfs[Te_std_av_lcfs>30e-3] = 30e-3

    # set ne to cm^-3 and Te in eV
    ne_av *= 1e14
    ne_std_av *= 1e14
    Te_av *= 1e3
    Te_std_av *= 1e3

    # set appropriate minima
    Te_min=10.0    # intended to be a better approximation overall than 3eV
    Te_av[Te_av<Te_min] = Te_min  
    ne_av[ne_av<1e12] = 1e12

    _out = assemble_emiss(shot, tmin, tmax, 'LYMID', tomo_inversion, zero_pos=zero_pos, tomo_err=tomo_err, lya_shift=lya_shift)
    _emiss, _emiss_unc, emiss_R = _out

    # output results in a dictionary, to allow us to add/subtract keys in the future
    res = {'fit':{}}
    out = res['fit']
    res['geqdsk'] = geqdsk  # also save geqdsk and shot number in results dictionary
    res['shot'] = shot
    
    # interpolate kinetic profiles on emissivity radial grid
    out['R'] = np.linspace(np.min(emiss_R), np.max(emiss_R), 200)

    #roa = aurora.rad_coord_transform(out['R'], 'Rmid', 'r/a', geqdsk)
    #out['r/a'] = eq.rho2rho('Rmid', 'r/a', out['R'], time)
    out['r/a'] = (out['R'] - geqdsk['RMAXIS'])/rminor
    
    #rhop = aurora.rad_coord_transform(out['R'], 'Rmid', 'rhop', geqdsk)
    #out['rhop'] = eq.rho2rho('Rmid', 'psinorm', out['R'], time, sqrt=True)
    out['rhop'] = aurora.get_rhop_RZ(out['R'], np.zeros_like(out['R']), geqdsk)

    print('Ly-a data extending r/a={:.3}-{:.3},rhop={:.3}-{:.3f}'.format(np.min(out['r/a']),np.max(out['r/a']),
                                                                         np.min(out['rhop']),np.max(out['rhop'])))

    # interpolate emissivity on a finer radial grid
    out['emiss'] = interp1d(emiss_R, _emiss)(out['R'])
    out['emiss_unc'] = interp1d(emiss_R, _emiss_unc)(out['R'])

    out['emiss'][out['emiss']<emiss_min] = np.nan

    # explore sensitivity of nn_prof within 1 sigma of ne and Te  
    out['ne'] = np.exp(interp1d(roa_kp,np.log(ne_av), bounds_error=False, fill_value=None)(out['r/a']))
    out['Te'] = interp1d(roa_kp,Te_av, bounds_error=False, fill_value=None)(out['r/a'])
  
    ne_up = np.exp(interp1d(roa_kp,np.log(ne_av+ne_std_av), bounds_error=False, fill_value=None)(out['r/a']))
    Te_up = interp1d(roa_kp,Te_av+Te_std_av, bounds_error=False, fill_value=None)(out['r/a'])
    
    ne_min = 1e12
    Te_min = 10

    ne_down = np.exp(interp1d(roa_kp,np.log(np.maximum(ne_av-ne_std_av,ne_min)),
                                   bounds_error=False, fill_value=None)(out['r/a']))
    Te_down = interp1d(roa_kp,np.maximum(Te_av-Te_std_av,Te_min),
                            bounds_error=False, fill_value=None)(out['r/a'])

    # useful for uncertainty estimation on nn/ne:
    out['ne_unc'] = interp1d(roa_kp,ne_std_av, bounds_error=False, fill_value=None)(out['r/a'])
    out['Te_unc'] = interp1d(roa_kp,Te_std_av, bounds_error=False, fill_value=None)(out['r/a'])

    _out = calc_neut_ion(out['ne'], out['Te'], out['emiss'], out['emiss_unc'], out['rhop'],
                        ne_up, ne_down, Te_up, Te_down)
    out['nn'], out['nn_unc'], out['S_ion'], out['S_ion_unc'] = _out

    # estimate uncertainty in nn/S_ion from EFIT
    efit_unc = estimate_efit_unc(out['R'], out['emiss'], lcfs_err=2) # guess err to be 2mm

    # add to other error
    out['nn_unc'] = np.sqrt(out['nn_unc']**2 + (efit_unc*out['nn'])**2)
    out['S_ion_unc'] = np.sqrt(out['S_ion_unc']**2 + (efit_unc*out['S_ion'])**2)

    # eliminate nan values
    mask = ~np.isnan(out['nn'])
    N = len(out['nn'])
    for key in out:
        if len(out[key])==N:
            out[key] = out[key][mask]

    return res
  
  
def get_lya_nn_raw(shot,tmin,tmax, p_ne, p_Te, geqdsk=None, lya_shift=0,
                   tomo_inversion=True, zero_pos=0.93, tomo_err=5,
                   SOL_exp_decay=True, decay_from_LCFS=False, emiss_min=5e-3):

    ''' Process Lyman-alpha data for a single shot/time interval. 
        Performs calculation by mapping Ly-a data onto "raw" TS/SP points.

    Parameters
    ----------
    shot : int
        CMOD shot number
    tmin : float
        Lower bound of time window
    tmax : float
        Upper bound of time window
    p_ne : profiletools object
        Electron density object containing experimental data from all loaded diagnostics.
    p_Te : profiletools object
        Electron temperature object containing experimental data from all loaded diagnostics.
    geqdsk : dict
        Dictionary containing processed geqdsk file
    lya_shift : float
        Allows for manual setting of relative shift between Ly-alpha data and kinetic profiles.
        Comes from uncertainty in EFIT position of separatrix.
    tomo_inversion : bool
        If True, use Tomas' tomographic inversion (else use data from tree - Matt's inversion)
    zero_pos : float
        Position at which 0 for brightness data is placed if using Tomas' tomographic inversion
    tomo_err : float
        Systematic error applied to brightness data for Tomas' tomographic inversion (in units of %)
    SOL_exp_decay : bool
        If True, apply an exponential decay (fixed to 1cm length scale) in the region outside of the 
        last radius covered by experimental data.
    decay_from_LCFS : bool
        If True, force exponential decay everywhere outside of the LCFS, ignoring any potential
        data in that region.
    emiss_min : float
        Minimum emissivity to be applied, in W/cm^3.

    Returns
    -------
    res: dict 
        Contains following keys: 

        R : 1D array
            Major radius grid
        roa : 1D array
            r/a grid
        rhop : 1D array
            rhop grid
        ne : 1D array
            "Raw" electron density in units of :math:`cm^{-3}`
        ne_unc : 1D array
            "Raw" electron density uncertaities in units of :math:`cm^{-3}`
        Te : 1D array
            "Raw" electron temperature in units of :math:`eV`
        Te_unc : 1D array
            "Raw" electron temperature uncertainties in units of :math:`eV`
        emiss : 1D array
            Emissivity in units of :math:`W/cm^3`
        emiss_unc : 1D array
            Emissivity uncertainties in units of :math:`W/cm^3`
        nn : 1D array
            Neutral density in units of :math:`cm^{-3}`
        nn_unc : 1D array
            Neutral density uncertainties in units of :math:`cm^{-3}`
        S_ion : 1D array
            Ionization rate in units of :math:`cm^{-3}s^{-1}`
        S_ion_unc : 1D array
            Ionization rate uncertainties in units of :math:`cm^{-3}s^{-1}`

    '''

    # transform coordinates:
    if geqdsk is None:
        geqdsk = lyman_data.get_geqdsk_cmod(
            shot, time*1e3, gfiles_loc = '/home/sciortino/EFIT/lya_gfiles/')

    Rsep = aurora.rad_coord_transform(1.0, 'r/a', 'Rmid', geqdsk)
    rminor = Rsep - geqdsk['RMAXIS'] # minor radius at the midplane

    # output results in a dictionary, to allow us to add/subtract keys in the future
    res = {'raw':{}}
    out = res['raw']

    ## calculate also from raw data points
    # set ne to cm^-3 and Te in eV

    ne = p_ne.y*1e14
    ne_unc = p_ne.err_y*1e14
    ne_roa = p_ne.X[:,0]

    Te = p_Te.y*1e3
    Te_unc = p_Te.err_y*1e3
    Te_roa = p_Te.X[:,0]
   
    # map onto whichever extends less far into SOL - want to sort points, interpolate, and then unsort them
    map_Te_on_ne = True if np.max(ne_roa) < np.max(Te_roa) else False

    ne_sort = np.argsort(ne_roa)
    ne_unsort = np.argsort(ne_sort)
    Te_sort = np.argsort(Te_roa)
    Te_unsort = np.argsort(Te_sort)


    if map_Te_on_ne:
        out['r/a'] = ne_roa
        out['ne'] = ne
        out['ne_unc'] = ne_unc

        Te_sorted = interp1d(Te_roa[Te_sort], Te[Te_sort], bounds_error=False, fill_value=Te[Te_sort][0])(ne_roa[ne_sort]) # for some reason using "extrapolate" gives error
                                                                                                                            # so fill in with value at smallest r
        Te_unc_sorted = interp1d(Te_roa[Te_sort], Te_unc[Te_sort], bounds_error=False, fill_value=Te_unc[Te_sort][0])(ne_roa[ne_sort])

        # save as unsorted
        out['Te'] = Te_sorted[ne_unsort]
        out['Te_unc'] = Te_unc_sorted[ne_unsort]

    else:
        out['r/a'] = Te_roa
        out['Te'] = Te
        out['Te_unc'] = Te_unc

        ne_sorted = interp1d(ne_roa[ne_sort], ne[ne_sort], bounds_error=False, fill_value=ne[ne_sort][0])(Te_roa[Te_sort]) # for some reason using "extrapolate" gives error
                                                                                                                            # so fill in with value at smallest r
        ne_unc_sorted = interp1d(ne_roa[ne_sort], ne_unc[ne_sort], bounds_error=False, fill_value=ne_unc[ne_sort][0])(Te_roa[Te_sort])

        # same
        out['ne'] = ne_sorted[Te_unsort]
        out['ne_unc'] = ne_unc_sorted[Te_unsort]


    # make sure ne/Te is not negative
    ne_min = 1e12
    Te_min = 10

    out['ne'] = np.maximum(out['ne'], ne_min)
    out['Te'] = np.maximum(out['Te'], Te_min)

    # map kps and emiss to midplane
    out['R'] = geqdsk['RMAXIS'] + out['r/a']*rminor
    out['rhop'] = aurora.get_rhop_RZ(out['R'], np.zeros_like(out['R']), geqdsk)
        
    _out = assemble_emiss(shot, tmin, tmax, 'LYMID', tomo_inversion, zero_pos=zero_pos, tomo_err=tomo_err, lya_shift=lya_shift)
    _emiss, _emiss_unc, emiss_R = _out                                          

    # interpolate emissivity on a finer radial grid
    out['emiss'] = interp1d(emiss_R, _emiss, bounds_error=False)(out['R'])
    out['emiss_unc'] = interp1d(emiss_R, _emiss_unc, bounds_error=False)(out['R'])

    out['emiss'][out['emiss']<emiss_min] = np.nan
    
    # repeat same procedure as above to quantiy errors in nn/ion_rate
    
    # don't need to interpolate these since calculating on kp grid
    ne_up = out['ne']+out['ne_unc']
    Te_up = out['Te']+out['Te_unc'] 

    Te_min=10.0    # intended to be a better approximation overall than 3eV
    ne_down = np.maximum(out['ne']-out['ne_unc'],ne_min)
    Te_down = np.maximum(out['Te']-out['Te_unc'],Te_min)
    
    _out = calc_neut_ion(out['ne'], out['Te'], out['emiss'], out['emiss_unc'], out['rhop'],
                        ne_up, ne_down, Te_up, Te_down)
    out['nn'], out['nn_unc'], out['S_ion'], out['S_ion_unc'] = _out

    # estimate uncertainty in nn/S_ion from EFIT
    efit_unc = estimate_efit_unc(out['R'], out['emiss'], lcfs_err=2) # guess err to be 2mm

    # add to other error
    out['nn_unc'] = np.sqrt(out['nn_unc']**2 + (efit_unc*out['nn'])**2)
    out['S_ion_unc'] = np.sqrt(out['S_ion_unc']**2 + (efit_unc*out['S_ion'])**2)

    return res


def assemble_emiss(shot, tmin, tmax, system, tomo_inversion, zero_pos=0.93, tomo_err=5, lya_shift=0):

    '''Assembles emissivity profile and its uncertainty by getting brightness and potentially 
        inverted emissivity from the tree or performing tomographic inversion.

    Parameters
    -----------
    shot : int
        CMOD shot number
    tmin : float
        Lower bound of time window
    tmax : float
        Upper bound of time window
    system : str
        Name of Lyman-alpha array in tree
    tomo_inversion : Bool
        Flag for determining whether tomographic inversion should be performed or just pulled 
        from tree
    zero_pos : float
        Position at which 0 for brightness data is placed if using Tomas' tomographic inversion
    tomo_err : float
        Systematic error applied to brightness data for Tomas' tomographic inversion (in units of %)
    lya_shift : float
        Allows for manual setting of relative shift between Ly-alpha data and kinetic profiles.

    Returns
    --------
    _emiss_prof : 1D array
        Averaged emissivity profile (pre-interpolation)
    _emiss_unc : 1D array
        Averaged uncertainty in emissivity profile
    edata.R.values : 1D array
        Radial coordinate of emissivity profile
    '''

    # load Lyman emissivity from MDS+
    if tomo_inversion:
        _out = lyman_data.fetch_tomo_emiss(shot, 'LYMID', r_end=zero_pos, sys_err=tomo_err, shift=lya_shift)
        edata, emiss_err = _out
    else:
        edata = lyman_data.fetch_emiss(shot, 'LYMID')

    # time average through 50 ms for each point
    interp = RectBivariateSpline(
        edata.R.values, edata.time.values, edata.values.T, s=0 # no smoothing
        )
    time_vec_av = np.linspace(tmin,tmax,100)
    _emiss_ = interp(edata.R.values,time_vec_av)

    _emiss_prof = np.mean(_emiss_, axis=1)*1e-6  # W/m^3 --> W/cm^3
    
    # tomo inversion returns error from inversion, so no need to recalculate it
    if tomo_inversion:
        interp_err = RectBivariateSpline(
            edata.R.values, edata.time.values, emiss_err.T, s=0 # no smoothing
            )
        _emiss_unc_ = interp_err(edata.R.values,time_vec_av)
        _emiss_unc = np.mean(_emiss_unc_, axis=1)*1e-6
    else:
        _emiss_unc = np.std(_emiss_, axis=1)*1e-6  # W/m^3 --> W/cm^3                                               

    return _emiss_prof, _emiss_unc, edata.R.values


def calc_neut_ion(ne, Te, emiss, emiss_unc, rhop, ne_up, ne_down, Te_up, Te_down):

    '''Calculates neutral density and ionization rate profiles from kinetic and deuterium Ly-alpha
    emissivity measurements.  

    Parameters
    -----------
    ne : 1D array
        Electron density in units of :math:`10^{20} m^{-3}`
    Te : 1D array
        Electron temperature in units of :math:`keV`
    emiss : 1D array
        Emissivity in units of :math: `W/cm^3`
    emiss_unc : 1D array
        Emissivity uncertainties in units of :math: `W/cm^3`
    rhop: 1D array
        rhop grid
    ne_up : 1D array
        Upper limit in electron density using its uncertainty
    ne_down : 1D array
        Lower limit in electron density using its uncertainty
    Te_up : 1D array
        Upper limit in electron temperature using its uncertainty
    Te_down : 1D array
        Lower limit in electron temperature using its uncertainty


    Returns
    --------
    nn : 1D array
        Neutral density in units of :math: `cm^{-3}`
    nn_unc : 1D array
        Neutral density uncertainties in units of :math: `cm^{-3}`
    S_ion : 1D array
        Ionization rate in units of :math: `cm^{-3}s^{-1}`
    S_ion_unc : 1D array
        Ionization rate uncertainties in units of :math: `cm^{-3}s^{-1}`
    '''

    ## neutral density

    # mean profile
    nn,_ = aurora.Lya_to_neut_dens(emiss, ne, Te, plot=False, rhop=rhop)
   
    # testing up and down shifts of ne,Te profiles within uncertainties
    nn_raw_uu,_ = aurora.Lya_to_neut_dens(emiss, ne_up, Te_up, plot=False, rhop=rhop)
    nn_raw_ud,_ = aurora.Lya_to_neut_dens(emiss, ne_up, Te_down, plot=False, rhop=rhop)
    nn_raw_du,_ = aurora.Lya_to_neut_dens(emiss, ne_down, Te_up, plot=False, rhop=rhop)
    nn_raw_dd,_ = aurora.Lya_to_neut_dens(emiss, ne_down, Te_down, plot=False, rhop=rhop)

    nn_raw_low = np.min([nn_raw_uu,nn_raw_ud,nn_raw_du,nn_raw_dd],axis=0)
    nn_raw_high = np.max([nn_raw_uu,nn_raw_ud,nn_raw_du,nn_raw_dd],axis=0)
    nn_raw_unc1 = (nn_raw_high-nn_raw_low)/2.

    # linear propagation of uncertainty:
    nn_unc = nn * np.sqrt((nn_raw_unc1/nn)**2+(emiss_unc/emiss)**2)

    ## ionization rate

    # mean profile 
    S_ion = lyman_data.Lya_to_ion_rate(emiss, ne, Te, rhop=rhop)

    # testing up and down shifts of ne,Te profiles within uncertainties
    ion_raw_uu = lyman_data.Lya_to_ion_rate(emiss, ne_up, Te_up, rhop=rhop)
    ion_raw_ud = lyman_data.Lya_to_ion_rate(emiss, ne_up, Te_down, rhop=rhop)
    ion_raw_du = lyman_data.Lya_to_ion_rate(emiss, ne_down, Te_up, rhop=rhop)
    ion_raw_dd = lyman_data.Lya_to_ion_rate(emiss, ne_down, Te_down, rhop=rhop)

    ion_raw_low = np.min([ion_raw_uu,ion_raw_ud,ion_raw_du,ion_raw_dd],axis=0)
    ion_raw_high = np.max([ion_raw_uu,ion_raw_ud,ion_raw_du,ion_raw_dd],axis=0)
    ion_raw_unc1 = (ion_raw_high-ion_raw_low)/2.

    # linear propagation of uncertainty:
    S_ion_unc = S_ion * np.sqrt((ion_raw_unc1/S_ion)**2+(emiss_unc/emiss)**2)
    
    return nn, nn_unc, S_ion, S_ion_unc


def estimate_efit_unc(R, emiss, lcfs_err=2):
    ''' Estimate error from EFIT propagated into nn/S_ion by calculating gradient
    of emissivity profile and using guess for error in position of LCFS
    '''
    
    # sort, calc deriv, unsort
    inds_sort = np.argsort(R)
    inds_unsort = np.argsort(inds_sort)
  
    # can't have np.diff(R) = 0
    diff_R = np.diff(R[inds_sort])
    no_diff = np.where(diff_R == 0)[0]

    diff_R[no_diff] = 1 # just a placeholder 

    demiss_dr = np.zeros_like(emiss[inds_sort])
    demiss_dr[:-1] = np.diff(emiss[inds_sort])/diff_R
    demiss_dr[-1] = demiss_dr[-2]
    
    prop_err = lcfs_err*1e-3*demiss_dr
    prop_err[no_diff] = prop_err[no_diff+1] # for placeholders, set to next value of error
    
    # unsort prop_err
    prop_err = prop_err[inds_unsort]

    return prop_err



def gaussian_shading(ax, x, y, y_unc, c='k', min_val=0.0):
    ''' Plot profile with uncertainties displayed as a shading whose color intensity represents a 
    gaussian PDF.
    '''
    norm_val = stats.norm.pdf(0)
    
    num=50  # discrete number of shades    
    for ij in np.arange(num):

        # below mean
        ax.fill_between(x,
                        np.maximum(y - 5*y_unc*(ij-1)/num, min_val),
                        np.maximum(y - 5*y_unc*ij/num, min_val),
                        alpha=stats.norm.pdf(5*ij/num)/norm_val,
                        linewidth=0.0,
                        color=c)

    # start looping from 2 to avoid overshading the same region
    for ij in np.arange(2,num):
        # above mean
        ax.fill_between(x, 
                        y + 5*y_unc*(ij-1.)/num,
                        y + 5*y_unc*ij/num,
                        alpha=stats.norm.pdf(5*ij/num)/norm_val,
                        linewidth=0.0,
                        color=c)


def gaussian_3quantiles(ax, a, y, y_unc, c='k'):
    '''Alternative style of plotting where we show the 1-99, 10-90 and 25-75 quantiles of the
    probability distribution, as in Francesco's impurity transport inferences. 

    This may be over-complicated and somewhat redundant if the uncertainties are actually gaussian
    (as assumed in this function). It makes more sense when distributions are arbitrary...
    but this function could be useful to keep consistency of plotting styles.
    '''
    alp=0.6
    ax[0].fill_between(x,
                       y+stats.norm.ppf(0.25)*ff*y_unc/y,
                       y+stats.norm.ppf(0.75)*ff*y_unc/y,
                       alpha=alp, color=c)
    ax[0].fill_between(x,
                       y+stats.norm.ppf(0.10)*ff*y_unc/y,
                       y+stats.norm.ppf(0.90)*ff*y_unc/y,
                       alpha=alp/2, color=c)    
    ax[0].fill_between(x,
                       y+stats.norm.ppf(0.01)*ff*y_unc/y,
                       y+stats.norm.ppf(0.99)*ff*y_unc/y,
                       alpha=alp/3, color=c)
    
    

def plot_nn_by_ne(res):
    '''Plot profiles of neutral density and its ratio to the electron density.
    '''
    ff = 1./np.log(10.)

    nn_by_ne = res['nn']/res['ne']
    # ne uncertainty in the SOL also goes to ne<0.... ignore it
    nn_by_ne_unc = np.sqrt((res['nn_unc']/res['ne'])**2) #+(nn_prof/ne_prof**2)**2*ne_prof_unc**2)  
    nn_by_ne[nn_by_ne<1e-8] = 1e-8
        
    # plot only n0 and n0/ne
    fig, ax = plt.subplots(1,2, figsize=(12,6), sharex=True)
    gaussian_shading(ax[0], res['rhop'], np.log10(res['nn']), ff*res['nn_unc']/res['nn'], c='k')
    gaussian_shading(ax[1], res['rhop'], np.log10(nn_by_ne), ff*np.abs(nn_by_ne_unc/nn_by_ne), c='k', min_val=-80.)

    ax[0].set_ylabel(r'$log_{10}(n_n$ [$cm^{-3}$])')
    ax[1].set_ylabel(r'$log_{10}(n_n/n_e)$')
    ax[1].set_xlabel(r'$\rho_p$')
    ax[0].set_xlabel(r'$\rho_p$')
    fig.suptitle(f'C-Mod shot {shot}')
    fig.tight_layout()


def plot_lyman_overview(res, overplot_raw=False, emiss_min = 5e-3, Te_min=10., num_SP=None):
    '''Plot overview of Lyman-alpha results, showing emissivity, electron density and temperature, 
    neutral density and ionization profile

    Parameters:
    res : dict
        Dictionary containing the processed result of :py:fun:`get_lya_nn_prof`.
    overplot_raw : bool
        If True, overplot data points mapped from TS and probe data all the way to the ionization
        rate. This isn't "raw data", but rather a mapping from the raw data locations.
    emiss_min : float
        Minimum Ly-a emissivity in W/cm^3 that was imposed during data processing. This is only
        used in this function to show it on the plot.
    '''

    fit = res['fit']
    if overplot_raw: raw = res['raw']
    geqdsk = res['geqdsk']
    Rsep = aurora.rad_coord_transform(1.0, 'r/a', 'Rmid', geqdsk)
    rminor = Rsep - geqdsk['RMAXIS'] # minor radius at the midplane

    ff = 1./np.log(10.)
   
    # complete analysis layout:
    fit['ne'][fit['ne']<1e10] = 1e10
    fit['nn'][fit['nn']<1.] = 1.0



    
    # -------------------------------

    fig,ax = plt.subplots(4,1, figsize=(8,12), sharex=True)

    set_xlim=False
    # set x-lims specifically based on where we have good Ly-a and TS coverage for FS shots
    if shot==1080416025:
        ax[0].set_xlim([0.955, 1.005])
    elif shot==1100308004:
        ax[0].set_xlim([0.96, 1.01])
    elif shot==1100305023:
        ax[0].set_xlim([0.94, 1.01]) # can go much further, very good Thomson data far out
    else:
        set_xlim=True
        pass
    
    # top plot is emissivity
    ax[0].axhline(emiss_min, c='r', ls='--')
    gaussian_shading(ax[0], fit['rhop'], fit['emiss'], fit['emiss_unc'], c='k', min_val=0.0)

    
    # Plot ne
    gaussian_shading(ax[1], fit['rhop'], fit['ne'], fit['ne_unc'], c='b', min_val = 1e11)
    if overplot_raw:
        ne_ts_mask = raw['r/a'] > np.min(fit['rhop']); ne_ts_mask[-num_SP['ne']:] = False
        ne_sp_mask = raw['r/a'] > np.min(fit['rhop']); ne_sp_mask[:-num_SP['ne']] = False
        p_ne_R = geqdsk['RMAXIS']+raw['r/a']*rminor
        p_ne_rhop = aurora.get_rhop_RZ(p_ne_R, np.zeros_like(p_ne_R), geqdsk)
        ax[1].errorbar(p_ne_rhop[ne_sp_mask], raw['ne'][ne_sp_mask], raw['ne_unc'][ne_sp_mask], color='blue', fmt='x')
        ax[1].errorbar(p_ne_rhop[ne_ts_mask], raw['ne'][ne_ts_mask], raw['ne_unc'][ne_ts_mask], color='blue', fmt='.')

    # find appropriate max of y scale
    _ne_range = raw['ne'][p_ne_rhop>np.min(fit['rhop'])]
    ax[1].set_ylim([0, np.max(_ne_range)+3e13])  # 3e13 cm^-3 above max
    
    # plot Te on twin (right) axis
    ax2 = ax[1].twinx()  # instantiate a second axis that shafit the same x-axis
    gaussian_shading(ax2, fit['rhop'], fit['Te'], fit['Te_unc'], c='r', min_val = Te_min)
    ax2.tick_params(axis='y', labelcolor='r')
    if overplot_raw:
        Te_ts_mask = raw['r/a'] > np.min(fit['rhop']); Te_ts_mask[-num_SP['Te']:] = False
        Te_sp_mask = raw['r/a'] > np.min(fit['rhop']); Te_sp_mask[:-num_SP['Te']] = False
        p_Te_R = geqdsk['RMAXIS']+raw['r/a']*rminor
        p_Te_rhop = aurora.get_rhop_RZ(p_Te_R, np.zeros_like(p_Te_R), geqdsk)
        ax2.errorbar(p_Te_rhop[Te_sp_mask], raw['Te'][Te_sp_mask], raw['Te_unc'][Te_sp_mask], color='red', fmt='x')
        ax2.errorbar(p_Te_rhop[Te_ts_mask], raw['Te'][Te_ts_mask], raw['Te_unc'][Te_ts_mask], color='red', fmt='.')

    # find appropriate max of y scale
    _Te_range = raw['Te'][p_Te_rhop>np.min(fit['rhop'])]
    ax2.set_ylim([0, np.max(_Te_range)+100])  # 100 eV above max

    
    # Plot neutral density
    if overplot_raw: 

        ax[2].errorbar(raw['rhop'][ne_sp_mask], np.log10(raw['nn'][ne_sp_mask]),
                       ff*raw['nn_unc'][ne_sp_mask]/raw['nn'][ne_sp_mask], color='black', fmt='x')
        ax[2].errorbar(raw['rhop'][ne_ts_mask], np.log10(raw['nn'][ne_ts_mask]),
                       ff*raw['nn_unc'][ne_ts_mask]/raw['nn'][ne_ts_mask], color='black', fmt='o')
        #ax[2].errorbar(raw['rhop'], np.log10(raw['nn']),
        #               ff*raw['nn_unc']/raw['nn'], color='salmon', fmt='.')
        
    mask = np.logical_and(fit['rhop']>ax[0].get_xlim()[0] , fit['rhop']<ax[0].get_xlim()[1])
    gaussian_shading(ax[2], fit['rhop'][mask], np.log10(fit['nn'][mask]), ###put mask in
                     ff*fit['nn_unc'][mask]/fit['nn'][mask], c='k', min_val = 5)

    #ax[2].set_ylim([8,12]) # reasonable range, may need to be adapted
    
    # Plot neutral density ratio with electron density
    #gaussian_shading(ax[3], fit['rhop'], nn_by_ne_prof, nn_by_ne_prof_unc, c='k', min_val = -15)


    # Ionization rate
    if overplot_raw: 

        ax[3].errorbar(raw['rhop'][ne_sp_mask], np.log10(raw['S_ion'][ne_sp_mask]),
                       ff*raw['S_ion_unc'][ne_sp_mask]/raw['S_ion'][ne_sp_mask], color='black', fmt='x')
        ax[3].errorbar(raw['rhop'][ne_ts_mask], np.log10(raw['S_ion'][ne_ts_mask]),
                       ff*raw['S_ion_unc'][ne_ts_mask]/raw['S_ion'][ne_ts_mask], color='black', fmt='o')
        #ax[3].errorbar(raw['rhop'], np.log10(raw['S_ion']),
        #               ff*raw['S_ion_unc']/raw['S_ion'], color='salmon', fmt='.')

    mask = np.logical_and(fit['rhop']>ax[0].get_xlim()[0] , fit['rhop']<ax[0].get_xlim()[1])      
    gaussian_shading(ax[3], fit['rhop'][mask], np.log10(fit['S_ion'][mask]),
                     ff*fit['S_ion_unc'][mask]/fit['S_ion'][mask], c='k', min_val = 0.0)

    
    ax[0].set_ylabel('Ly-a emiss [$W/cm^3$]')
    ax[1].set_ylabel(r'$n_e$ [$cm^{-3}$]', color='b')
    ax[1].tick_params(axis='y', colors='b')
    ax2.set_ylabel(r'$T_e$ [$eV$]', color='r')
    ax[2].set_ylabel(r'$log_{10}(n_n$ [$cm^{-3}$])')
    #ax[2].set_ylim([5,15])
    #ax[3].set_ylabel(r'$log_{10}(n_n/n_e)$')
    ax[3].set_ylabel(r'$log_{10}(S_{ion}$ [$cm^{-3}s^{-1}$])')


    # limit radial range to where we have some kinetic profile data and some lya data
    te_mask = np.logical_and(p_Te_rhop>np.min(fit['rhop']), p_Te_rhop<np.max(fit['rhop']))
    ne_mask = np.logical_and(p_ne_rhop>np.min(fit['rhop']), p_ne_rhop<np.max(fit['rhop']))

    if set_xlim:
        xmin_val = np.min([np.min(p_ne_rhop[ne_mask]),np.min(p_Te_rhop[te_mask]), np.min(fit['rhop'])]) #-0.01
        xmax_val = np.max([np.max(p_ne_rhop[ne_mask]),np.max(p_Te_rhop[te_mask])])  # lya data always goes further out
        ax[0].set_xlim([xmin_val, xmax_val])
    
    ax[-1].set_xlabel(r'$\rho_p$')

    #ax[0].set_xlim([np.min(fit['rhop'][~np.isnan(fit['emiss'])]),np.max(fit['rhop'][~np.isnan(fit['emiss'])])])
   
    fig.suptitle(f'C-Mod shot {res["shot"]}')
    fig.tight_layout()




    

if __name__=='__main__':
    # I-mode:
    #shot=1080416025   # 24/25
    #tmin=0.8 #1.2 # gas puffs around 1.0s
    #tmax=1.0 #1.2

    # L-mode:
    #shot=1100308004
    #tmin=0.9 #0.7
    #tmax=1.0 #1.1 #1.4

    # EDA H-mode:
    shot=1100305023
    tmin=0.85
    tmax=1.3
    
    # Other L-modes:
    #shot=1080110005
    #tmin=1.0
    #tmax=1.3

    # L-mode:
    #shot=1070511002
    #tmin=0.7
    #tmax=1.4

    # great probe data:
    #shot = 1070511002
    #tmin = 0.9
    #tmax = 1.1

    # case with good probe data
    #shot = 1070511010
    #tmin = 0.9
    #tmax = 1.1

    # snypes shots
    #shot = 1080430003
    #tmin=0.9
    #tmax=1.1

    # lipschultz
    #shot=1080416026   # 24/25
    #tmin=1.1 #1.2 # gas puffs around 1.0s
    #tmax=1.2

    # marmar -- very nice
    #shot=1080416028   # 24/25
    #tmin=0.8 #1.1 #1.2 # gas puffs around 1.0s
    #tmax=0.9 #1.2

    # test
    #shot = 1080124015
    #tmin = 0.6
    #tmax = 0.7

    # for Richard
    #shot = 1080416025 
    #tmin = 0.9
    #tmax = 1.1

    # for Abhi
    #shot = 1120711021
    #tmin = 0.6
    #tmax = 1.2
    
    # SOLPS
    #shot = 1120917011
    #tmin = 0.6
    #tmax = 1.2

    # more shots
    #shot = 1070614006
    #tmin = 1.0
    #tmax = 1.2
   
    # pedestal transport shots
    #shot = 1080411025
    #tmin = 0.5
    #tmax = 1.5

    ############
    emiss_min = 5e-3  # W/cm^3
    ne_min = 1e12
    Te_min = 10. #eV
    force_to_zero = True if shot not in [1080416025, 1100308004, 1100305023] else False
    
    #gfiles_loc = '/home/sciortino/EFIT/gfiles/' #'.'
    gfiles_loc = '.'

    geqdsk = lyman_data.get_geqdsk_cmod(shot,(tmin+tmax)/2.*1e3, gfiles_loc=gfiles_loc)
    
    kp_out = lyman_data.get_cmod_kin_profs(shot, tmin, tmax, geqdsk = geqdsk, probes=['A','F'],
                                           apply_final_sep_stretch=True, force_to_zero=force_to_zero) 
    roa_kp, ne, ne_std, Te, Te_std, p_ne, p_Te, num_SP = kp_out

    time = (tmin+tmax)/2

    # gather Lyman-alpha profiles:
    res = get_lya_nn_prof(shot,tmin,tmax,roa_kp, ne, ne_std, Te, Te_std,
                          p_ne,p_Te, geqdsk=geqdsk, tomo_inversion=True,
                          SOL_exp_decay=False)
     
    # gather Lyman-alpha raw mapped data:
    res_raw = get_lya_nn_raw(shot,tmin,tmax, p_ne,p_Te, geqdsk=geqdsk, tomo_inversion=True,
                             SOL_exp_decay=False)
                   
    res.update(res_raw)
    
    # plot analysis results
    plot_nn_by_ne(res['fit'])
    
    plot_lyman_overview(res, overplot_raw=True, emiss_min=emiss_min, Te_min=Te_min, num_SP=num_SP)

        
    # store profiles in a pickle file (warning: changing this data structure will break things)
    fit = res['fit']
    out2 = [fit['rhop'],fit['r/a'],fit['R'],fit['nn'],fit['nn_unc'],fit['S_ion'],
            fit['S_ion_unc'],fit['ne'],fit['ne_unc'],fit['Te'],fit['Te_unc'],
            fit['emiss'], fit['emiss_unc']]
    
    with open(f'Dicts/lyman_data_{shot}.pkl','wb') as f:
        pkl.dump(out2,f)

    
