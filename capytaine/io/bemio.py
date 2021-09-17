import logging

import numpy as np
import pandas as pd
import xarray as xr
from scipy.optimize import newton

from capytaine import __version__


LOG = logging.getLogger(__name__)

#######################
#  Import from Bemio  #
#######################

def dataframe_from_bemio(bemio_obj, wavenumber, wavelength):
    """Transform a :class:`bemio.data_structures.bem.HydrodynamicData` into a
        :class:`pandas.DataFrame`.

        Parameters
        ----------
        bemio_obj: Bemio data_stuctures.bem.HydrodynamicData class
            Loaded NEMOH, AQWA, or WAMIT data created using `bemio.io.nemoh.read`,
            `bemio.io.aqwa.read`, or `bemio.io.wamit.read` functions, respectively.
        wavenumber: bool
            If True, the coordinate 'wavenumber' will be added to the output dataset.
        wavelength: bool
            If True, the coordinate 'wavelength' will be added to the output dataset.
        """


    dofs = np.array(['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'])
    df = pd.DataFrame()
    for i in range(bemio_obj.body[0].num_bodies):
        difr_dict = []
        rad_dict = []

        rho = bemio_obj.body[0].rho
        g = bemio_obj.body[0].g

        if bemio_obj.body[i].water_depth == 'infinite':
            bemio_obj.body[i].water_depth = np.infty

        if bemio_obj.body[i].bem_code == 'WAMIT': # WAMIT coefficients need to be dimensionalized
            from_wamit = True

        for omega_idx, omega in enumerate(np.sort(bemio_obj.body[i].w)):

            # DiffractionProblem variable equivalents
            for dir_idx, dir in enumerate(bemio_obj.body[i].wave_dir):
                temp_dict = {}
                temp_dict['body_name'] = bemio_obj.body[i].name
                temp_dict['water_depth'] = bemio_obj.body[i].water_depth
                temp_dict['omega'] = omega
                temp_dict['rho'] = rho
                temp_dict['g'] = g
                temp_dict['wave_direction'] = np.radians(dir)
                temp_dict['convention'] = bemio_obj.body[i].bem_code
                temp_dict['influenced_dof'] = dofs
                
                if wavenumber or wavelength:
                    if temp_dict['water_depth'] == np.infty or omega**2*temp_dict['water_depth']/temp_dict['g'] > 20:
                        k = omega**2/temp_dict['g']
                    else:
                        k = newton(lambda x: x*np.tanh(x) - omega**2*temp_dict['water_depth']/temp_dict['g'], x0=1.0)/temp_dict['water_depth']
                    
                    if wavenumber:
                        temp_dict['wavenumber'] = k

                    if wavelength:
                        if k == 0.0:
                            temp_dict['wavelength'] = np.infty
                        else:
                            temp_dict['wavelength'] = 2*np.pi/k

                Fexc = np.empty(shape=bemio_obj.body[i].ex.re[:, dir_idx, omega_idx].shape, dtype=np.complex128)
                if from_wamit:
                    Fexc.real = bemio_obj.body[i].ex.re[:, dir_idx, omega_idx] * rho * g
                    Fexc.imag = bemio_obj.body[i].ex.im[:, dir_idx, omega_idx] * rho * g
                else:
                    Fexc.real = bemio_obj.body[i].ex.re[:, dir_idx, omega_idx]
                    Fexc.imag = bemio_obj.body[i].ex.im[:, dir_idx, omega_idx]
                temp_dict['diffraction_force'] = Fexc.flatten()
            
                try:
                    Fexc_fk = np.empty(shape=bemio_obj.body[i].ex.fk.re[:, dir_idx, omega_idx].shape, dtype=np.complex128)
                    if from_wamit:
                        Fexc_fk.real = bemio_obj.body[i].ex.fk.re[:, dir_idx, omega_idx] * rho * g
                        Fexc_fk.imag = bemio_obj.body[i].ex.fk.im[:, dir_idx, omega_idx] * rho * g
                    else:
                        Fexc_fk.real = bemio_obj.body[i].ex.fk.re[:, dir_idx, omega_idx]
                        Fexc_fk.imag = bemio_obj.body[i].ex.fk.im[:, dir_idx, omega_idx]
                    temp_dict['Froude_Krylov_force'] = Fexc_fk.flatten()

                except AttributeError:
                        # LOG.warning('\tNo Froude-Krylov forces found for ' + bemio_obj.body[i].name + ' at ' + str(dir) + \
                        #       ' degrees (omega = ' + str(omega) + '), replacing with zeros.')
                        temp_dict['Froude_Krylov_force'] = np.zeros((bemio_obj.body[i].ex.re[:, dir_idx, omega_idx].size,), dtype=np.complex128)

                difr_dict.append(temp_dict)

            # RadiationProblem + Hydrostatics variable equivalents
            for radiating_dof_idx, radiating_dof in enumerate(dofs):
                temp_dict = {}
                temp_dict['body_name'] = bemio_obj.body[i].name
                temp_dict['water_depth'] = bemio_obj.body[i].water_depth
                temp_dict['omega'] = omega
                temp_dict['rho'] = rho
                temp_dict['g'] = g
                temp_dict['influenced_dof'] = dofs
                temp_dict['radiating_dof'] = radiating_dof
                temp_dict['added_mass'] = bemio_obj.body[i].am.all[radiating_dof_idx, :, omega_idx].flatten()
                temp_dict['radiation_damping'] = bemio_obj.body[i].rd.all[radiating_dof_idx, :, omega_idx].flatten()

                if from_wamit:
                    temp_dict['added_mass'] = temp_dict['added_mass'] * rho
                    temp_dict['radiation_damping'] = temp_dict['radiation_damping'] * rho * omega

                if wavenumber or wavelength:
                    if temp_dict['water_depth'] == np.infty or omega**2*temp_dict['water_depth']/temp_dict['g'] > 20:
                        k = omega**2/temp_dict['g']
                    else:
                        k = newton(lambda x: x*np.tanh(x) - omega**2*temp_dict['water_depth']/temp_dict['g'], x0=1.0)/temp_dict['water_depth']
                    
                    if wavenumber:
                        temp_dict['wavenumber'] = k

                    if wavelength:
                        if k == 0.0:
                            temp_dict['wavelength'] = np.infty
                        else:
                            temp_dict['wavelength'] = 2*np.pi/k

                rad_dict.append(temp_dict)

    df = df.append(pd.DataFrame.from_dict(difr_dict).explode(['influenced_dof', 'diffraction_force', 'Froude_Krylov_force']))
    df = df.append(pd.DataFrame.from_dict(rad_dict).explode(['influenced_dof', 'added_mass', 'radiation_damping']))
    df = df.astype({'added_mass': np.float64, 'radiation_damping': np.float64, 'diffraction_force': np.complex128, 'Froude_Krylov_force': np.complex128})

    return df