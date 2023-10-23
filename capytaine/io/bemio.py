import logging

import numpy as np
import pandas as pd
from scipy.optimize import newton

from capytaine.tools.optional_imports import import_optional_dependency

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
                temp_dict['period'] = 2*np.pi/omega
                temp_dict['rho'] = rho
                temp_dict['g'] = g
                temp_dict['wave_direction'] = np.radians(dir)
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

    df = pd.concat([
        pd.DataFrame.from_dict(difr_dict).explode(['influenced_dof', 'diffraction_force', 'Froude_Krylov_force']),
        pd.DataFrame.from_dict(rad_dict).explode(['influenced_dof', 'added_mass', 'radiation_damping'])
        ])
    df = df.astype({'added_mass': np.float64, 'radiation_damping': np.float64, 'diffraction_force': np.complex128, 'Froude_Krylov_force': np.complex128})

    return df

#####################
#  Export to Bemio  #
#####################

def save_as_bemio_hdf5(filename, dataset):
    h5py = import_optional_dependency("h5py")
    with h5py.File(filename, "w") as f:
        key = 0

        # name = f.create_dataset('body' + str(key+1) + '/properties/name', data=dataset.body_names.values[0])
        # name.attrs['description'] = 'Name of rigid body'

        forces = [dict(bemio_name='excitation', capytaine_name='excitation_force', plain_name='excitation force'),
                  dict(bemio_name='scattering', capytaine_name='diffraction_force', plain_name='scattering force'),
                  dict(bemio_name='froud_krylof', capytaine_name='Froude_Krylov_force', plain_name='Froude-Krylov force')]

        for force in forces:

            mag = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/{}/mag'.format(force['bemio_name']), data=abs(dataset[force['capytaine_name']].values))
            mag.attrs['units'] = ''
            mag.attrs['description'] = 'Magnitude of {}'.format(force['plain_name'])

            phase = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/{}/phase'.format(force['bemio_name']), data=np.angle(dataset[force['capytaine_name']].values))
            phase.attrs['units'] = 'rad'
            phase.attrs['description'] = 'Phase angle of {}'.format(force['plain_name'])

            re = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/{}/re'.format(force['bemio_name']), data=np.real(dataset[force['capytaine_name']].values))
            re.attrs['units'] = ''
            re.attrs['description'] = 'Real component of {}'.format(force['plain_name'])

            im = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/{}/im'.format(force['bemio_name']), data=np.imag(dataset[force['capytaine_name']].values))
            im.attrs['units'] = ''
            im.attrs['description'] = 'Imaginary component {}'.format(force['plain_name'])


        # # Write added mass information
        # am = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/added_mass/all',data=bemio_obj.body[key].am.all)
        # am.attrs['units for translational degrees of freedom'] = 'kg'
        # am.attrs['units for rotational degrees of freedom'] = 'kg-m^2'
        # am.attrs['description'] = 'Added mass. Frequency is the third dimension of the data structure.'
        #
        # rad = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/all', data=bemio_obj.body[key].rd.all)
        # rad.attrs['units'] = ''
        # rad.attrs['description'] = 'Radiation damping. Frequency is the third dimension of the data structure.'
        #
        # for m in range(bemio_obj.body[key].am.all.shape[0]):
        #
        #     for n in range(bemio_obj.body[key].am.all.shape[1]):
        #
        #         amComp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/added_mass/components/' + str(m+1) + '_' + str(n+1),data=np.array([bemio_obj.body[key].T, bemio_obj.body[key].am.all[m,n,:]]).transpose())
        #         amComp.attrs['units'] = ''
        #         amComp.attrs['description'] = 'Added mass components as a function of frequency'
        #
        #         radComp = f.create_dataset('body' + str(key+1) + '/hydro_coeffs/radiation_damping/components/' + str(m+1) + '_' + str(n+1),data=np.array([bemio_obj.body[key].T, bemio_obj.body[key].rd.all[m,n,:]]).transpose())
        #         radComp.attrs['units'] = ''
        #         radComp.attrs['description'] = 'Radiation damping components as a function of frequency'

    # Simulation parameters
    g = f.create_dataset('simulation_parameters/g', data=dataset.g)
    g.attrs['units'] = 'm/s^2'
    g.attrs['description'] = 'Gravitational acceleration'

    rho = f.create_dataset('simulation_parameters/rho', data=dataset.rho)
    rho.attrs['units'] = 'kg/m^3'
    rho.attrs['description'] = 'Water density'

    T = f.create_dataset('simulation_parameters/T', data=dataset.period)
    T.attrs['units'] = 's'
    T.attrs['description'] = 'Wave periods'

    w = f.create_dataset('simulation_parameters/w', data=dataset.omega)
    w.attrs['units'] = 'rad/s'
    w.attrs['description'] = 'Wave frequencies'

    water_depth = f.create_dataset('simulation_parameters/water_depth', data=dataset.water_depth)
    water_depth.attrs['units'] = 'm'
    water_depth.attrs['description'] = 'Water depth'

    wave_dir = f.create_dataset('simulation_parameters/wave_dir', data=dataset.wave_direction)
    wave_dir.attrs['units'] = 'rad'
    wave_dir.attrs['description'] = 'Wave direction'

    scaled = f.create_dataset('simulation_parameters/scaled', data=False)
    scaled.attrs['description'] = 'True: The data is scaled by rho*g, False: The data is not scaled by rho*g'

    code = f.create_dataset('bem_data/code', data="Capytaine")
    code.attrs['description'] = 'BEM code'
