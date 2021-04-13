import unittest
from os.path import abspath, dirname, join, isfile
import os
import numpy as np
# import pandas as pd
# import json
import matplotlib.pylab as plt
from scipy.interpolate import interp1d
import pygmsh
import numpy as np
import capytaine as cpt
import logging
import matplotlib.pyplot as plt
from capytaine.io.mesh_loaders import load_mesh


testdir = dirname(abspath(__file__))
datadir = join(testdir, 'data')

class TestPygmshIo(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.fbname = 'cyl'
        geom = pygmsh.opencascade.Geometry()
        geom.add_cylinder(x0=[0.0, 0.0, 0],
                         axis=[0.0, 0.0, -2],
                         radius=1,
                         angle=2 * np.pi,
                         char_length=1)
        mshRefFactor = 0.25
        mshArgs = ['-clscale', str(mshRefFactor),
           '-clcurv', str(360/50)]
        self.mesh = pygmsh.generate_mesh(geom,
                            dim=2,
                            extra_gmsh_arguments=mshArgs,
                            remove_lower_dim_cells=True,
                            geo_filename=self.fbname + '.geo',
                            msh_filename=self.fbname + '.stl',
                            mesh_file_type='stl')

    @classmethod
    def tearDownClass(self):
        pass # TODO

    def testReadStl(self):
        """
        Tests that Capytain can load STL files
        """
        cpt_mesh = load_mesh(filename=self.fbname + '.stl')
        fb = cpt.FloatingBody(mesh=cpt_mesh, name="my_body")
        expVal = 6.243988962054469
        self.assertAlmostEqual(expVal, cpt_mesh.volume)

    def testDirectLoadPygmsh(self):
        """
        Tests that mesh objects created by pygmsh can be loaded directly into
        Capytaine
        """

        def all_faces_as_tetra(cells):
            all_faces = []
            if 'tetra' in cells:
                all_faces.append(cells['tetra'])
            if 'triangle' in cells:
                triangles_as_tetra = np.empty((cells['triangle'].shape[0], 4), 
                    dtype=np.int)
                triangles_as_tetra[:, :3] = cells['triangle'][:, :]
                triangles_as_tetra[:, 3] = cells['triangle'][:, 2]  # Repeat one node to make a tetra
                all_faces.append(triangles_as_tetra)
            return np.concatenate(all_faces)

        cpt_mesh = cpt.Mesh(vertices=self.mesh.points, 
            faces=all_faces_as_tetra(self.mesh.cells_dict))
        fb = cpt.FloatingBody(mesh=cpt_mesh, name="my_body")
        expVal = 6.243988962054469
        self.assertAlmostEqual(expVal, cpt_mesh.volume)

    def test_all_io(self):
        # vertices, faces = load_VTP('data/SEAREV.vtp')
    
        for (loader, writer) in list(extension_dict.values()):
            try:
                writer('meshfile', vertices, faces)
                can_try_to_load = True
            except NotImplementedError:
                can_try_to_load = False
            
            if can_try_to_load:
                try:
                    loader('meshfile')
                except NotImplementedError:
                    pass
        
        os.remove('meshfile')


if __name__ == '__main__':
    unittest.main() 

# class TestResourceSpectrum(unittest.TestCase):

#     @classmethod
#     def setUpClass(self):
#         omega = np.arange(0.1,3.5,0.01)
#         self.f = omega/(2*np.pi)
#         self.Hs = 2.5
#         self.Tp = 8
#         df = self.f[1] - self.f[0]
#         Trep = 1/df
#         self.t = np.arange(0, Trep, 0.05)
            
#     @classmethod
#     def tearDownClass(self):
#         pass
    
#     def test_pierson_moskowitz_spectrum(self):
#         S = wave.resource.pierson_moskowitz_spectrum(self.f,self.Tp)
#         Tp0 = wave.resource.peak_period(S).iloc[0,0]
        
#         error = np.abs(self.Tp - Tp0)/self.Tp
        
#         self.assertLess(error, 0.01)
        
#     def test_bretschneider_spectrum(self):
#         S = wave.resource.bretschneider_spectrum(self.f,self.Tp,self.Hs)
#         Hm0 = wave.resource.significant_wave_height(S).iloc[0,0]
#         Tp0 = wave.resource.peak_period(S).iloc[0,0]
        
#         errorHm0 = np.abs(self.Tp - Tp0)/self.Tp
#         errorTp0 = np.abs(self.Hs - Hm0)/self.Hs
        
#         self.assertLess(errorHm0, 0.01)
#         self.assertLess(errorTp0, 0.01)

#     def test_surface_elevation_moments(self):
#         S = wave.resource.jonswap_spectrum(self.f, self.Tp, self.Hs)
#         eta = wave.resource.surface_elevation(S, self.t)
#         dt = self.t[1] - self.t[0]
#         Sn = wave.resource.elevation_spectrum(eta, 1/dt, len(eta.values), 
#                                               detrend=False, window='boxcar',
#                                               noverlap=0)

#         m0 = wave.resource.frequency_moment(S,0).m0.values[0]
#         m0n = wave.resource.frequency_moment(Sn,0).m0.values[0]
#         errorm0 = np.abs((m0 - m0n)/m0)

#         self.assertLess(errorm0, 0.01)

#         m1 = wave.resource.frequency_moment(S,1).m1.values[0]
#         m1n = wave.resource.frequency_moment(Sn,1).m1.values[0]
#         errorm1 = np.abs((m1 - m1n)/m1)

#         self.assertLess(errorm1, 0.01)

#     def test_surface_elevation_rmse(self):
#         S = wave.resource.jonswap_spectrum(self.f, self.Tp, self.Hs)
#         eta = wave.resource.surface_elevation(S, self.t)
#         dt = self.t[1] - self.t[0]
#         Sn = wave.resource.elevation_spectrum(eta, 1/dt, len(eta), 
#                                               detrend=False, window='boxcar',
#                                               noverlap=0)

#         fSn = interp1d(Sn.index.values, Sn.values, axis=0)
#         rmse = (S.values - fSn(S.index.values))**2
#         rmse_sum = (np.sum(rmse)/len(rmse))**0.5

#         self.assertLess(rmse_sum, 0.02)
    
#     def test_jonswap_spectrum(self):
#         S = wave.resource.jonswap_spectrum(self.f, self.Tp, self.Hs)
#         Hm0 = wave.resource.significant_wave_height(S).iloc[0,0]
#         Tp0 = wave.resource.peak_period(S).iloc[0,0]
        
#         errorHm0 = np.abs(self.Tp - Tp0)/self.Tp
#         errorTp0 = np.abs(self.Hs - Hm0)/self.Hs
        
#         self.assertLess(errorHm0, 0.01)
#         self.assertLess(errorTp0, 0.01)
    
#     def test_plot_spectrum(self):            
#         filename = abspath(join(testdir, 'wave_plot_spectrum.png'))
#         if isfile(filename):
#             os.remove(filename)
        
#         S = wave.resource.pierson_moskowitz_spectrum(self.f,self.Tp)
        
#         plt.figure()
#         wave.graphics.plot_spectrum(S)
#         plt.savefig(filename, format='png')
#         plt.close()
        
#         self.assertTrue(isfile(filename))

#     def test_plot_chakrabarti(self):            
#         filename = abspath(join(testdir, 'wave_plot_chakrabarti.png'))
#         if isfile(filename):
#             os.remove(filename)
        
#         D = 5
#         H = 10
#         L = 200

#         wave.graphics.plot_chakrabarti(H, L, D)
#         plt.savefig(filename)

#     def test_plot_chakrabarti_np(self):            
#         filename = abspath(join(testdir, 'wave_plot_chakrabarti_np.png'))
#         if isfile(filename):
#             os.remove(filename)
        
#         D = np.linspace(5, 15, 5)
#         H = 10 * np.ones_like(D)
#         L = 200 * np.ones_like(D)

#         wave.graphics.plot_chakrabarti(H, L, D)
#         plt.savefig(filename)
        
#         self.assertTrue(isfile(filename))

#     def test_plot_chakrabarti_pd(self):            
#         filename = abspath(join(testdir, 'wave_plot_chakrabarti_pd.png'))
#         if isfile(filename):
#             os.remove(filename)
        
#         D = np.linspace(5, 15, 5)
#         H = 10 * np.ones_like(D)
#         L = 200 * np.ones_like(D)
#         mvals = pd.DataFrame([H.flatten(),L.flatten(),D.flatten()],
#                          index=['H','L','D']).transpose()

#         wave.graphics.plot_chakrabarti(H, L, D)
#         plt.savefig(filename)
        
#         self.assertTrue(isfile(filename))
        

# class TestResourceMetrics(unittest.TestCase):

#     @classmethod
#     def setUpClass(self):
#         file_name = join(datadir, 'ValData1.json')
#         with open(file_name, "r") as read_file:
#             self.valdata1 = pd.DataFrame(json.load(read_file))
        
#         self.valdata2 = {}

#         file_name = join(datadir, 'ValData2_MC.json')
#         with open(file_name, "r") as read_file:
#             data = json.load(read_file)
#         self.valdata2['MC'] = data
#         for i in data.keys():
#             # Calculate elevation spectra
#             elevation = pd.DataFrame(data[i]['elevation'])
#             elevation.index = elevation.index.astype(float)
#             elevation.sort_index(inplace=True)
#             sample_rate = data[i]['sample_rate']
#             NFFT = data[i]['NFFT']
#             self.valdata2['MC'][i]['S'] = wave.resource.elevation_spectrum(elevation, 
#                          sample_rate, NFFT)

#         file_name = join(datadir, 'ValData2_AH.json')
#         with open(file_name, "r") as read_file:
#             data = json.load(read_file)
#         self.valdata2['AH'] = data
#         for i in data.keys():
#             # Calculate elevation spectra
#             elevation = pd.DataFrame(data[i]['elevation'])
#             elevation.index = elevation.index.astype(float)
#             elevation.sort_index(inplace=True)
#             sample_rate = data[i]['sample_rate']
#             NFFT = data[i]['NFFT']
#             self.valdata2['AH'][i]['S'] = wave.resource.elevation_spectrum(elevation, 
#                          sample_rate, NFFT)
        
#         file_name = join(datadir, 'ValData2_CDiP.json')
#         with open(file_name, "r") as read_file:
#             data = json.load(read_file)
#         self.valdata2['CDiP'] = data
#         for i in data.keys():
#             temp = pd.Series(data[i]['S']).to_frame('S')
#             temp.index = temp.index.astype(float)
#             self.valdata2['CDiP'][i]['S'] = temp
        

            
#     @classmethod
#     def tearDownClass(self):
#         pass

#     def test_kfromw(self):
#         for i in self.valdata1.columns:
#             f = np.array(self.valdata1[i]['w'])/(2*np.pi)
#             h = self.valdata1[i]['h']
#             rho = self.valdata1[i]['rho']
            
#             expected = self.valdata1[i]['k']
#             calculated = wave.resource.wave_number(f, h, rho).loc[:,'k'].values
#             error = ((expected-calculated)**2).sum() # SSE
            
#             self.assertLess(error, 1e-6)
    
#     def test_moments(self):
#         for f in self.valdata2.keys(): # for each file MC, AH, CDiP
#             datasets = self.valdata2[f]
#             for s in datasets.keys(): # for each set
#                 #print(f, s)
#                 data = datasets[s]
#                 for m in data['m'].keys():
#                     expected = data['m'][m]
#                     S = data['S']
                    
#                     if s == 'CDiP1' or s == 'CDiP6':
#                         f_bins=pd.Series(data['freqBinWidth'])
#                     else: 
#                         f_bins = None
                    
#                     calculated = wave.resource.frequency_moment(S, int(m),frequency_bins=f_bins).iloc[0,0]
#                     error = np.abs(expected-calculated)/expected
                    
#                     self.assertLess(error, 0.01) 
    

#     def test_metrics(self):
#        for f in self.valdata2.keys(): # for each file MC, AH, CDiP
#             datasets = self.valdata2[f]
            
#             for s in datasets.keys(): # for each set
                
                
#                 data = datasets[s]
#                 S = data['S']
#                 if f == 'CDiP':
#                     f_bins=pd.Series(data['freqBinWidth'])
#                 else: 
#                     f_bins = None
                
#                 # Hm0
#                 expected = data['metrics']['Hm0']
#                 calculated = wave.resource.significant_wave_height(S,frequency_bins=f_bins).iloc[0,0]
#                 error = np.abs(expected-calculated)/expected
#                 #print('Hm0', expected, calculated, error)
#                 self.assertLess(error, 0.01) 

#                 # Te
#                 expected = data['metrics']['Te']
#                 calculated = wave.resource.energy_period(S,frequency_bins=f_bins).iloc[0,0]
#                 error = np.abs(expected-calculated)/expected
#                 #print('Te', expected, calculated, error)
#                 self.assertLess(error, 0.01) 
                
#                 # T0
#                 expected = data['metrics']['T0']
#                 calculated = wave.resource.average_zero_crossing_period(S,frequency_bins=f_bins).iloc[0,0]
#                 error = np.abs(expected-calculated)/expected
#                 #print('T0', expected, calculated, error)
#                 self.assertLess(error, 0.01) 

#                 # Tc
#                 expected = data['metrics']['Tc']
#                 calculated = wave.resource.average_crest_period(S,frequency_bins=f_bins).iloc[0,0]**2 # Tc = Tavg**2
#                 error = np.abs(expected-calculated)/expected
#                 #print('Tc', expected, calculated, error)
#                 self.assertLess(error, 0.01) 

#                 # Tm
#                 expected = np.sqrt(data['metrics']['Tm'])
#                 calculated = wave.resource.average_wave_period(S,frequency_bins=f_bins).iloc[0,0]
#                 error = np.abs(expected-calculated)/expected
#                 #print('Tm', expected, calculated, error)
#                 self.assertLess(error, 0.01) 
                
#                 # Tp
#                 expected = data['metrics']['Tp']
#                 calculated = wave.resource.peak_period(S).iloc[0,0]
#                 error = np.abs(expected-calculated)/expected
#                 #print('Tp', expected, calculated, error)
#                 self.assertLess(error, 0.001) 
                
#                 # e
#                 expected = data['metrics']['e']
#                 calculated = wave.resource.spectral_bandwidth(S,frequency_bins=f_bins).iloc[0,0]
#                 error = np.abs(expected-calculated)/expected
#                 #print('e', expected, calculated, error)
#                 self.assertLess(error, 0.001) 

#                 # v
#                 if f == 'CDiP': # this should be updated to run on other datasets
#                     expected = data['metrics']['v']
#                     calculated = wave.resource.spectral_width(S,frequency_bins=f_bins).iloc[0,0]
#                     error = np.abs(expected-calculated)/expected

                       
#                     self.assertLess(error, 0.01) 

#                 if f == 'MC':
#                     expected = data['metrics']['v']
#                     calculated = wave.resource.spectral_width(S).iloc[0,0] # testing that default uniform frequency bin widths works 
#                     error = np.abs(expected-calculated)/expected

                       
#                     self.assertLess(error, 0.01)

                
#     def test_plot_elevation_timeseries(self):            
#         filename = abspath(join(testdir, 'wave_plot_elevation_timeseries.png'))
#         if isfile(filename):
#             os.remove(filename)
        
#         data = self.valdata2['MC']
#         temp = pd.DataFrame(data[list(data.keys())[0]]['elevation'])
#         temp.index = temp.index.astype(float)
#         temp.sort_index(inplace=True)
#         eta = temp.iloc[0:100,:]
        
#         plt.figure()
#         wave.graphics.plot_elevation_timeseries(eta)
#         plt.savefig(filename, format='png')
#         plt.close()
        
#         self.assertTrue(isfile(filename))
        
# class TestPerformance(unittest.TestCase):

#     @classmethod
#     def setUpClass(self):
#         np.random.seed(123)
#         Hm0 = np.random.rayleigh(4, 100000)
#         Te = np.random.normal(4.5, .8, 100000)
#         P = np.random.normal(200, 40, 100000)
#         J = np.random.normal(300, 10, 100000)
        
#         self.data = pd.DataFrame({'Hm0': Hm0, 'Te': Te, 'P': P,'J': J})
#         self.Hm0_bins = np.arange(0,19,0.5)
#         self.Te_bins = np.arange(0,9,1)

#     @classmethod
#     def tearDownClass(self):
#         pass
    
#     def test_capture_length(self):
#         L = wave.performance.capture_length(self.data['P'], self.data['J'])
#         L_stats = wave.performance.statistics(L)
        
#         self.assertAlmostEqual(L_stats['mean'], 0.6676, 3)
        
#     def test_capture_length_matrix(self):
#         L = wave.performance.capture_length(self.data['P'], self.data['J'])
#         LM = wave.performance.capture_length_matrix(self.data['Hm0'], self.data['Te'], 
#                         L, 'std', self.Hm0_bins, self.Te_bins)
        
#         self.assertEqual(LM.shape, (38,9))
#         self.assertEqual(LM.isna().sum().sum(), 131)
        
#     def test_wave_energy_flux_matrix(self):
#         JM = wave.performance.wave_energy_flux_matrix(self.data['Hm0'], self.data['Te'], 
#                         self.data['J'], 'mean', self.Hm0_bins, self.Te_bins)
        
#         self.assertEqual(JM.shape, (38,9))
#         self.assertEqual(JM.isna().sum().sum(), 131)
        
#     def test_power_matrix(self):
#         L = wave.performance.capture_length(self.data['P'], self.data['J'])
#         LM = wave.performance.capture_length_matrix(self.data['Hm0'], self.data['Te'], 
#                         L, 'mean', self.Hm0_bins, self.Te_bins)
#         JM = wave.performance.wave_energy_flux_matrix(self.data['Hm0'], self.data['Te'], 
#                         self.data['J'], 'mean', self.Hm0_bins, self.Te_bins)
#         PM = wave.performance.power_matrix(LM, JM)
        
#         self.assertEqual(PM.shape, (38,9))
#         self.assertEqual(PM.isna().sum().sum(), 131)
        
#     def test_mean_annual_energy_production(self):
#         L = wave.performance.capture_length(self.data['P'], self.data['J'])
#         maep = wave.performance.mean_annual_energy_production_timeseries(L, self.data['J'])

#         self.assertAlmostEqual(maep, 1754020.077, 2)
    
#     def test_ac_power_three_phase(self):
#         current = pd.DataFrame([[1,2,3],[4,5,6],[7,8,9],[10,11,12]], columns=['A1', 'A2', 'A3'])
#         voltage = pd.DataFrame([[1,5,9],[2,6,10],[3,7,11],[4,8,12]], columns=['V1', 'V2', 'V3'])
        
#         P1 = wave.performance.ac_power_three_phase(current, voltage, 1, False)
#         P1b = wave.performance.ac_power_three_phase(current, voltage, 0.5, False)
#         P2 = wave.performance.ac_power_three_phase(current, voltage, 1, True)
#         P2b = wave.performance.ac_power_three_phase(current, voltage, 0.5, True)
        
#         self.assertEqual(P1.sum()[0], 584)
#         self.assertEqual(P1b.sum()[0], 584/2)
#         self.assertAlmostEqual(P2.sum()[0], 1011.518, 2)
#         self.assertAlmostEqual(P2b.sum()[0], 1011.518/2, 2)
        
#     def test_dc_power(self):
#         current = pd.DataFrame([[1,2,3],[4,5,6],[7,8,9],[10,11,12]], columns=['A1', 'A2', 'A3'])
#         voltage = pd.DataFrame([[1,5,9],[2,6,10],[3,7,11],[4,8,12]], columns=['V1', 'V2', 'V3'])
        
#         P = wave.performance.dc_power(current, voltage)
    
#         self.assertEqual(P.sum()['Gross'], 584)
        
#     def test_plot_matrix(self):
#         filename = abspath(join(testdir, 'wave_plot_matrix.png'))
#         if isfile(filename):
#             os.remove(filename)
        
#         M = wave.performance.wave_energy_flux_matrix(self.data['Hm0'], self.data['Te'], 
#                         self.data['J'], 'mean', self.Hm0_bins, self.Te_bins)
        
#         plt.figure()
#         wave.graphics.plot_matrix(M)
#         plt.savefig(filename, format='png')
#         plt.close()
        
#         self.assertTrue(isfile(filename))
    
# class TestIO(unittest.TestCase):

#     @classmethod
#     def setUpClass(self):
#         self.expected_columns_metRT = ['WDIR', 'WSPD', 'GST', 'WVHT', 'DPD', 
#             'APD', 'MWD', 'PRES', 'ATMP', 'WTMP', 'DEWP', 'VIS', 'PTDY', 'TIDE']
#         self.expected_units_metRT = {'WDIR': 'degT', 'WSPD': 'm/s', 'GST': 'm/s', 
#             'WVHT': 'm', 'DPD': 'sec', 'APD': 'sec', 'MWD': 'degT', 'PRES': 'hPa', 
#             'ATMP': 'degC', 'WTMP': 'degC', 'DEWP': 'degC', 'VIS': 'nmi', 
#             'PTDY': 'hPa', 'TIDE': 'ft'}
        
#         self.expected_columns_metH = ['WDIR', 'WSPD', 'GST', 'WVHT', 'DPD', 
#             'APD', 'MWD', 'PRES', 'ATMP', 'WTMP', 'DEWP', 'VIS', 'TIDE']
#         self.expected_units_metH = {'WDIR': 'degT', 'WSPD': 'm/s', 'GST': 'm/s', 
#             'WVHT': 'm', 'DPD': 'sec', 'APD': 'sec', 'MWD': 'deg', 'PRES': 'hPa', 
#             'ATMP': 'degC', 'WTMP': 'degC', 'DEWP': 'degC', 'VIS': 'nmi', 
#             'TIDE': 'ft'}
        
#     @classmethod
#     def tearDownClass(self):
#         pass
    
#     ### Realtime data
#     def test_read_NDBC_realtime_met(self):
#         data, units = wave.io.read_NDBC_file(join(datadir, '46097.txt'))
#         expected_index0 = pd.datetime(2019,4,2,13,50)
#         self.assertSetEqual(set(data.columns), set(self.expected_columns_metRT))
#         self.assertEqual(data.index[0], expected_index0)
#         self.assertEqual(data.shape, (6490, 14))
#         self.assertEqual(units,self.expected_units_metRT)
            
#     ### Historical data
#     def test_read_NDBC_historical_met(self):
#         # QC'd monthly data, Aug 2019
#         data, units = wave.io.read_NDBC_file(join(datadir, '46097h201908qc.txt'))
#         expected_index0 = pd.datetime(2019,8,1,0,0)
#         self.assertSetEqual(set(data.columns), set(self.expected_columns_metH))
#         self.assertEqual(data.index[0], expected_index0)
#         self.assertEqual(data.shape, (4464, 13))
#         self.assertEqual(units,self.expected_units_metH)
        
#     ### Spectral data
#     def test_read_NDBC_spectral(self):
#         data, units = wave.io.read_NDBC_file(join(datadir, 'data.txt'))
#         self.assertEqual(data.shape, (743, 47))
#         self.assertEqual(units, None)

# if __name__ == '__main__':
#     unittest.main() 
