# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 14:46:43 2020

@author: Richard Reksoatmodjo

Parameter Quantity Dictionary
"""

EireneDict = {'PDENA':{'FileName':'EirAtom', 'Label': r'Neutral D Density $(m^{-3})$'},
              'PDENM':{'FileName':'EirMol', 'Label': r'Molecular D2 Density $(m^{-3})$'},
              'EDENA':{'FileName':'AtomEnergy', 'Label': r'Neutral D Energy Density $(eV*m^{-3})$'},
              'EDENM':{'FileName':'MolEnergy', 'Label': r'Molecular D2 Energy Density $(eV*m^{-3})$'}}


B2Dict = {'Ne': r'Electron Density $n_e\;(m^{-3})$',
          'Te': r'Electron Temperature $T_e\;(eV)$',
          'NI' : r'Ion (D+) Density $n_i\;(m^{-3})$',
          'Ti': r'Ion (D+) Temperature $T_i\;(eV)$',
          'DN': r'Particle Density Diffusivity $D\;(m^2/s)$',
          'KYE': r'Electron Thermal Diffusivity $\chi_e (m^2/s)$',
          'KYI': r'Ion Thermal Diffusivity $\chi_i (m^2/s)$',
          'NeuDen': r'Neutral Atom (D) Density $(m^{-3})$',
          'MolDen': r'Neutral Molecule (D2) Density $(m^{-3})$',
          'NeuTemp': r'Neutral Atom (D) Temperature (eV)',
          'MolTemp': r'Neutral Molecule (D2) Temperature (eV)',
          'IonFlx': r'Radial Particle Flux $s^{-1}$',
          'MolFlx': r'Radial Molecule Flux $m^{-2}s^{-1}$',
          'RadFlx': r'Radial Atomic Flux $m^{-2}s^{-1}$',
          'IonPol': r'Poloidal Particle Flux $m^{-2}s^{-1}$',
          'VLY': r'Radial Pinch $v_y (m^2/s)$',
          'SX': r'Poloidal Contact Area SX $(m^{2})$',
          'SY': r'Radial Contact Area SY $(m^{2})$',
          'SZ': r'Poloidal Cross-Sectional Area SZ $(m^{2})$',
          'VOL': r'Cell Volume VOL $(m^{2})$',
          'RadPinch': r'Anomalous Radial Pinch Velocity $v_y\;(m/s)$',
          'AHAL': r'Atomic $H_\alpha$ emissivity $(photons*m^{-3}s^{-1})$',
          'MHAL': r'Molecular $H_\alpha$ emissivity $(photons*m^{-3}s^{-1})$',
          'NTI' : r'Test Ion (D2+) Density $n_{ti}\;(m^{-3})$',
          'TestIonTemp' : r'Test Ion (D2+) Temperature $T_{testion}\;(eV)$'}