#! /usr/bin/env python

'''
import unittest2 as unittest
import numpy as np
import json
from datetime import datetime
import os
import sys
sys.path.insert(0, '/Users/chris/Workspace/nwb/api-python')
sys.path.insert(0, '/Users/chris/Workspace/datajoint-python')
os.chdir('/Users/chris/Workspace/nwb/ALM-1')
dj.config['database.host'] = 'localhost'
dj.config['database.user'] = 'chris'
dj.config['database.password'] = ''
dj.config['display.limit'] = 5
dj.config['ingest.database'] = 'tutorial_alm1_ingest'
dj.config['production.database'] = 'catalog_alm1_dimitri'
'''
import os

import re
import code
from contextlib import contextmanager

import yaml

import datajoint as dj

from nwb import nwb_file
from nwb import nwb_utils

import h5py


[re, code, contextmanager, yaml, dj, nwb_file, nwb_utils, h5py, schema]

# 23456789_123456789_123456789_123456789_123456789_123456789_123456789_12345678

# OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD #
# OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD #
# OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD #

@contextmanager
def kv(dct, key):
    yield dct[key]


def test_open(infile):
    read_via_nwb = True
    if read_via_nwb:
        return nwb_file.open(infile, None, 'r').file_pointer
    else:
        return h5py.File(infile, 'r')


def dump_experiment_data(nwb):

    with kv(nwb, 'general') as gen:

        exp = gen['experimenter']
        print('experimenter: ', exp[exp.shape].decode())

        with kv(gen, 'extracellular_ephys') as eep:
            print('extracellular_ephys:', eep)

            # em.shape -> (32,3); eg.shape -> (32,); therefore combine
            with kv(eep, 'electrode_map') as em:
                for (m, g) in zip(em, eep['electrode_group']):
                    [m[0], m[1], m[2], g.decode()]  # 0.0, 0.0, 0.0, 'shankN'

            for shank in [k for k in eep if 'shank' in k]:
                print('shank: ', shank)
                # eep[shank][{'location','device','description'}][()].decode()

    with kv(nwb, 'processing') as proc:

        print('processing: ', proc)

        with kv(proc, 'extracellular_units') as ecu:

            print('extracellular_units: ', ecu)
            '''
            unittimes/unit_NN[30] == eventwaveform/unit_NN[30]
            unittimes/unit_NN:
            - source
            - unit_description
            - times[eventwaveform/unit_NN/num_samples]
            - trial_ids[eventwaveform/unit_NN/num_samples] -> epoch/trial_NNN
            eventwaveform/unit_NN:
            - electrode_idx
            - num_samples
            - sample_length
            - data[num_samples, sample_length]
            - timestamps[num_samples]
            # correlate w/r/t epochs?
            '''
            # code.interact(banner='ecu', local=locals())
            u_re = re.compile('^unit_[0-9]+')  # filter out unit_list
            with kv(ecu, 'UnitTimes') as utimes:
                ewf = ecu['EventWaveform']
                for unit_k in [u for u in utimes if u_re.match(u)]:
                    t_unit = utimes[unit_k]
                    e_unit = ewf[unit_k]
                    print('unit_k: ', unit_k)
                    # print('t_unit: ', t_unit)
                    # print('e_unit: ', e_unit)

    # needs correlation with stimulus/epochs
    with kv(nwb, 'epochs') as ep:
        # extract 'trial'
        # note: 'trial' is not NWB specified member
        [ep]

    with kv(nwb, 'stimulus') as stim:

        print('stimulus', stim)

        with kv(stim, 'presentation') as pres:

            print('presentation', pres)

            with kv(pres, 'auditory_cue') as cue:
                print('auditory_cue', cue)

            with kv(pres, 'pole_in') as pin:
                print('pole_in', pin)

            with kv(pres, 'pole_out') as pout:
                print('pole_out', pout)


if __name__ == '__old_main__':
    infile = ('/Users/chris/Workspace'
              '/crcns/download/tools/download'
              '/alm-1/datafiles/nwb_files/'
              'data_structure_ANM210861_20130701.nwb')
    fh = test_open(infile)
    dump_experiment_data(fh)
    code.interact(banner="nwb repl w/'fh' var as an hdf5 filehandle",
                  local=locals())

# OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD #
# OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD #
# OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD # OLD #

# NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW #
# NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW #
# NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW # NEW #


def bogus():
    ''' bogus function to make my navigation shortcuts work '''
    pass


schema = dj.schema(dj.config['ingest.database'], locals())
schema.drop(force=True)
schema = dj.schema(dj.config['ingest.database'], locals())

# 'data_structure_ANM210861_20130701.nwb'
nwbfiledir = 'data'


@schema
class InputFiles(dj.Lookup):
    definition = '''
    nwb_file: varchar(255)
    '''

    contents = [[os.path.join(nwbfiledir, f)]
                for f in os.listdir(nwbfiledir) if f.endswith('.nwb')]


@schema
class Animal(dj.Manual):
    definition = '''
    animal_id		: int		# Janelia institution mouse IDs
    ---
    species		: varchar(30)
    date_of_birth	: date
    '''


@schema
class Session(dj.Imported):
    definition = '''
    -> Animal
    session		: int
    ---
    -> InputFiles
    '''

    @property
    def key_source(self):
        return InputFiles()

    def _make_tuples(self, key):
        print(yaml.dump(key))


if __name__ == '__main__':
    Session().poplulate()
    code.interact(banner="alm-1 datajoint repl", local=locals())
