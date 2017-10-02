#! /usr/bin/env python

'''
import unittest2 as unittest
import numpy as np
import json
# from datetime import datetime
import time
import os
import sys
sys.path.insert(0, '/Users/chris/Workspace/nwb/api-python')
sys.path.insert(0, '/Users/chris/Workspace/datajoint-python')
os.chdir('/Users/chris/Workspace/nwb/ALM-1')
dj.config['database.host'] = 'localhost'
dj.config['database.user'] = 'chris'
dj.config['database.password'] = ''
dj.config['display.limit'] = 5
dj.config['safemode']= False
dj.config['ingest.database'] = 'tutorial_alm1_ingest'
dj.config['production.database'] = 'catalog_alm1_dimitri'
'''
import os

import re
import code
from decimal import Decimal
from datetime import datetime
from contextlib import contextmanager

import yaml

import datajoint as dj

from nwb import nwb_file
from nwb import nwb_utils

import h5py

from pymysql.err import IntegrityError


[re, code, contextmanager, yaml, dj, nwb_file, nwb_utils, h5py]
# 23456789_123456789_123456789_123456789_123456789_123456789_123456789_12345678


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
class Keyword(dj.Lookup):
    definition = """
    # Tag of study types
    keyword : varchar(24)
    """
    contents = zip(['behavior', 'extracellular', 'photostim'])

    def _make_tuples(self, key):
        TODO = True


@schema
class Study(dj.Manual):
    definition = """
    # Study
    study		: varchar(8)    # short name of the study
    ---
    study_description	: varchar(255)	#
    institution		: varchar(255)  # institution conducting the study
    lab			: varchar(255)	#  lab conducting the study
    reference_atlas	: varchar(255)	# e.g. "paxinos"
    """

    def _make_tuples(self, key):
        TODO = True


@schema
class StudyKeyword(dj.Manual):
    definition = """
    # Study keyword (see general/notes)
    -> Study
    -> Keyword
    """

    def _make_tuples(self, key):
        TODO = True


@schema
class Publication(dj.Manual):
    definition = """
    # Publication
    doi			: varchar(60)	# publication DOI
    ----
    full_citation	: varchar(4000)
    authors=''		: varchar(4000)
    title=''		: varchar(1024)
    """

    def _make_tuples(self, key):
        TODO = True


@schema
class RelatedPublication(dj.Manual):
    definition = """
    -> Study
    -> Publication
    """

    def _make_tuples(self, key):
        TODO = True


@schema
class Animal(dj.Manual):
    definition = '''
    animal_id		: int		# Janelia institution mouse IDs
    ---
    species		: varchar(30)
    date_of_birth	: date
    '''


@schema
class Surgery(dj.Manual):
    definition = """
    -> Animal
    ---
    surgery		: varchar(4000)	# description of surgery
    """

    def _make_tuples(self, key):
        TODO = True


@schema
class Virus(dj.Manual):

    TODO = True

    definition = """
    -> Animal
    """

    class InfectionSite(dj.Part):

        TODO = True

        definition = """
        -> Virus
        site		: tinyint	# virus infection site
        ---
        infection_x	: decimal(3,2)	# (mm)
        infection_y	: decimal(3,2)	# (mm)
        infection_z	: decimal(3,2)	# (mm)
        """


@schema
class BrainArea(dj.Lookup):

    definition = """
    brain_area		: varchar(12)
    """
    contents = zip(['M2', 'V1'])


@schema
class Session(dj.Imported):

    definition = '''
    -> Animal
    session		: int		# session
    ---
    session_date	: date		# session date
    session_suffix	: char(1)	# suffix for session disambiguation
    experimenter	: varchar(60)	# experimenter
    session_start_time	: datetime
    raw_data_path	: varchar(255)	# raw data file path
    recording_type	: varchar(8)	# e.g. actute
    -> InputFiles
    '''

    @property
    def key_source(self):
        return InputFiles()

    def _make_tuples(self, key):

        fname = key['nwb_file']

        use_nwb_file = False  # slower due to validation; more memory use
        if use_nwb_file:
            f = nwb_file.open(key['nwb_file'], None, 'r').file_pointer
        else:
            f = h5py.File(key['nwb_file'], 'r')

        # fname: /path/to/file/data_structure_ANM210861_20130701.nwb
        (__, __, aid, sdate) = os.path.split(fname)[-1:][0].split('_')

        # animal id / session id / session_date / session_suffix
        aid = aid[3:]  # ANMNNN -> NNN
        sdate = sdate.split('.')[:-1][0]  # drop '.nwb'
        if re.match('[a-zA-Z]', sdate):
            (sdate, sfx) = (sdate[:-1], sdate[-1:],)
        else:
            (sdate, sfx) = (sdate, '',)

        # Animal
        if not (Animal() & {'animal_id': aid}):
            Animal().insert1((aid, 'squeaky', '1970-01-01',))

        key['animal_id'] = aid
        key['session'] = int(sdate + '00')  # TODO: sfx -> NN{00:26}
        key['session_date'] = sdate
        key['session_suffix'] = sfx

        # various
        key['experimenter'] = f['general']['experimenter'][()].decode()
        key['raw_data_path'] = f['acquisition']['timeseries']['extracellular_traces']['ephys_raw_data'][()].decode()

        key['recording_type'] = 'TODO'  # TODO: recording_type

        key['session_start_time'] = f['session_start_time'][()].decode()
        key['session_start_time'] = '1970-01-01'  # TODO: strptime()

        f.close()
        self.insert1(key)


@schema
class Optogenetics(dj.Manual):

    TODO = True

    definition = """
    # Optogenetic stimulation information for the sesssion
    -> Session
    """

    class Site(dj.Part):

        TODO = True

        definition = """
        # Optogenetic stimulation site
        -> Optogenetics
        site_number		: tinyint	# optogenetic site number
        ---
        description		: varchar(255)	# optogenetic site description
        stimulation_method	: varchar(255)
        device			: varchar(60)
        location_x		: decimal(4,2)  # mm
        location_y		: decimal(4,2)  # mm
        excitation_lambda	: decimal(4,1)  # nm
        """


@schema
class TrialType(dj.Lookup):

    TODO = '''
    - split into expected / observed
    - implies updating Acquisition.TrialTypes to match
    '''

    definition = """
    # Trial Type  'response/correct response'
    trial_type		: varchar(16)
    """
    contents = zip([
        'HitR',
        'HitL',
        'ErrL',
        'ErrR',
        'NoLickR',
        'NoLickL',
        'LickEarly',
        'StimTrials'
    ])


@schema
class Ephys(dj.Computed):

    definition = """
    -> Session
    ---
    ad_bits		: tinyint	# A/D converter bits
    recording_marker	: varchar(30)	# e.g. "stereotaxic"
    ground_x		: decimal(4,2)	# (mm)
    ground_y		: decimal(4,2)	# (mm)
    ground_z		: decimal(4,2)	# (mm)
    """

    class Shank(dj.Part):

        definition = """
        -> Ephys
        shank		: tinyint	# shank of probe
        ---
        -> BrainArea
        posterior	: decimal(3,2)	# (mm)
        lateral		: decimal(3,2)	# (mm)
        """

    class Electrode(dj.Part):

        definition = """
        -> Ephys
        electrode	: tinyint	# electrode on probe
        ---
        -> Ephys.Shank
        electrode_x	: decimal(6,4)	# (mm) electrode map
        electrode_y	: decimal(6,4)	# (mm) electrode map
        electrode_z	: decimal(6,4)	# (mm) electrode map
        """

        # loc: ['M2'], P: 2.5, Lat: -1.2
        re = re.compile("loc: \['(.*?)'\], P: (.*), Lat: (.*)")

    def _make_tuples(self, key):

        key['nwb_file'] = (Session() & key).fetch1()['nwb_file']
        f = h5py.File(key['nwb_file'], 'r')

        g_gen = f['general']
        g_ephys = g_gen['extracellular_ephys']

        key['ad_bits'] = g_ephys['ADunit'][()].decode().split(' ')[0]  # 'NN '
        key['recording_marker'] = g_ephys['recording_marker'][()].decode()

        key['ground_x'] = Decimal(float(g_ephys['ground_coordinates'][0]))
        key['ground_y'] = Decimal(float(g_ephys['ground_coordinates'][1]))
        key['ground_z'] = Decimal(float(g_ephys['ground_coordinates'][2]))

        self.insert1(key, ignore_extra_fields=True)

        for shank_k in [k for k in g_ephys if 'shank_' in k]:

            shank = g_ephys[shank_k]
            shank_id = shank_k.split('_')[1]

            shank_loc = shank['location'][()].decode()

            m = Ephys.Electrode.re.match(shank_loc)
            (brain_area, shank_posterior, shank_lateral) = m.groups()

            key['shank'] = shank_id
            key['brain_area'] = brain_area
            key['posterior'] = Decimal(float(shank_posterior))
            key['lateral'] = Decimal(float(shank_lateral))

            self.Shank().insert1(key, ignore_extra_fields=True)

        electrode_group = g_ephys['electrode_group']
        electrode_map = g_ephys['electrode_map']

        for i in range(len(electrode_group)):

            key['electrode'] = i
            key['shank'] = electrode_group[i].decode()[-1:]
            key['electrode_x'] = electrode_map[i][0]
            key['electrode_y'] = electrode_map[i][1]
            key['electrode_z'] = electrode_map[i][2]

            self.Electrode().insert1(key, ignore_extra_fields=True)


@schema
class CellType(dj.Lookup):
    definition = """
    cell_type  : varchar(12)
    """
    contents = zip(['pyramidal', 'FS'])


@schema
class SpikeSorting(dj.Computed):

    definition = """
    -> Ephys
    ---
    identification_method	: varchar(60)
    """

    class Unit(dj.Part):

        definition = """
        -> SpikeSorting
        unit			: smallint	# unit number
        """

    class Type(dj.Part):

        definition = """
        -> SpikeSorting.Unit
        ---
        -> CellType
        """

    class Spikes(dj.Part):

        definition = """
        -> SpikeSorting.Unit
        ---
        spike_times		: longblob
        """

    class Waveform(dj.Part):

        definition = """
        -> SpikeSorting.Unit
        -> Ephys.Electrode
        ---
        waveform		: longblob	# uV
        """

    def _make_tuples(self, key):
        key['nwb_file'] = (Session() & key).fetch1()['nwb_file']
        f = h5py.File(key['nwb_file'], 'r')

        g_xlu = f['processing']['extracellular_units']

        key['identification_method'] = \
            g_xlu['identification_method'][()].decode()

        self.insert1(key, ignore_extra_fields=True)

        ure = re.compile('unit_[0-9]')
        for ukey in [k for k in g_xlu['UnitTimes'] if ure.match(k)]:

            # Unit
            unit = ukey.split('_')[1]  # 'unit_NN' -> NN
            key['unit'] = unit
            self.Unit().insert1(key, ignore_extra_fields=True)

            # Type
            cell_type = g_xlu['UnitTimes']['cell_types']

            (c_unit, c_str) = cell_type[int(unit)-1].decode().split(' - ')
            c_unit = c_unit.split('_')[1]

            try:
                assert(unit == c_unit)
            except:
                print('SpikeSorting.Unit mismatch:',
                      'unit:', unit, 'cell_unit:', c_unit)

            key['cell_type'] = c_str
            if c_str != '[]':
                self.Type().insert1(key, ignore_extra_fields=True)

            # Spikes
            key['spike_times'] = g_xlu['UnitTimes'][ukey]['times']
            self.Spikes().insert1(key, ignore_extra_fields=True)

            # Waveform
            g_wv = g_xlu['EventWaveform']
            key['waveform'] = g_wv[ukey]['data'][()]
            key['electrode'] = g_wv[ukey]['electrode_idx'][()][0]
            self.Waveform().insert1(key, ignore_extra_fields=True)


@schema
class Acquisition(dj.Computed):

    definition = """
    -> Session
    ---
    """

    class LickTrace(dj.Part):

        definition = """
        -> Acquisition
        ---
        lick_trace		: longblob
        timestamps		: longblob
        """

    class Trial(dj.Part):

        definition = """
        -> Acquisition
        trial			: smallint	# trial within a session
        ---
        start_time		: float
        stop_time		: float
        """

    class TrialTypes(dj.Part):

        definition = """
        -> Acquisition.Trial
        -> TrialType
        """

    class UnitInTrial(dj.Part):

        definition = """
        -> Acquisition.Trial
        -> SpikeSorting.Unit
        """

    class StimulusPresentation(dj.Part):

        definition = """
        -> Acquisition.Trial
        ---
        auditory_timestamp	: decimal(7,3)		# (s)
        auditory_cue		: tinyint		# -1 or 1
        pole_in_timestamp	: decimal(7,2)		# (s)
        pole_out_timestamp	: decimal(7,2)		# (s)
        """

    def _make_tuples(self, key):

        print('yoyo')
        key['nwb_file'] = (Session() & key).fetch1()['nwb_file']
        f = h5py.File(key['nwb_file'], 'r')

        g_acq = f['acquisition']
        self.insert1(key, ignore_extra_fields=True)

        # Acquisition.LickTrace
        key['lick_trace'] = g_acq['timeseries']['lick_trace']['data']
        key['timestamps'] = g_acq['timeseries']['lick_trace']['timestamps']

        # FIXME: ok? 21M rows.. hmm.
        # self.LickTrace().insert1(key, ignore_extra_fields=True)

        # Acquisition.Trial
        '''
        alternately:
        load TrialTypes, UnitInTrial, StimulusPresentation
        in single top-level trial iteration
        '''
        g_epochs = f['epochs']
        g_analysis = f['analysis']
        g_pres = f['stimulus']['presentation']

        for tkey in [k for k in g_epochs if 'trial_' in k]:

            # Acquisition.Trial
            tno = int(tkey.split('_')[1])  # TODO: dup var? (vs string)
            key['trial'] = tkey.split('_')[1]
            key['start_time'] = g_epochs[tkey]['start_time'][()]
            key['stop_time'] = g_epochs[tkey]['stop_time'][()]
            self.Trial().insert1(key, ignore_extra_fields=True)

            # Acquisition.UnitInTrial
            d_units = g_epochs[tkey]['units_present']
            if d_units[(0,)].decode() != 'NA':  # XXX: assuming convention
                for u in d_units:
                    key['unit'] = u.decode().split('_')[1]
                    self.UnitInTrial().insert1(key, ignore_extra_fields=True)

            # Acquisition.TrialTypes
            # XXX: ignoring trial_start_times yet doesn't match epoch start
            # trial_start_times *seems* like unrounded epoch[n]['start_time']?

            ttmat = g_analysis['trial_type_mat']
            ttstr = g_analysis['trial_type_string']
            for i in range(len(ttmat)):  # 8 trial types
                if ttmat[i][(tno - 1)]:  # off-by-1;
                    key['trial_type'] = ttstr[i].decode()
                    self.TrialTypes().insert1(key, ignore_extra_fields=True)

            # Acquisition.StimulusPresentation
            cuestamps = g_pres['auditory_cue']['timestamps']
            cuedat = g_pres['auditory_cue']['data']

            ipoletamps = g_pres['pole_in']['timestamps']
            opolestamps = g_pres['pole_out']['timestamps']

            tidx = tno - 1
            key['auditory_timestamp'] = Decimal(float(cuestamps[tidx]))
            key['auditory_cue'] = cuedat[tidx]
            key['pole_in_timestamp'] = ipoletamps[tidx]
            key['pole_out_timestamp'] = opolestamps[tidx]
            try:
                self.StimulusPresentation().insert1(
                    key, ignore_extra_fields=True)
            except IntegrityError:  # TODO: handle NaN in timestamps
                print('.StimulusPresentation error: session', key['session'],
                      'trial:', key['trial'])

        f.close()


if __name__ == '__main__':
    Session().populate()
    Ephys().populate()
    SpikeSorting().populate()
    Acquisition().populate()
    # code.interact(banner="alm-1 datajoint repl", local=locals())


'''
# misc repl stuffs:
fn = InputFiles().fetch()[0][0]
fh = h5py.File(fh, 'r')
g_gen = fh['general']
# misc hdf5 stuffs:
exp[exp.shape].decode()
'''

'''
SpikeSorting
XXX: how to query waveform by trial?
e.g.:
['UnitTimes'][unit]['times'][event (797)] -> timestamp
['UnitTimes'][unit]['trial_ids'][event (797)] -> trial
'''


#
# _make_tuple_wip
#
def _make_tuple_wip(self, key):

    print('yoyo')
    key['nwb_file'] = (Session() & key).fetch1()['nwb_file']
    f = h5py.File(key['nwb_file'], 'r')

    g_acq = f['acquisition']
    self.insert1(key, ignore_extra_fields=True)

    # Acquisition.LickTrace
    key['lick_trace'] = g_acq['timeseries']['lick_trace']['data']
    key['timestamps'] = g_acq['timeseries']['lick_trace']['timestamps']

    # FIXME: ok? 21M rows.. hmm.
    # self.LickTrace().insert1(key, ignore_extra_fields=True)

    # Acquisition.Trial
    '''
    alternately:
    load TrialTypes, UnitInTrial, StimulusPresentation
    in single top-level trial iteration
    '''
    g_epochs = f['epochs']
    g_analysis = f['analysis']
    g_pres = f['stimulus']['presentation']

    for tkey in [k for k in g_epochs if 'trial_' in k]:

        # Acquisition.Trial
        tno = int(tkey.split('_')[1])  # TODO: dup var? (vs string) - cleanup?
        key['trial'] = tkey.split('_')[1]
        key['start_time'] = g_epochs[tkey]['start_time'][()]
        key['stop_time'] = g_epochs[tkey]['stop_time'][()]
        self.Trial().insert1(key, ignore_extra_fields=True)

        # Acquisition.UnitInTrial
        d_units = g_epochs[tkey]['units_present']
        if d_units[(0,)].decode() != 'NA':  # XXX: assuming convention
            for u in d_units:
                key['unit'] = u.decode().split('_')[1]
                self.UnitInTrial().insert1(key, ignore_extra_fields=True)

        # Acquisition.TrialTypes
        # XXX: ignoring trial_start_times yet doesn't match epoch start
        # trial_start_times *seems* like unrounded epoch[n]['start_time']?

        ttmat = g_analysis['trial_type_mat']
        for i in range(len(ttmat)):  # 8 trial types
            if ttmat[i][(tno - 1)]:  # off-by-1;
                key['trial_type'] = g_analysis['trial_type_string'][i].decode()
                self.TrialTypes().insert1(key, ignore_extra_fields=True)

        # Acquisition.StimulusPresentation
        cuestamps = g_pres['auditory_cue']['timestamps']
        cuedat = g_pres['auditory_cue']['data']

        ipoletamps = g_pres['pole_in']['timestamps']
        opolestamps = g_pres['pole_out']['timestamps']

        tidx = tno - 1
        key['auditory_timestamp'] = Decimal(float(cuestamps[tidx]))
        key['auditory_cue'] = cuedat[tidx]
        key['pole_in_timestamp'] = ipoletamps[tidx]
        key['pole_out_timestamp'] = opolestamps[tidx]
        try:
            self.StimulusPresentation().insert1(key, ignore_extra_fields=True)
        except IntegrityError:  # TODO: handle NaN in timestamps
            print('.StimulusPresentation error: session', key['session'],
                  'trial:', key['trial'])

    f.close()


def _make_tuple_test(cls, func=None):
    if func is None:
        if '_orig_make_tuples' in cls.__dict__:
            print('reset')
            cls._make_tuples = cls._orig_make_tuples

    else:
        if '_orig_make_tuples' not in cls.__dict__:
            print('backup')
            cls._orig_make_tuples = cls._make_tuples

        print('replace')
        cls._make_tuples = func

    cls().delete()
    cls().populate()
    return cls()


def _make_tuple_loop():
    # doh: need 'force=true'alike for auto looping;
    # plus auto loop prevents redefining function in repl...
    pass

