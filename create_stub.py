from pathlib import Path

import numpy as np
from hdmf.common.table import DynamicTableRegion
from joblib import Parallel, delayed
from pynwb import NWBFile
from pynwb import NWBHDF5IO
from pynwb.base import TimeSeries
from pynwb.ecephys import ElectricalSeries
from pynwb.file import Subject

STUB_percentage = 0.01


def create_stub(data_path, n_jobs=10, **kwargs):
    data_path = Path(data_path)
    if data_path.suffix == '.nwb':
        out_path = data_path.with_name(f'stub_{data_path.name}')
        copy_nwb(data_path, out_path, **kwargs)
        print(f'saved {data_path} to {out_path}')
    elif data_path.is_dir():
        def convert_to_nwb(nwb_file, data_path):
            out_path = data_path.parent/'stub'/data_path.name
            save_path = out_path/nwb_file.relative_to(data_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            copy_nwb(nwb_file, save_path, **kwargs)
            print(f'saved {nwb_file} to {save_path}')

        Parallel(n_jobs=n_jobs)(delayed(convert_to_nwb)(nwb_file, data_path)
                                for nwb_file in data_path.glob('**/*.nwb'))
    else:
        raise ValueError('data_path should be path to .nwb or folder containing nwb')


def copy_nwb(nwb_file_path, nwb_save_path, **kwargs):
    with NWBHDF5IO(str(nwb_file_path), 'r') as io1:
        nwbfile = io1.read()
        with NWBHDF5IO(str(nwb_save_path), 'w') as io2:
            subject_info = extract_subject(nwbfile)

            nwbfile_stub = NWBFile(session_description=nwbfile.session_description,
                                   identifier=nwbfile.identifier,
                                   session_start_time=nwbfile.session_start_time,
                                   experimenter=nwbfile.experimenter,
                                   experiment_description=nwbfile.experiment_description,
                                   session_id=nwbfile.session_id,
                                   institution=nwbfile.institution,
                                   lab=nwbfile.lab,
                                   subject=subject_info)
            # Electrical series copy:
            es_name = kwargs.get('ElectricalSeries', 'ElectricalSeries')
            stub_perc = kwargs.get('stub', STUB_percentage)
            es = create_electricalseries(nwbfile, nwbfile_stub, series_name=es_name, stub=stub_perc)
            nwbfile_stub.add_acquisition(es)
            copy_trials(nwbfile, nwbfile_stub)
            # create processing:
            # TODO
            io2.write(nwbfile_stub)


def copy_electrodes_table(nwbfile_in, nwbfile_out):
    e_table = nwbfile_in.electrodes
    # create device:
    for device_name in nwbfile_in.devices:
        nwbfile_out.create_device(name=device_name)
    # create electroce groups:
    for e_group_name, e_group in nwbfile_in.electrode_groups.items():
        nwbfile_out.create_electrode_group(e_group_name,
                                           description=e_group.description,
                                           location=e_group.location,
                                           device=nwbfile_out.devices.get(e_group.device.name, None))
    if e_table is not None:
        default_electrode_colnames = ['x', 'y', 'z', 'group', 'group_name', 'imp', 'location', 'filtering',
                                      'id', 'rel_x', 'rel_y', 'rel_z', 'reference']
        for electrode_no in range(len(e_table)):
            in_dict = {}
            for colname in e_table.colnames:
                if colname in default_electrode_colnames:
                    if colname == 'group':
                        in_dict.update(
                            {colname: nwbfile_out.electrode_groups.get(e_table[colname].data[electrode_no].name, None)})
                    else:
                        in_dict.update({colname: e_table[colname].data[electrode_no]})
            nwbfile_out.add_electrode(**in_dict)
        for custom_e_column in set(e_table.colnames) - set(default_electrode_colnames):
            nwbfile_out.add_electrode_column(name=e_table[custom_e_column].name,
                                             description=e_table[custom_e_column].description,
                                             data=e_table[custom_e_column].data[()])
    return nwbfile_out.electrodes


def create_electricalseries(nwbfile_in, nwbfile_out, stub=STUB_percentage, series_name='ElectricalSeries', region=None):
    es_input = nwbfile_in.acquisition.get(series_name, None)

    if es_input is not None:
        electrodes_table = copy_electrodes_table(nwbfile_in, nwbfile_out)
        if region is None:
            region = list(range(es_input.data.shape[1]))
        electrodes_table_region = DynamicTableRegion(name=electrodes_table.name,
                                                     description=electrodes_table.description,
                                                     data=region,
                                                     table=electrodes_table)
        stub_length = np.round(es_input.data.shape[0]*stub).astype('int')
        stub_electrical_series = ElectricalSeries(name=es_input.name,
                                                  rate=es_input.rate,
                                                  conversion=es_input.conversion,
                                                  description=es_input.description,
                                                  starting_time=es_input.starting_time,
                                                  electrodes=electrodes_table_region,
                                                  data=es_input.data[:stub_length, :])
    else:
        stub_electrical_series = None
    return stub_electrical_series


def extract_subject(nwbfile):
    if nwbfile.subject is not None:
        sub = nwbfile.subject
        sub_ = Subject(age=sub.age, description=sub.description, genotype=sub.genotype, sex=sub.sex,
                       species=sub.species, subject_id=sub.subject_id, weight=sub.weight,
                       date_of_birth=sub.date_of_birth)
    else:
        sub_ = None
    return sub_


def copy_trials(nwbfile_in, nwbfile_out, stub=STUB_percentage):
    default_trial_columns = ['start_time', 'stop_time', 'tags', 'timeseries']
    trials_table = nwbfile_in.trials
    if trials_table is not None:
        for trial_no in range(len(trials_table)):
            in_dict = {}
            for colname in trials_table.colnames:
                if colname in default_trial_columns:
                    if 'timeseries' == colname:
                        ts_in = trials_table[colname][trial_no]
                        ts_kwargs = {i: j for i, j in ts_in.fields if i not in ['data', 'timestamps']}
                        stub_length = np.round(ts_in.data.shape[0]*stub).astype('int')
                        ts_kwargs.update(data=ts_in.data[: stub_length], timestamps=ts_in.timestamps[:stub_length])
                        in_dict.update(timeseries=TimeSeries(ts_kwargs))
                    else:
                        in_dict.update({colname: trials_table[colname][trial_no]})
            nwbfile_out.add_trial(**in_dict)

        for custom_e_column in set(trials_table.colnames) - set(default_trial_columns):
            nwbfile_out.add_trial_column(name=trials_table[custom_e_column].name,
                                         description=trials_table[custom_e_column].description,
                                         data=trials_table[custom_e_column].data[()])
