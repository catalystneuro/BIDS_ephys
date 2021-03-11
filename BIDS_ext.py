import json
from pathlib import Path

import pandas as pd
from pynwb import NWBHDF5IO
from pynwb.ecephys import ElectricalSeries

REQ_DATASETS = ['dataset_description.json', 'participants.tsv', 'sessions.tsv']


def bep_organize(dataset_path, output_path=None, move_nwb=False):
    """
    organize data according to teh BIDS extention proposal
    Parameters
    ----------
    dataset_path : [str, Path]
        path to the folder containing all the nwb datasets that need organization.
    """
    dataset_path = Path(dataset_path)
    if output_path is None:
        output_path = dataset_path.parent/'BIDSExt'/dataset_path.name

    participants_df = pd.DataFrame(
        columns=['Species', 'ParticipantID', 'Sex', 'Birthdate', 'Age', 'Genotype', 'Weight'])
    conversion_dict = dict(milli=1e-3, micro=1e-6)
    dataset_desc_json = None

    dataset_path = Path(dataset_path)
    file_count = 0
    sessions_count = 0
    for nwb_file in dataset_path.glob('**/*.nwb'):
        file_count += 1
        channels_df = pd.DataFrame(columns=['channel_id', 'Contact_id', 'type', 'units', 'sampling_frequency'])
        contacts_df = pd.DataFrame(columns=['x', 'y', 'z', 'impedance', 'contact_id', 'probe_id', 'Location'])
        probes_df = pd.DataFrame(columns=['probeID', 'type'])

        with NWBHDF5IO(str(nwb_file), 'r') as io:
            nwbfile = io.read()

            # subject info:
            if nwbfile.subject is not None:
                sb = nwbfile.subject
                if sb.subject_id is not None:
                    subject_label = f'sub-{sb.subject_id}'
                else:
                    subject_label = f'sub-{sb.date_of_birth.strftime("%Y%m%dT%H%M")}'
                if not participants_df['ParticipantID'].str.contains(subject_label).any():
                    participants_df.loc[len(participants_df.index)] = \
                        [sb.species, subject_label, sb.sex, sb.date_of_birth, sb.age, sb.genotype, sb.weight]
            else:
                subject_label = 'sub-none'
                if not participants_df['ParticipantID'].str.contains(subject_label).any():
                    participants_df.loc[len(participants_df.index)] = \
                        [None, subject_label, None, None, None, None, None]

            # dataset info:
            if dataset_desc_json is None:
                dataset_desc_json = dict(InstitutionName=nwbfile.institution, InstitutionalDepartmentName=nwbfile.lab,
                                         Name='Electrophysiology', BIDSVersion='1.0.X',
                                         Licence='CC BY 4.0',
                                         Authors=[
                                             list(nwbfile.experimenter) if nwbfile.experimenter is not None else None])
            # sessions info:
            subject_path = output_path/subject_label
            bep_sessions_path = subject_path/f'{subject_label}_sessions.tsv'
            if not bep_sessions_path.exists():
                print(f'writing for subject: {subject_label}')
                sessions_df = pd.DataFrame(columns=['session_id', '#_trials', 'comment'])
            else:
                sessions_df = pd.read_csv(bep_sessions_path,sep='\t')
            if nwbfile.session_id is not None:
                session_label = f'ses-{nwbfile.session_id}'
                # label_count = sessions_df['session_id'].str.contains(session_label).sum()
                # run_label = f'_run-{label_count}'
                # if label_count>0:
            else:
                session_label = f'ses-{nwbfile.session_start_time.strftime("%Y%m%dT%H%M")}'
            trials_len = len(nwbfile.trials) if nwbfile.trials is not None else None
            if not sessions_df['session_id'].str.contains(session_label).any():
                sessions_count += 1
                sessions_df.loc[len(sessions_df.index)] = \
                    [session_label, trials_len, nwbfile.session_description]

            # channels_info:
            es = [i for i in nwbfile.children if isinstance(i, ElectricalSeries)]
            if len(es) > 0:
                es = es[0]
                no_channels = es.data.shape[1]
                sampling_frequency = es.rate
                unit = es.unit
                conversion_factor = [i if j == es.conversion else ''
                                     for i, j in conversion_dict.items()][0]
                for chan_no in range(no_channels):
                    channels_df.loc[len(channels_df.index)] = [chan_no, chan_no, 'neural signal',
                                                               conversion_factor + unit,
                                                               sampling_frequency]

            #update ephys json:
            ephys_desc_json = dict(PowerLineFrequency=60.0)
            # contacts/probes info:
            e_table = nwbfile.electrodes
            if e_table is not None:
                for contact_no in range(len(e_table)):
                    contacts_df.loc[len(contacts_df.index)] = [e_table.x[contact_no],
                                                               e_table.y[contact_no],
                                                               e_table.z[contact_no],
                                                               e_table.imp[contact_no],
                                                               contact_no,
                                                               e_table.group_name[contact_no],
                                                               e_table.location[contact_no]]
            for probe_id in contacts_df['probe_id'].unique():
                probes_df.loc[len(probes_df.index)] = [probe_id, 'acute']#TODO: acute/chronic

        # construct the folders:
        generic_ephys_name = f'{subject_label}_{session_label}_'

        ses_path = subject_path/session_label
        data_path = ses_path/'ephys'
        data_path.mkdir(parents=True, exist_ok=True)

        # move nwbfile
        bep_nwbfile_path = data_path/(generic_ephys_name + 'ephys.nwb')
        if move_nwb:
            if not bep_nwbfile_path.exists():
                nwb_file.replace(bep_nwbfile_path)
        else:
            if not bep_nwbfile_path.exists():
                bep_nwbfile_path.symlink_to(nwb_file)

        # channels.tsv:
        bep_channels_path = data_path/(generic_ephys_name + 'channels.tsv')
        if not bep_channels_path.exists():
            channels_df.dropna(axis='columns', how='all', inplace=True)
            channels_df.to_csv(bep_channels_path, sep='\t', index=False)

        # probes/contacts:
        bep_probes_path = data_path/(generic_ephys_name + 'probes.tsv')
        if not bep_probes_path.exists():
            probes_df.dropna(axis='columns', how='all', inplace=True)
            probes_df.to_csv(bep_probes_path, sep='\t', index=False)
        bep_contacts_path = data_path/(generic_ephys_name + 'contacts.tsv')
        if not bep_contacts_path.exists():
            contacts_df.dropna(axis='columns', how='all', inplace=True)
            if len(contacts_df) > 0:
                contacts_df.to_csv(bep_contacts_path, sep='\t', index=False)

        # ephys description:
        bep_ephysdesc_path = data_path/(generic_ephys_name + 'ephys.json')
        with open(bep_ephysdesc_path, 'w') as j:
            json.dump(ephys_desc_json, j)
        # create sessions.json
        # sessions_df.dropna(axis='columns', how='all', inplace=True) # TODO, remove none fields
        sessions_df.to_csv(bep_sessions_path, sep='\t', index=False)

    # create participants.tsv:
    participants_df.dropna(axis='columns', how='all', inplace=True)
    participants_df.to_csv(output_path/'participants.tsv', sep='\t', index=False)

    # create dataset_desrciption.json
    with open(output_path/'dataset_description.json', 'w') as j:
        if all([True for au in dataset_desc_json['Authors'] if au is None]):
            _ = dataset_desc_json.pop('Authors')
        dataset_desc_tosave = {k: v for k, v in dataset_desc_json.items() if v is not None}
        json.dump(dataset_desc_tosave, j)
    print(f'total nwbfiles orgainzed {file_count}, sessions count {sessions_count}')

def bep_check(folder_path):
    pass
