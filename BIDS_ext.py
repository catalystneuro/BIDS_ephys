import json
from pathlib import Path

import pandas as pd
from pynwb import NWBHDF5IO

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
        output_path = dataset_path/'BIDSExt'

    participants_df = pd.DataFrame(
        columns=['Species', 'ParticipantID', 'Sex', 'Birthdate', 'Age', 'Genotype', 'Weight'])
    conversion_dict = dict(milli=1e-3, micro=1e-6)
    dataset_desc_json = None

    dataset_path = Path(dataset_path)

    for nwb_file in dataset_path.glob('**/*.nwb'):
        channels_df = pd.DataFrame(columns=['channel_id', 'Contact_id', 'type', 'units', 'sampling_frequency'])
        sessions_df = pd.DataFrame(columns=['session_id', '#_trials', 'comment'])
        with NWBHDF5IO(str(nwb_file), 'r') as io:
            nwbfile = io.read()
            # subject info:
            if nwbfile.subject is not None:
                sb = nwbfile.subject
                participants_df.loc[len(participants_df.index)] = \
                    [sb.species, sb.subject_id, sb.sex, sb.date_of_birth, sb.age, sb.genotype, sb.weight]
                if sb.subject_id is not None:
                    subject_label = f'sub-{sb.subject_id}'
                else:
                    subject_label = f'sub-{sb.date_of_birth.strftime("%Y%m%dT%X")}'
            # dataset info:
            if dataset_desc_json is None:
                dataset_desc_json = dict(InstitutionName=nwbfile.institution, InstitutionalDepartmentName=nwbfile.lab,
                                         Name='Electrophysiology', BIDSVersion='1.0.X',
                                         Licence='CC BY 4.0',
                                         Authors=[list(nwbfile.experimenter) if nwbfile.experimenter is not None else None])
            # sessions info:
            trials_len = [len(nwbfile.trials) if nwbfile.trials is not None else None]
            sessions_df.loc[len(sessions_df.index)] = \
                [nwbfile.session_id, trials_len, nwbfile.session_description]
            session_label = f'ses-{nwbfile.session_id}'
            # channels_info:
            no_channels = nwbfile.acquisition['ElectricalSeries'].data.shape[1]
            sampling_frequency = nwbfile.acquisition['ElectricalSeries'].rate
            unit = nwbfile.acquisition['ElectricalSeries'].unit
            conversion_factor = [i if j == nwbfile.acquisition['ElectricalSeries'].conversion else ''
                                 for i, j in conversion_dict.items()][0]
            for chan_no in range(no_channels):
                channels_df.loc[len(channels_df.index)] = [chan_no, 'n.a.', 'neural signal', conversion_factor + unit,
                                                            sampling_frequency]

        # construct the folders:
        generic_ephys_name = f'{subject_label}_{session_label}_'
        subject_path = output_path/subject_label
        ses_path = subject_path/session_label
        data_path = ses_path/'ephys'
        data_path.mkdir(parents=True,exist_ok=True)
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
            channels_df.to_csv(bep_channels_path, sep='\t')
        # create sessions.json
        bep_sessions_path = subject_path/f'sub-{subject_label}_sessions.tsv'
        if not bep_sessions_path.exists():
            print(f'writing for subject: {subject_label}')
            sessions_df.to_csv(bep_sessions_path, sep='\t')

    # create participants.tsv:
    participants_df.to_csv(output_path/'participants.csv')
    # create dataset_desrciption.json
    with open(output_path/'dataset_description.json', 'w') as j:
        json.dump(dataset_desc_json, j)
