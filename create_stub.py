from pynwb import NWBFile
from pynwb import NWBHDF5IO
from pynwb.ecephys import ElectricalSeries
from pynwb.file import Subject
from copy import deepcopy

def create_stub(data_path):
    out_path = data_path.parent/'stub'/data_path.name
    print(out_path)
    for nwb_file in data_path.glob('**/*.nwb'):
        save_path = out_path/nwb_file.relative_to(data_path)
        save_path.parent.mkdir(parents=True,exist_ok=True)
        print(f'saving {nwb_file} to {save_path}')
        with NWBHDF5IO(str(nwb_file), 'r') as io1:
            nwbfile = io1.read()
            with NWBHDF5IO(str(save_path), 'w') as io2:
                es_input = nwbfile.acquisition.get('ElectricalSeries', None)
                if es_input is not None:
                    stub_electrical_series = ElectricalSeries(name=es_input.name,
                                                              rate=es_input.rate,
                                                              conversion=es_input.conversion,
                                                              description=es_input.description,
                                                              starting_time=es_input.starting_time,
                                                              electrodes=deepcopy(es_input.electrodes),
                                                              data=es_input.data[:10, :])
                else:
                    stub_electrical_series = None

                if nwbfile.subject is not None:
                    sub = nwbfile.subject
                    sub_ = Subject(age=sub.age,description=sub.description,genotype=sub.genotype,sex=sub.sex,
                                   species=sub.species,subject_id=sub.subject_id,weight=sub.weight,
                                   date_of_birth=sub.date_of_birth)
                else:
                    sub_ = None

                nwbfile_stub = NWBFile(session_description=nwbfile.session_description,
                                       identifier=nwbfile.identifier,
                                       session_start_time=nwbfile.session_start_time,
                                       experimenter=nwbfile.experimenter,
                                       experiment_description=nwbfile.experiment_description,
                                       session_id=nwbfile.session_id,
                                       institution=nwbfile.institution,
                                       lab=nwbfile.lab,
                                       acquisition=[stub_electrical_series],
                                       trials=nwbfile.trials,
                                       subject=sub_)
                io2.write(nwbfile_stub)