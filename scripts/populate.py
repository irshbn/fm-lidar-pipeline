from pathlib import Path
import datetime as dt
import logging as log
from tqdm import tqdm
import numpy as np
from src import my_schema as schema

def do():
    # Logger config
    log.basicConfig(format='[%(asctime)s][%(levelname)s]: %(message)s',
                level=log.INFO)

    msk = dt.timezone(dt.timedelta(hours=+3))

    expdir = Path(__file__).resolve().parent.parent / "data" / "LIDAR_test_5mm_by_100um"
    log.info('Current experiment: %s', expdir.name.replace('_',' '))

    # Experiment data setup
    expdata = dict(exp_id=0,
                   exp_desc='Sensitivity test: 5 mm length, 100 um step',
                   exp_df=58.5e9,
                   exp_fmod=1e3)
    schema.Experiment.insert1(expdata, skip_duplicates=True)

    wfmdata_array=[]

    for filepath in tqdm(expdir.glob('*.csv'),
                         desc='Extracting data from files',
                         total=sum(1 for file in expdir.glob('*.csv'))):

        # Metadata extraction
        id_str, date, time = filepath.name.removeprefix('Wfm').removesuffix('.csv').split(sep='_')

        id = int(id_str.removesuffix('.'))

        dd, mm, yyyy = date.split(sep='.')
        h, m, s = time.split(sep='-')
        timestamp=dt.datetime(int(yyyy), int(mm), int(dd),
                              int(h), int(m), int(s),
                              tzinfo=msk)

        # Data extraction and packing

        blob = np.genfromtxt(filepath.resolve(),delimiter=',',missing_values='')
        delta = blob[1,1]
        data = blob[3:,1]

        wfmdata = dict(wfm_id=id, exp_id=0, wfm_timestamp=timestamp, wfm_delta=delta, wfm_data=data)
        wfmdata_array.append(wfmdata)

    # Pushing wfm data to DB
    log.info('Pushing data to database')
    schema.Waveform.insert(wfmdata_array, skip_duplicates=True)
    log.info('Preconditioning data')
    schema.WaveformPreconditioned.populate(display_progress=True)
    log.info('Performing least squares fit')
    schema.LeastSquaresFit.populate(display_progress=True)
    log.info('Population complete')

if __name__ == '__main__':
    do()