import datajoint as dj
from .preconditioning import precondition
from .processing import lsq_cosine_fit

schema = dj.Schema(f"{dj.config['database.user']}_schema")

@schema
class Experiment(dj.Manual):
    definition = """
        # LIDAR experiment metadata.
        exp_id      : int unsigned	# Unique expriment ID
        ---
        exp_desc=''	: varchar(255)  # Experiment description
        exp_df      : float			# Optical frequency range
        exp_fmod	: float			# Modulation frequency
        """


@schema
class Waveform(dj.Manual):
    definition = """
        # Waveform data blobs and metadata.
        wfm_id          : int           # Unqiue waveform ID
        ---
        -> Experiment
        wfm_timestamp   : timestamp     # acquisition datetime
        wfm_delta       : float         # wfm time step (sec)
        wfm_data        : longblob      # raw data vector
        """


@schema
class WaveformPreconditioned(dj.Computed):
    definition = """
        # Container for preconditioned waveform data blobs.
        -> Waveform
        ---
        pwfm_data	: longblob
        pwfm_delta	: float
        """

    def make(self, key):
        query = Waveform * Experiment & key
        data, delta, fmod = query.fetch1(
            'wfm_data', 'wfm_delta', 'exp_fmod')
        key['pwfm_data'] = precondition(
            data, delta, 20*fmod, 2/5, 121)
        key['pwfm_delta'] = delta
        self.insert1(key)


@schema
class LeastSquaresFit(dj.Computed):
    definition = """
    # Computed LIDAR beat frequency and fit quality parameters.
    -> WaveformPreconditioned
    ---
    lsq_frequency	: float	# Beat frequency (Hz)
    lsq_sigma		: float	# Frequency standard deviation (Hz)
    lsq_r2			: float	# Adjusted R-squared parameter 
    """

    def make(self, key):
        query = WaveformPreconditioned * Waveform & key
        pdata, dt = query.fetch1('pwfm_data', 'pwfm_delta')

        freq, sigma, r2 = lsq_cosine_fit(pdata)

        key['lsq_frequency'] = freq / dt
        key['lsq_sigma'] = sigma / dt
        key['lsq_r2'] = r2
        self.insert1(key)
