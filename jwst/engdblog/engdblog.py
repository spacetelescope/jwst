from ..lib.engdb_tools import ENGDB_Service
from ..stpipe import Step


class EngDBLogStep(Step):
    """Pipeline step to retrieve selected engineering mnemonic values

    Notes
    -----
    This is primarily developed to satisfy the following DMS requirements.
    However, it has been developed in hopes that it may provide useful
    diagnostics and/or as an example of step writing and engineering DB
    access.

    DMS-256 Calibration software shall access the DMS engineering database
        in order to derive the engineering information needed for science
        instrument calibration.

    DMS-458 The Calibration software shall retrieve relevant engineering
        data from the Engineering Database.
    """

    spec = """
    stime = string(default='2021-01-25')  # Start time
    etime = string(default='2021-01-27')  # End time
    verbosity = option('initial', 'all', default='initial')  # How much to report.
    engdb_url = string(default='http://localhost')  # Mock url
    """

    def process(self, mnemonics):
        """
        Step processing

        Parameters
        ----------
        mnemonics: str or [str (, ...)]
            The list of mnemonics to retrieve

        Returns
        -------
        {mnemonic: values[, ...]}
            `dict` of the specified mnemonics and their values
        """

        result = {}
        stime = self.stime
        etime = self.etime
        verbosity = self.verbosity
        edb = ENGDB_Service(base_url=self.engdb_url)

        if isinstance(mnemonics, str):
            mnemonics = [mnemonics]
        for mnemonic in mnemonics:
            try:
                values = edb.get_values(mnemonic, stime, etime)
            except Exception:
                self.log.info(
                    'Cannot retrieve info for {}'.format(
                        mnemonic
                    )
                )
                continue

            if len(values) < 1:
                self.log.info(
                    '{} has no entries in time range {}:{}'.format(
                        mnemonic, stime, etime
                    )
                )
                continue

            if verbosity == 'initial':
                result[mnemonic] = values[0]
                self.log.info(
                    '{}[{}:{}] = {}'.format(
                        mnemonic,
                        stime,
                        etime,
                        values[0]
                    )
                )
            elif verbosity == 'all':
                result[mnemonic] = values
                self.log.info(
                    '{}[{}:{}] = {}'.format(
                        mnemonic,
                        stime,
                        etime,
                        values
                    )
                )
                next

        return result
