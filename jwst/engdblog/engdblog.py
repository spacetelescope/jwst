from jwst.lib.engdb_tools import ENGDB_Service
from jwst.lib.engdb_mast import EngdbMast
from jwst.stpipe import Step


class EngDBLogStep(Step):
    """
    Pipeline step to retrieve selected engineering mnemonic values.

    Notes
    -----
    This is primarily developed to satisfy the following DMS requirements.
    However, it has been developed in hopes that it may provide useful
    diagnostics and/or as an example of step writing and engineering DB
    access.

    * DMS-256 Calibration software shall access the DMS engineering database
      in order to derive the engineering information needed for science
      instrument calibration.
    * DMS-458 The Calibration software shall retrieve relevant engineering
      data from the Engineering Database.
    """

    spec = """
    stime = string(default='2022-01-25 02:00:00')  # Start time
    etime = string(default='2022-01-26 02:10:00')  # End time
    verbosity = option('initial', 'all', default='initial')  # How much to report.
    engdb_url = string(default=None)  # Mock url
    """  # noqa: E501

    def process(self, mnemonics):
        """
        Retrieve selected engineering mnemonic values.

        Parameters
        ----------
        mnemonics : str or list of str
            The list of mnemonics to retrieve.

        Returns
        -------
        result : dict
            The specified mnemonics and their values in the form of
            ``{mnemonic: values[, ...]}``.
        """
        result = {}
        stime = self.stime
        etime = self.etime
        verbosity = self.verbosity
        if self.engdb_url is not None:
            edb = ENGDB_Service(base_url=self.engdb_url)
        else:
            edb = EngdbMast()

        if isinstance(mnemonics, str):
            mnemonics = [mnemonics]
        for mnemonic in mnemonics:
            try:
                values = edb.get_values(mnemonic, stime, etime)
            except Exception:
                self.log.info("Cannot retrieve info for %s", mnemonic)
                continue

            if len(values) < 1:
                self.log.info("%s has no entries in time range %s:%s", mnemonic, stime, etime)
                continue

            if verbosity == "initial":
                result[mnemonic] = values[0]
                self.log.info("%s[%s:%s] = %s", mnemonic, stime, etime, str(values[0]))
            elif verbosity == "all":
                result[mnemonic] = values
                self.log.info("%s[%s:%s] = %s", mnemonic, stime, etime, str(values))

        return result
