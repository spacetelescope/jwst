class DictWithAttributes(dict):
    """dict where valid keys are also attributes.

    Parameters
    ---------
    Basically all that dict takes.
    """
    def __init__(self, *args, **kwargs):
        super(DictWithAttributes, self).__init__(*args, **kwargs)
        self.__dict__ = self
