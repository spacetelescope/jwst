from .dictwithattrs import DictWithAttributes


class Member(DictWithAttributes):
    """Create a new member.

    Create a member based on a dict, remapping
    attributes as necessary.

    Parameters
    ---------
    initial_params : dict
        Initiial dict of key/values to populate the member.

    member_map : {input_param: member_param,}
        A dict where the keys are names of in input key names
        and the values are the names of the member key names.

    Notes
    -----
    If a map is given, and an input_param does not appear in the map, it
    is simply copyied using the same key.
    """
    def __init__(self, initial_params=None, member_map=None):
        if initial_params is None:
            initial_params = {}
        if member_map is None:
            super(Member, self).__init__(initial_params)
        else:
            super(Member, self).__init__()
            for key, value in initial_params.items():
                self.__setitem__(member_map.get(key, key), value)
