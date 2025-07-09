def attrs_to_group_id(obs_meta):
    """
    Combine a number of file metadata values into a ``group_id`` string.

    Parameters
    ----------
    obs_meta : dict or ObjectNode
        A dictionary or ObjectNode containing meta.observation metadata, either
        model.meta.observation or meta["meta"]["observation"] from read_metadata.

    Returns
    -------
    str
        The group_id string.
    """
    obs_meta = dict(obs_meta.items())  # a bit circular to do to a dict, but needed for ObjectNode
    for key in [
        "program_number",
        "observation_number",
        "visit_number",
        "visit_group",
        "sequence_id",
        "activity_id",
        "exposure_number",
    ]:
        if key not in obs_meta:
            raise KeyError(f"Missing required keyword in meta.observation for group_id: {key}")
    return (
        f"jw{obs_meta['program_number']}"
        f"{obs_meta['observation_number']}"
        f"{obs_meta['visit_number']}"
        f"_{obs_meta['visit_group']}"
        f"{obs_meta['sequence_id']}"
        f"{obs_meta['activity_id']}"
        f"_{obs_meta['exposure_number']}"
    )
