from jwst.cube_build.cube_build_step import CubeBuildStep


def test_cube_pa_overrides_ifualign_coordinate_system():
    assert CubeBuildStep._resolve_coord_system("ifualign", 37.5) == "skyalign"


def test_ifualign_is_preserved_without_cube_pa():
    assert CubeBuildStep._resolve_coord_system("ifualign", None) == "ifualign"


def test_world_alias_is_normalized_to_skyalign():
    assert CubeBuildStep._resolve_coord_system("world", None) == "skyalign"
