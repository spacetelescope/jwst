"""Helpers for tests."""

from collections import namedtuple
from contextlib import contextmanager
import re
from tempfile import TemporaryDirectory
from pathlib import Path

from astropy.table import Table, vstack

from jwst.associations import AssociationRegistry, AssociationPool, generate
from jwst.associations.lib.counter import Counter
from jwst.associations.lib.utilities import is_iterable


class PoolParams(
    namedtuple(
        "PoolParams", ["path", "n_asns", "n_orphaned", "candidates", "valid_suffixes", "kwargs"]
    )
):
    """Class to hold parameters that define pool initial conditions."""

    __slots__ = ()

    def __new__(
        cls, path="", n_asns=0, n_orphaned=0, candidates=None, valid_suffixes=None, kwargs=None
    ):
        """Generate new PoolParams from parameters."""
        if not kwargs:
            kwargs = {}
        if candidates is None:
            candidates = []
            if valid_suffixes is None:
                valid_suffixes = ["cal", "calints", "cat"]
        return super(PoolParams, cls).__new__(
            cls, path, n_asns, n_orphaned, candidates, valid_suffixes, kwargs
        )


class BasePoolRule:
    """Basic Pool/Rule test class."""

    # Define the pools and testing parameters related to them.
    # Each entry is a tuple starting with the path of the pool.
    pools: list = []

    # Define the rules that SHOULD be present.
    # Each entry is the class name of the rule.
    valid_rules: list = []

    def test_rules_exist(self):
        """Test that level 3 registry rules are valid."""
        rules = registry_level3_only()
        assert len(rules) >= len(self.valid_rules)  # noqa: S101
        rule_names = get_rule_names(rules)
        for rule in self.valid_rules:
            assert rule in rule_names  # noqa: S101

    def test_run_generate(self):
        """Test association generation for list of pools."""
        rules = registry_level3_only()
        for ppars in self.pools:
            pool = combine_pools(ppars.path, **ppars.kwargs)
            asns = generate(pool, rules)
            assert len(asns) == ppars.n_asns, (  # noqa: S101
                ppars.path + f": n_asns not expected {len(asns)} {ppars.n_asns}"
            )
            for asn, candidates in zip(asns, ppars.candidates, strict=False):
                assert set(asn.candidates) == set(candidates)  # noqa: S101
            file_regex = re.compile(r".+_(?P<suffix>.+)\..+")
            for asn in asns:
                for product in asn["products"]:
                    for member in product["members"]:
                        if member["exptype"] == "science":
                            match = file_regex.match(member["expname"])
                            assert match is not None, (  # noqa: S101
                                ppars.path + ": No suffix match for {}".format(member["expname"])
                            )
                            assert match.groupdict()["suffix"] in ppars.valid_suffixes, (  # noqa: S101
                                ppars.path
                                + ": Suffix {} not valid".format(match.groupdict()["suffix"])
                            )


def make_megapool():
    """
    Combine the individual test pools into one.

    Notes
    -----
    This is meant to be run in the source tree in
    the folder the test pools reside in. The package
    does need to have been installed.
    `python -c 'import jwst.associations.tests.helpers as helpers; helpers.make_megapool()'`
    """
    pool_files = Path.glob("pool_*.csv")
    pool_files.sort()
    pool = combine_pools(pool_files)
    pool.write("mega_pool.csv", format="ascii", delimiter="|", overwrite=True)


# Basic utilities.
def t_path(partial_path):
    """
    Construct the full path for test files.

    Because the output of this utility is frequently passed
    as a command-line argument, a str object is required for
    playing nicely with argparse.

    Parameters
    ----------
    partial_path : Path or str
        The partial path for a test file.

    Returns
    -------
    str
        The string-converted full path.
    """
    test_dir = Path(__file__).parent

    return str(test_dir / partial_path)


def combine_pools(pools, **pool_kwargs):
    """
    Combine pools into a single pool.

    Parameters
    ----------
    pools: str, astropy.table.Table, [str|Table, ...]
        The pools to combine. Either a singleton is
        passed or and iterable can be passed.
        The entries themselves can be either a file path
        or an astropy.table.Table-like object.

    pool_kwargs: dict
        Other keyword arguments to pass to AssociationPool.read.

    Returns
    -------
    AssociationPool|astropy.table.Table
        The combined pool.
    """
    if not is_iterable(pools):
        pools = [pools]
    just_pools = []
    for pool in pools:
        if not isinstance(pool, Table):
            pool = AssociationPool.read(pool, **pool_kwargs)
        just_pools.append(pool)
    if len(just_pools) > 1:
        mega_pool = vstack(just_pools, metadata_conflicts="silent")
    else:
        mega_pool = just_pools[0].copy(copy_data=True)

    # Replace OBS_NUM and ASN_CANDIDATE_ID with actual numbers, if
    # necessary
    expnum = Counter(start=0)
    obsnum = Counter(start=0)
    acid = Counter(start=999)
    local_env = locals()
    global_env = globals()
    for row in mega_pool:
        mega_pool[row.index] = [
            parse_value(v, global_env=global_env, local_env=local_env) for v in row
        ]

    return mega_pool


def parse_value(v, global_env=None, local_env=None):
    """Evaluate if indicated."""
    if global_env is None:
        global_env = globals()
    if local_env is None:
        local_env = locals()

    result = v
    try:
        m = re.match("@!(.+)", v)
    except TypeError:
        pass
    else:
        if m:
            result = eval(m.group(1), global_env, local_env)  # noqa: S307
    return result


def fmt_cand(candidate_list):
    """
    Format the candidate field.

    Parameters
    ----------
    candidate_list: iterator
        An iterator with each element a 2-tuple of:
            cid: int
                Candidate ID
            ctype: str
                Candidate type

    Returns
    -------
    candidate_list_field: str
        A string of the list of candidates, with any evaluation
        performed.
    """
    evaled_list = []
    for cid, ctype in candidate_list:
        if isinstance(cid, int):
            if ctype == "observation" and cid < 1000:
                cid_format = "o{:0>3d}"
            elif cid >= 1000 and cid < 3000:
                cid_format = "c{:0>4d}"
            else:
                cid_format = "r{:0>4d}"
        else:
            cid_format = cid

        cid_str = cid_format.format(cid)
        evaled_list.append((cid_str, ctype))

    return str(evaled_list)


def fmt_fname(expnum):
    """Format the filename."""
    return f"jw_{expnum:0>5d}_uncal.fits"


def generate_params(request):
    """Generate param for pytest.fixtures."""
    return request.param


@contextmanager
def mkstemp_pool_file(pools, **pool_kwargs):
    """Make an actual pool file."""
    pool = combine_pools(pools, **pool_kwargs)
    with TemporaryDirectory() as path:
        pool_path = path + "pool"
        pool.write(
            pool_path,
            format="ascii",
            delimiter="|",
        )
        yield pool_path


def generate_pool_paths(request):
    """Fixture to create temporary files for pools."""
    pool_file = t_path(request.param)
    with mkstemp_pool_file(pool_file) as pool_path:
        yield pool_path


def get_rule_names(rules):
    """Return rule names found in a registry.

    Parameters
    ----------
    rules : AssociationRegistry
        The registry to look through.

    Returns
    -------
    rule_names : list
        The list of rule names.
    """
    return [rule.rule_name() for rule_name, rule in rules.items()]


def level3_rule_path():
    """Return the path to the level 3 rules."""
    return t_path("../lib/rules_level3.py")


def level2_rule_path():
    """Return the path to the level 2 rules."""
    return t_path("../lib/rules_level2b.py")


def registry_level3_only(global_constraints=None):
    """Get registry with only Level3 rules."""
    return AssociationRegistry(
        definition_files=[level3_rule_path()],
        include_default=False,
        global_constraints=global_constraints,
    )


def registry_level2_only(global_constraints=None):
    """Get registry with only Level2 rules."""
    return AssociationRegistry(
        definition_files=[level2_rule_path()],
        include_default=False,
        global_constraints=global_constraints,
    )
