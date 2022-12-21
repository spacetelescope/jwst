#!/usr/bin/env python

import importlib
import pkgutil

import pytest

import jwst


def dependencies(package, exclude: [str]):
    return [
        module[1]
        for module in pkgutil.walk_packages(
            package.__path__, prefix=package.__name__ + "."
        )
        if not any(exclude_module in module[1] for exclude_module in exclude)
    ]


MODULES = dependencies(jwst, exclude=["test", "time"])


@pytest.mark.parametrize(
    "module_name",
    MODULES,
)
def test_module_import(module_name):
    importlib.import_module(module_name)
