#!/usr/bin/env python
"""
Regression test okifying script.  Download the raw Jenkins output (right click
on "View as plain text" to save locally).  This script will attempt to parse
that output for failed fitsdiff comparisons, and prompt the user to okify or
skip each file.

Requires JFrog CLI (https://jfrog.com/getcli/) configured with credentials
that have write access to the jwst-pipeline repository.
"""

import sys
import re
import subprocess
import os
from argparse import ArgumentParser


RESULT_PATH_RE = re.compile(r"^E\s+(a:|---) (.*)$")
TRUTH_PATH_RE = re.compile(r"^E\s+(b:|\+\+\+) (.*)$")
DIFF_RE = re.compile(r"E\s+(fitsdiff|@@.*@@).*$")
DEPLOYING_RE = re.compile(r"^\[consumer_[0-9]+\] Deploying artifact:.*/jwst-pipeline-results/([^/]+)/.*$")
NO_DIFFERENCES_RE = re.compile(r"^E\s+No differences found.*$")


def parse_args():
    parser = ArgumentParser(description="Process regression test logs and prompt to okify results")
    parser.add_argument("log_path", help="path to regression test log file", metavar="log-path")
    parser.add_argument("--dry-run", action="store_true", help="pass the --dry-run flag to JFrog CLI")

    return parser.parse_args()


def artifactory_copy(source_path, target_path, dry_run=False):
    jfrog_args = ["--flat"]

    if dry_run:
        jfrog_args.append("--dry-run")

    subprocess.check_call(["jfrog", "rt", "cp"] + jfrog_args + [source_path, target_path])

    print("")


def get_artifactory_result_path(file_path, session_name):
    file_parts = file_path.split("/")
    assert len(file_parts) > 1

    artifactory_parts = ["jwst-pipeline-results", session_name]

    i = len(file_parts) - 2
    while i >= 0:
        if file_parts[i].startswith("test_"):
            artifactory_parts.extend(file_parts[i:])
            break
        i -= 1

    assert len(artifactory_parts) > 2, "Could not determine artifactory path of result file"

    return "/".join(artifactory_parts)


def get_artifactory_truth_path(file_path, input_path):
    file_parts = file_path.split("/")
    assert len(file_parts) > 1

    filename = file_parts[-1]

    i = len(file_parts) - 2
    subpath_parts = None
    while i >= 0:
        if file_parts[i].startswith("test_"):
            subpath_parts = file_parts[i+1:-1]
            break
        i -= 1

    if input_path[-1] == Ellipsis:
        artifactory_parts = input_path[:-1] + subpath_parts + [filename]
    else:
        artifactory_parts = input_path + [filename]

    return "/".join(artifactory_parts)


def main():
    args = parse_args()

    with open(args.log_path) as f:
        lines = f.readlines()

    session_name = None
    for line in reversed(lines):
        match = DEPLOYING_RE.match(line)
        if match:
            session_name = match.group(1)
            break

    assert session_name, "Could not determine test session name (e.g., 2019-07-08_jenkins-RT-JWST-314_0)"

    i = 0
    current_input_path = None
    while i < len(lines):
        if lines[i].startswith("input_path = ["):
            loc = {}
            exec(lines[i], globals(), loc)
            current_input_path = loc['input_path']
            i += 1
        elif DIFF_RE.match(lines[i]):
            block = [lines[i]]
            j = i + 1
            result_found = False
            truth_found = False
            while lines[j].startswith("E   ") and not DIFF_RE.match(lines[j]):
                if RESULT_PATH_RE.match(lines[j]):
                    if result_found:
                        # Sometimes if there is no diff, the fitsdiff output will
                        # omit the "fitsdiff: x.y.z" line that we're using to delineate
                        # file comparisons.  This just says: stop if we have already seen
                        # a result path and another one appears.
                        break
                    else:
                        result_found = True
                elif TRUTH_PATH_RE.match(lines[j]):
                    if truth_found:
                        break
                    else:
                        truth_found = True
                block.append(lines[j])
                j += 1

            no_differences = False
            for line in block:
                if NO_DIFFERENCES_RE.match(line):
                    no_differences = True

            if not no_differences:
                print("Diff output:\n")
                print("".join(block))

                result_match = RESULT_PATH_RE.match(block[1])
                truth_match = TRUTH_PATH_RE.match(block[2])

                if result_match and truth_match:
                    result_path = result_match.group(2).strip()
                    truth_path = truth_match.group(2).strip()
                else:
                    result_path = ""
                    truth_path = ""

                if not (result_path.startswith("/") and truth_path.startswith("/")):
                    k = i - 1
                    result_match = None
                    truth_match = None
                    while not (result_match and truth_match):
                        if not result_match:
                            result_match = RESULT_PATH_RE.match(lines[k])
                        if not truth_match:
                            truth_match = TRUTH_PATH_RE.match(lines[k])
                        k -= 1
                        assert i - k < 5, "Could not locate full result and truth file paths"

                    new_result_path = result_match.group(2).strip()
                    new_truth_path = truth_match.group(2).strip()

                    assert new_result_path.endswith(result_path)
                    assert new_truth_path.endswith(truth_path)

                    result_path = new_result_path
                    truth_path = new_truth_path

                artifactory_result_path = get_artifactory_result_path(result_path, session_name)
                artifactory_truth_path = get_artifactory_truth_path(truth_path, current_input_path)

                print(f"Artifactory result path: {artifactory_result_path}")
                print(f"Artifactory truth path: {artifactory_truth_path}\n")

                while True:
                    result = input("Enter 'o' to okify, 's' to skip: ")
                    if result not in ['o', 's']:
                        print(f"Unrecognized command '{result}', try again")
                    else:
                        break

                if result == 's':
                    print("Skipping\n")
                else:
                    print(f"Copying {artifactory_result_path} to {artifactory_truth_path}\n")
                    artifactory_copy(artifactory_result_path, artifactory_truth_path, dry_run=args.dry_run)

            i = j
        else:
            i += 1


if __name__ == "__main__":
    main()
