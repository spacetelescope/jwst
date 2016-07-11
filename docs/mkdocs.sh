#!/bin/sh
SOURCE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RSYNC_OPT="-aiH --delete"
OUTDIR=$(pwd)/jwst_git
ERRORS=0

function usage
{
    printf "%s [-h|-o]\n\n\
    Not specifiying -o will place generated files into:\n    $OUTDIR
    \n
    -c  --clean  Runs 'make clean' on all projects\n\
    -o  --output {DIR} (default: \$PWD/jwst_git)\n\
    -h  --help   This message\n\n" "$(basename $0)"
}

function clean
{
    for doc in $SOURCE_DIR/*
    do
        [[ -f $doc ]] && continue
        printf "cleaning %-20s %s " "$(basename $doc)" "..."
        pushd $doc 2>&1 >/dev/null
            make clean &>/dev/null
            _ok $?
        popd 2>&1 >/dev/null
        echo
    done
}

function _ok
{
    if [[ $1 > 0 ]]; then
        /bin/echo -n "(failed: $1) "
        ERRORS=$((ERRORS + 1))
    fi
    return $1
}

function build_parent
{
    pushd "$SOURCE_DIR/.." &>/dev/null
    set -e
        python setup.py develop
    set +e
    popd &>/dev/null
    echo
    echo '----'
    echo
}

function install_travis_deps
{
    # because travis is a mess
    conda install -c $CONDA_CHANNELS stsci.sphinxext
}

while [[ $# -ge 1 ]]
do
    arg="$1"
    case $arg in
    -h|--help)
        usage
        exit 0
    ;;
    -c|--clean)
        clean
        exit 0
    ;;
    -o|--output)
        OUTDIR="$2"
        if [[ -z $2 ]]; then
            echo "No output directory specified..."
            exit 1
        fi

        # Replace with absolute path...
        mkdir -p $OUTDIR
        OUTDIR=$(cd $OUTDIR && pwd)
        shift
    ;;
    *)
    ;;
    esac
    shift
done

if [[ ! -e `which sphinx-build` ]]; then
    echo "Please install sphinx."
    exit 1
fi

python -c 'import stsci.sphinxext' &>/dev/null
_ok $? &>/dev/null

if [[ $? > 0 ]]; then
    echo "Please install stsci.sphinxext."
    exit 1
fi

echo "
!!
!! Sphinx 1.3.5 is required to build this documentation
!!             (downgrade if necessary)
!!
"

if [[ ! -e `which rsync` ]]; then
    echo "Please install rsync."
    exit 1
fi

LOGS=$OUTDIR/logs
FINAL=$OUTDIR/docs
TASK_CURRENT=1
TASK_TOTAL=0

mkdir -pv $OUTDIR
mkdir -pv $LOGS
mkdir -pv $FINAL

# A portable way of determining the top-level directory count
for d in $SOURCE_DIR/*
do
    if [[ -d $d ]]; then
        TASK_TOTAL=$((TASK_TOTAL + 1))
    fi
done

if [[ -n $TRAVIS ]]; then
    echo "Inside travis-ci environment, building in develop mode..."
    install_travis_deps
    build_parent
fi

for docs in $SOURCE_DIR/*
do
    [[ -f $docs ]] && continue

    name=$(basename $docs)
    logname="$LOGS/$name"

    printf "[%2d/%d] building %-20s %s " "$TASK_CURRENT" "$TASK_TOTAL" "$name" "..."
    pushd $docs 2>&1 >/dev/null
        /bin/echo -n "html "
        make html &>$logname-build-html.txt
        _ok $?
        if [[ $? == 0 ]]; then
            rsync $RSYNC_OPT $docs/build/html \
                $FINAL/$name/ &>$logname-rsync-html.txt
            _ok $?
        fi

        /bin/echo -n "latexpdf "
        yes 'q' | make latexpdf &>$logname-build-pdf.txt
        _ok $?
        if [[ $? == 0 ]]; then
            rsync $RSYNC_OPT $docs/build/latex/*.pdf \
                $FINAL/$name &>$logname-rsync-pdf.txt
            _ok $?
        fi
    popd 2>&1 >/dev/null
    echo
    TASK_CURRENT=$((TASK_CURRENT + 1))
done

if [[ $ERRORS > 0 ]]; then
    echo
    echo "Failures detected: $ERRORS"
    echo "Investigate '$LOGS' for more information..."
    echo
fi

exit $ERRORS
