/eng/ssb/auto/astroconda/bin/tests_run.sh --test-name test_miri_fixed \
    --test-context example \
    --test-root /data1/jwst_rt/ \
    --tests "general/miri/test_miri_fixed.py" \
    --test-files-only \
    --channel http://ssb.stsci.edu/conda-dev \
    --python 2.7 \
    jwst


