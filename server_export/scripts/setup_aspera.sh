#!/usr/bin/env bash
# One-time setup: install aspera-cli (ascli) + the ascp transfer engine via
# a dedicated conda env, so it never touches the base env's torch install.
# No root needed -- bioconda's aspera-cli package + `ascli conf ascp install`
# both install into user-writable locations.
set -euo pipefail

ENV_NAME="aspera"

if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "=== creating conda env '${ENV_NAME}' ==="
    if ! conda create -y -n "${ENV_NAME}" -c bioconda -c conda-forge aspera-cli; then
        echo "conda create failed -- if the error mentions 'InvalidArchiveError' or"
        echo "'seeking backwards is not allowed', that's a conda tarball-extraction"
        echo "bug, not a real problem with the package. Fix: use mamba instead"
        echo "(handles archive extraction differently), e.g.:"
        echo "  conda install -n base -c conda-forge mamba -y"
        echo "  mamba create -y -n ${ENV_NAME} -c bioconda -c conda-forge aspera-cli"
        exit 1
    fi
else
    echo "conda env '${ENV_NAME}' already exists, skipping creation"
fi

echo "=== activating env and installing ascp transfer engine ==="
eval "$(conda shell.bash hook)"
conda activate "${ENV_NAME}"

ascli --version
ascli conf ascp install

echo
echo "=== done. verify ==="
ascli conf ascp info
