#!/usr/bin/env bash
# One-time env setup on the HPC. Base conda already has torch 2.6+cu124 per
# prior verification (gpu02 = L20 48G). This only adds the WSI-specific bits.
set -euo pipefail
pip install --user openslide-python wsidicom timm huggingface_hub
echo "If openslide-python import fails with a libopenslide.so error, the C"
echo "library itself is missing (openslide-python only wraps it):"
echo "  conda install -c conda-forge openslide"
