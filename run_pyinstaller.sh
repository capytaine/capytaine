#!/usr/bin/env sh

# Run in a environment created with
# conda create -n capy_install -c conda-forge capytaine pyinstaller

CAPY_BINARIES_DIR="$CONDA_PREFIX/lib/python3.10/site-packages/capytaine/green_functions/libs/"

BINARIES=""
for f in $CAPY_BINARIES_DIR/*; do
    BINARIES="--add-binary=$f:capytaine/green_functions/libs/ $BINARIES";
done

pyinstaller --clean --onefile --noconfirm $BINARIES examples/custom_dofs.py
