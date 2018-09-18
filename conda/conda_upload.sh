if [[ $1 == "develop" ]]; then
	anaconda -t $CONDA_UPLOAD_TOKEN upload -u mancellin -l nightly $HOME/miniconda/conda-bld/linux-64/capytaine*.tar.bz2 --force
else
	anaconda -t $CONDA_UPLOAD_TOKEN upload -u mancellin $HOME/miniconda/conda-bld/linux-64/capytaine*.tar.bz2 || true
fi
