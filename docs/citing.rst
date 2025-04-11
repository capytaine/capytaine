======
Citing
======

If you are using Capytaine, please cite both following references:

Ancellin and Dias (2019), **Capytaine: a Python-based linear potential flow solver**, Journal of Open Source Software, 4(36), 1341, https://doi.org/10.21105/joss.01341

Babarit and Delhommeau (2015), **Theoretical and numerical aspects of the open source BEM solver NEMOH**, Proceedings of the 11th European Wave and Tidal Energy Conference (EWTEC2015),

or in bibtex form:

.. code-block:: bib

    @article{ancellin_capytaine_2019,
      author  = {Matthieu Ancellin and Fr{\'{e}}d{\'{e}}ric Dias},
      title   = {Capytaine: a {Python-based} linear potential flow solver},
      journal = {Journal of Open Source Software},
      year    = 2019,
      month   = {apr},
      volume  = {4},
      number  = {36},
      pages   = {1341},
      doi     = {10.21105/joss.01341},
      url     = {https://doi.org/10.21105%2Fjoss.01341},
    }

    @inproceedings{babarit_theoretical_2015,
      author    = {Babarit, Aur{\'e}lien and Delhommeau, G{\'e}rard},
      title     = {Theoretical and numerical aspects of the open source {BEM} solver {NEMOH}},
      year      = {2015},
      booktitle = {{Proceedings of the 11th European Wave and Tidal Energy Conference (EWTEC2015)}},
      address   = {Nantes, France}
    }


Forward speed
-------------

If you are using the forward speed feature of Capytaine, please cite the following paper:

Donatini et al., **Implementation of forward speed effects on an open source seakeeping solver**. In 6th MASHCON: International Conference on Ship Manoeuvring in Shallow and Confined Water, 2022

.. code-block:: bib

    @inproceedings{donatini2022implementation,
        title={Implementation of forward speed effects on an open source seakeeping solver},
        author={Donatini, Luca and Herdayanditya, Ivandito and Verao Fernandez, Gael and Pribadi, Ajie Brama Krishna and Lataire, Evert and Delefortrie, Guillaume},
        booktitle={6th MASHCON: International Conference on Ship Manoeuvring in Shallow and Confined Water},
        pages={20--33},
        year={2022},
        organization={Knowledge Centre for Manoeuvring in Shallow and Confined Water}
    }


Alternative Green functions
---------------------------

If you use alternative Green function implementations included in Capytaine, please cite the following papers:

HAMS' finite depth Green function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yingyi Liu, Shigeo Yoshida, Changhong Hu, Makoto Sueyoshi, Liang Sun, Junliang Gao, Peiwen Cong, Guanghua He. **A reliable open-source package for performance evaluation of floating renewable energy systems in coastal and offshore regions**. Energy Conversion and Management, 174 (2018): 516-536.

Yingyi Liu **HAMS: A Frequency-Domain Preprocessor for Wave-Structure Interactions—Theory, Development, and Application.** Journal of Marine Science and Engineering, (2019) 7: 81.

.. code-block:: bib

   @article{liu2018reliable,
      title={A reliable open-source package for performance evaluation of floating renewable energy systems in coastal and offshore regions},
      author={Liu, Yingyi and Yoshida, Shigeo and Hu, Changhong and Sueyoshi, Makoto and Sun, Liang and Gao, Junliang and Cong, Peiwen and He, Guanghua},
      journal={Energy conversion and Management},
      volume={174},
      pages={516--536},
      year={2018},
      publisher={Elsevier}
    }

   @article{liu2019hams,
      title={HAMS: A frequency-domain preprocessor for wave-structure interactions—Theory, development, and application},
      author={Liu, Yingyi},
      journal={Journal of Marine Science and Engineering},
      volume={7},
      number={3},
      pages={81},
      year={2019},
      publisher={MDPI}
    }


HAMS' infinite depth Green function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

H. Wu, C. Zhang, Y. Zhu, W. Li, D. Wan, F. Noblesse, **A global approximation to the Green function for diffraction radiation of water waves**, Eur. J. Mech. B Fluids 65 (2017) 54-64.

H. Liang, H. Wu, F. Noblesse, **Validation of a global approximation for wave diffraction-radiation in deep water**, Appl. Ocean Res. 74 (2018) 80-86.

.. code-block:: bib

   @article{wu2017global,
      title={A global approximation to the Green function for diffraction radiation of water waves},
      author={Wu, Huiyu and Zhang, Chenliang and Zhu, Yi and Li, Wei and Wan, Decheng and Noblesse, Francis},
      journal={European Journal of Mechanics-B/Fluids},
      volume={65},
      pages={54--64},
      year={2017},
      publisher={Elsevier}
    }

    @article{liang2018validation,
      title={Validation of a global approximation for wave diffraction-radiation in deep water},
      author={Liang, Hui and Wu, Huiyu and Noblesse, Francis},
      journal={Applied Ocean Research},
      volume={74},
      pages={80--86},
      year={2018},
      publisher={Elsevier}
    }


Source code
-----------

To cite the source code itself, you can use the DOI :code:`10.5281/zenodo.1426306`.
It represents all versions of Capytaine and always points to the most recent version.
If you want to cite a specific version, you can find its own DOI on `Zenodo <http://doi.org/10.5281/zenodo.1426306>`_.
