;; test with guix shell -f capytaine.scm python -- python3 -c 'import capytaine; print(capytaine.__version__)'

(define-module (gnu packages python-capytaine)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix utils)
  #:use-module (guix packages)
  #:use-module (guix git-download)
  #:use-module (guix build-system python)
  #:use-module (guix build-system pyproject)
  #:use-module (gnu packages python)
  #:use-module (gnu packages build-tools)
  #:use-module (gnu packages check)
  #:use-module (gnu packages commencement)
  #:use-module (gnu packages python-xyz)
  #:use-module (gnu packages python-science)
  )

(package
  (name "python-capytaine")
  (version "2.0")
  (source (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/capytaine/capytaine")
             (commit "v2.0")))
       (file-name (git-file-name name version))
       (sha256
        (base32 "0gcv6l771xz5rksl14bv6h48ba0002rrcsmvbc5d3xxr5vlzw3aa"))
       ))
  (build-system pyproject-build-system)
  (native-inputs (list python-toolchain gfortran-toolchain meson-python python-pytest))
  (propagated-inputs (list python-numpy python-scipy python-pandas python-xarray))
  (home-page "https://github.com/capytaine/capytaine")
  (synopsis "Python BEM solver for linear potential flow, based on Nemoh")
  (description "Python BEM solver for linear potential flow, based on Nemoh")
  (license license:gpl3))

(package
  (name "python-capytaine")
  (version "2.1.dev")
  (source (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/capytaine/capytaine")
             (commit "c69dac76cb786e9bdf945402d93895cbf163a290")))
       (file-name (git-file-name name version))
       (sha256
        (base32 "1ypsbkb3z4kakgwvk46263hgr80m0i76js5vxfjkqzlzkm9hhcjw"))
       ))
  (build-system pyproject-build-system)
  (native-inputs (list python-toolchain gfortran-toolchain meson-python python-pytest))
  (propagated-inputs (list python-numpy python-scipy python-pandas python-xarray python-rich))
  (home-page "https://github.com/capytaine/capytaine")
  (synopsis "Python BEM solver for linear potential flow, based on Nemoh")
  (description "Python BEM solver for linear potential flow, based on Nemoh")
  (license license:gpl3))
