
                                                        FINGREEN3D

FINGREEN3D is an open-source package for computation of the free-surface Green's function under a finite water depth, which is the core-part of the boundary integral equations in the potential flow theory for analysis of wave-structure interactions. It is currently written in FORTRAN 90. The code is based on a state-of-the-art method developed in recent years by the authors. The main strategy of this method is to use a set of different series expansions in appropriate sub-regions. Therefore, we can avoid calculating the Green's function by direct time-consuming integrations.

FINGREEN3D is freely distributed under the LGPL License, Version 3.0, https://www.gnu.org/licenses/lgpl-3.0.html, and may be modified and
extended by researchers who intend to enhance its capabilities and port the code to other platforms.

It should be noted that, any modified version should be licensed under the LGPL License and be released open-publicly as well. The contributors can add their names in the "contributors list" ahead of the modified subroutine(s).

...

Note that the input & output parameters should keep the correct form as described in the accompanying ECM paper:

Yingyi Liu, Shigeo Yoshida, Changhong Hu, Makoto Sueyoshi, Liang Sun, Junliang Gao, Peiwen Cong, Guanghua He. A reliable open-source package for performance evaluation of floating renewable energy systems in coastal and offshore regions. Energy Conversion and Management, 174 (2018): 516-536.

Please cite the above paper in your relevant publications if the FINGREEN3D code or its executable (dynamic link library) has contributed to your work.

NO WARRANTY is given to FINGREEN3D.

The authors are not liable for any misuse or real and/or presumed damage which can arise from an improper use or removing/deleting some
parts of this software.

Adapted for Capytaine from https://github.com/YingyiLiu/FinGreen3D by Matthieu Ancellin, 2025.
