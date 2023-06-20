Software generating the Daejeon16 NN interaction described in:
A. M. Shirokov, I. J. Shin, Y. Kim, M. Sosonkina, P. Maris,
and J. P. Vary, Phys. Lett. B761 (2016) 87;
arXiv:1605.00413 [nucl-th] (2016).

___________________
IMPORTANT NOTICE:
Please note, there was a mistake in the previous version of the
fortran file, Daejeon16_public_v1.f90, which provided all
matrix elements multiplied by the factor of 25.
Use the current version, the file Daejeon16_public_v2.f90,
where this bug is corrected and which works much faster.
__________________


Authored by A. M. Shirokov, I. J. Shin, Y. Kim, M. Sosonkina,
P. Maris, and J. P. Vary with funding from the US
Department of Energy under Grants No. DESC0008485
(SciDAC/NUCLEI) and No. DE-FG02-87ER40371, from
the US National Science Foundation under Grant No.
1516096, from the Russian Foundation for Basic Research
Grant No. 15-02-06604-a, and the Ames Laboratory, operated
by Iowa State University under contract No.
DE-AC02-07CH11358. This work was also partially supported
by the Rare Isotope Science Project of Institute
for Basic Science funded by Ministry of Science, ICT and
Future Planning and National Research Foundation of
Korea (2013M7A1A1075764).

Copyright © 2016, A. M. Shirokov, I. J. Shin, Y. Kim, M. Sosonkina,
P. Maris, and J. P. Vary. All rights reserved.

If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from Iowa State University.
Additionally, redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
(1) Redistribution of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
(2) Redistribution in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with distribution.
(3) Neither the name of Iowa State University, the
U.S. Government, nor the names of its contributors may be used to
endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY IOWA STATE UNIVERSITY,
AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL IOWA STATE UNIVERSITY
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.


Description

The FORTRAN file Daejeon16_public_v2.f90 generates Daejeon16 NN interaction
matrix elements in the oscillator basis with hw=25 MeV. See the beginning
of the file Daejeon16_public_v2.f90 for explanations of how to use it.
