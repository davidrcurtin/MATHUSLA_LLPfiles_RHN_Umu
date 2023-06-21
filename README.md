# MATHUSLA_LLPfiles_RHN
LLP production and decay files for simulation of heavy neutral lepton LLP production (Ue, Umu or Utau dominance) at the HL-LHC and decay in the MATHUSLA (or other transverse LLP) detector 


By David Curtin and Jaipratap Grewal, Jun 2023

BENCHMARK MODEL:

heavy neutral leptons (HNL) (also called right-handed neutrinos RHN) added to SM. Our main reference for model, production, and decay is 1805.08567. 

In generality, there are three HNL N1, N2, N3 which each can mix with active neutrinos v_e, v_mu, v_tau neutrinos, via a 3x3 mixing matrix U_ij. We consider three simplified scenarios corresponding to the PBC benchmark points BC6,7,8 in 1901.09966, where we consider one RHN LLP, N, which has mixing Ue with v_e only, or mixing U_mu with v_mu only, or mixing U_tau with v_tau only. In the below, UX refers to U_e, U_mu or U_tau depending on which benchmark model is being used. The files for the three benchmark models are identically formatted, so usage is the same. 

Each of these three simplified models has two free parameters:
- LLP mass mN
- mixing angle with corresponding active neutrino U, written as U^2 = Usq.

We concentrate on mN < 5 GeV. 

The most important production modes above 0.1 GeV are exotic D and B meson decays, but DY production (pp>W/Z>N+X) and tau-lepton decays to N are also important. 

RHN LLPs inherit all SM couplings of the active neutrino it mixes with, suppressed by mixing U. Therefore, N decays via weak decays N -> l+ W* or v Z*, with W*/Z* > SM fermions. 


AVAILABLE LLP MASSES:
==================

The easiest way to get a list of the available mN values in each benchmark model is to look at the numerical suffixes of the RHN_Umu_LLPweight4vectorWZlist_mN_YYY files. These production modes are available for all considered mN, and hence these files are always present for all considered mN. 




LLP PRODUCTION:
=================

Four production modes for the RHN N are separately considered: B meson decay, D meson decay, tau decay, and DY production via W*/Z*. B/D-decay is the most important, but the EW production modes are important to extend sensitivity for higher masses/shorter lifetimes. 


(1) B/D decay to RHN

FONLL was used to calculate double-differential production cross sections ds/dpTdy for B and D mesons at the HL-LHC, which are used to simulate unweighted B-meson momentum vectors from the IP, see Appendix A1. The BtoNbr_UX.pdf and DtoNBr_UX.pdf figures show the branching ratio of B and D mesons to N calculated following 1805.08567. Each exclusive decay was simulated using full 2- and 3-body decay kinematics to turn the simulated B/D-meson 4-vectors into LLP 4-vectors for each different LLP mass mN. 

Each of the files 
RHN_UX_LLPweight4vectorDmesonlist_mN_YYY.csv, 
RHN_UX_LLPweight4vectorBmesonlist_mN_YYY.csv 
in the All_RHN_UX_LLPweight4vectors_from_BDWZtau subdirectory contain a set of WEIGHTED LLP 4-vectors for that mass mN = YYYY (in GeV). (This is because each B/D 4-vector is decayed separately depending on whether it could represent a B0, B+, Bs, Bc, D0, D+, which all have  different decays to N.)
Note that the phi-angles around the beam axis of these LLPs are not physically distributed, since it is assumed that LLP 4-vectors will be rotated randomly around the beam axis prior to being shot towards MATHUSLA (see instructions below). These files have format (weight, E, pX, pY, pZ) with Z = beam axis. 

For each B/D decay event, the weight represents (0.5) * (the number of LLPs that a single 4-vector represents for Usq = 1), so simply multiply that weight by (2 * Usq) to get the real number of LLP Production events represented by the sample. 

WARNING 1: note the required factor of 2. This is because the factor of 2 needed to obtain the number of produced B-mesons from the FONLL cross section was not included in these event weights (sorry). This factor of 2 is NOT required in the weights for electorweak production events. 

WARNING 2: Note that this weight is set assuming you use all 4-vectors in the file. If you use some number N_MC events < the total number N_total in the file, the weights should be multiplied by 
(2 * Usq * N_total/N_MC) 
so they still represent the actual number of LLPs produced at the LHC. 
(This is OK since the different weights just correspond to different B/D-species, and hence there are only a limited number of different 'kinds' of weights that keep repeating in the file. If the different weights were representing kinematic properties, we should not just take the N_total/N_MC ratio but instead the ratio of total vs used weights.)


NB: If mN is too large to permit N production in D decay, then there will be no corresponding decay file. For all considered mN, B decay to N is always open, except for a few higher mass points in the Utau scenario. In this case, those B decay files are blank with a size of 1 byte, rather than missing.


(2) Tau decay

Madgraph was used to simulate tau production at the HL-LHC, see Appendix A2. 

The TautoNbr_UX.pdf figure shows the branching ratio(s) of tau to N calculated following 1805.08567. Each exclusive decay was simulated using full 2- and 3-body decay kinematics (but not angular dependencies in the matrix element) to turn the simulated tau-lepton 4-vectors into LLP 4-vectors for each different LLP mass mN. 

Each of the files 
RHN_UX_LLPweight4vectorTaulist_mN_YYYY.csv
in the All_RHN_UX_LLPweight4vectors_from_BDWZtau subdirectory contain a set of weighted LLP 4-vectors for that mass mN = YYYY (in GeV) IN THE SAME FORMAT AS THE LLP 4-VECTORS FROM B/D DECAY. USE IN ALMOST THE SAME WAY, EXCEPT FOR A DIFFERENT FACTOR OF 2 IN THE WEIGHT:
If you use N_MC of the sample of size N_total, each event represents
weight * Usq * N_total/N_MC 
LLP production events. (The correction factor needed for the FONLL cross sections is not needed here). 

NB: If mN is too large to permit N production in tau decay, then there will be no corresponding decay file. 


(2) DY production

Madgraph was used to simulate N production via s-channel Z/W the HL-LHC, see Appendix A2. 
Each of the files 
RHN_UX_LLPweight4vectorWZlist_mN_YYY.csv 
contains a set of (identically weighted) LLP 4-vectors for that mass mN = YYY (in GeV), in the IDENTICAL format to the LLP 4-vectors from tau decay. 




LLP DECAY LENGTH:
=============

Computed following 1805.08567, see RHNctauUeUmuUtau.pdf. 

For different mN in GeV, the decay length ctau in meters for Usq = 1 is supplied as a space-separated tablein

RHNctauUX.dat

Simply devide the second column by Usq to get the decay length for your mixing angle. 



LLP DECAY:
=========== 

Somewhat tricky for mN ~ GeV due to hadronization uncertainties, though the effect on a DV analysis is likely not as significant as for the long-lived scalar in the SM+S model, due to the absence of hierarchical Yukawa-ordered couplings to SM quarks. 

For a displaced vertex analysis, it is important to approximately capture the higher multiplicity of charged tracks resulting from LLP decays above ~ GeV, compared to exlusive hadron final state calculations that are typically limited to a very small multiplicity. To accurately treat the low-mass regime while also allowing for hadronization to produce multi-track vertices at higher masses, we adopt the following strategy:

- For low masses below the multi-hadron threshold (0.42 GeV for Ue and Utau dominated, 0.53 GeV for Umu dominated) we compute exclusive decay branching fractions following 1805.08567. See RHNBrUX_exclusive.pdf for the exclusive branching fractions. N decay events are generated using 2- and 3-body decay kinematics in the N restframe. 

- For higher masses above the multi-hadron threshold, Madgraph5 and the model SM_HeavyN_CKM_AllMasses_LO are used to decay N's via the respective mixing, then hadronize in Pythia8, extract detector-stable N decay products and save in N restframe. Care is taken to define the parton-level decay processes to only allow hadronic final states that are kinematically accessible for a given N mass mN. (Above ~ GeV, the partonic decay branching ratios simulated by the madgraph model, shown in in RHNBrUX_partonic.pdf which was also computed following 1805.08567, should be reasonably reliable.) More details on how these decays were generated is below in A.3. 


Unweighted decay events are supplied in the LLP restframe in the RHN_UX_hadronic_decays_geant subdirectory. There is one file for each LLP mass mN indicated in the file name suffix, with the file name text indicating the method of generating these decay events. Each file contains all decay modes in the appropriate proportions. In each file, different decay events are separated by two blank lines. 
The first row of each event block is the original LLP 4-vector used in the simulation. This can be ignored (since the simulation used for LLP decays does not match the kinematic distribution of actual LLP production at the HL-LHC.)
Subsequent rows contain lists of LLP decay final states IN THE LLP RESTFRAME, in format 
E, pX, pY, pZ, mass, PID, particle name. Note not all PIDs have a 'particle name' entry, that column is just for convenience. 



Note: The precise N branching fractions in the ~ 0.5 -  2 GeV mN window are subject to significant uncertainties, but the actual effect on DV reconstruction efficiencies should not be major, you can estimate it by comparing efficiency at 0.5 vs 2 GeV. Since reconstruction efficiency should increase roughly monotonically with mass, the difference between the start and end of the hadronic uncertainty window should give an idea of the relative uncertainty within the window. 





USAGE FOR MATHUSLA GEANT SIMULATIONS
=======================================

Select some LLP mass mN to study, say mN = 3.08517 GeV. 

You can consider each production mode separately, or lump them all together (which you can easily do if you keep track of the individual event weights). For each production mode, you will need the corresponding LLP 4-vector file, say 
RHN_Ue_LLPweight4vectorBmesonlist_mN_3.08517.csv.

For your chosen mass, you will also need the N decay file
vN_Ntoall_inclDs_3.08517.txt

as well as the lifetime from RHNctauUX.dat. 


1. Define rapidity and azimuthal angle ranges ((eta_min, eta_max), (phi_min, phi_max)) that roughly bracket the solid angle acceptance of MATHUSLA. (Overestimate so that no part of MATHUSLA falls outside of that range.)

2. Select a subset of LLP 4-vectors (preferably sequential subset of file, but for sure avoid duplicates) called {MCsample} from RHN_Ue_LLPweight4vectorBmesonlist_mN_3.08517.csv that you want to use in your simulations, or use all of them. Let the number of 4-vectors in this sample be N_MC. 

Each LLP 4-vector i in this sample {MCsample} has weight
w_i = (2 * Usq * N_total/N_MC) 
where N_total is the total number of events in the file. 
This weight corresponds to the number of LLPs that this simulated 4-vector represents. Note you do not need to specify Usq right now, since it won't affect kinematics. Later, it will set overall production rate and LLP lifetime. 

NB: for Tau/DY production modes, OMIT the factor of 2. 


3. From {MCsample}, select the subset {MCsample_etacut} that has rapidity in the range (eta_min, eta_max).

Furthermore, due to the symmetry around the beam axis, the phi-angle can be freely chosen for each event. For each surviving LLP 4-vector in {MCsample_etacut}, choose a random phi-angle in the range (phi_min, phi_max) and rotate the 4-vector around the beam axis to have that angle. The chance that a random LLP has a phi angle in the MATHUSLA acceptance is 
f_phi = (phi_max - phi_min)/2pi

Give each rotated LLP 4-vector i in our sample {MCsample_etacut} a weight
w'_i = w_i * f_phi. 
This weight corresponds to the number of LLPs flying through (or close to) MATHUSLA that this simulated 4-vector represents. 


4. Now we have to choose Usq. This will numerically set all the weights wi and set the lifetime in accordance with the lifetime curve in RHNctauUX.dat.

For each LLP 4-vector i in {MCsample_etacut}, find lengths along its trajectory (L1, L2) where it enters and exits the MATHUSLA detector volume that we include in our analysis. The chance of the LLP to decay in the detector is then

P_decay_i  = Exp(-L1/(b*ctau)) - Exp(-L2/(b*ctau))

where b = sqrt(px^2+py^2+pz^2)/mX is the boost of the LLP. 

Next, choose an  LLP decay position along the LLP trajectory between L1 and L2 from an exponential distribution
dPdecay/dl = 1/(b*ctau) Exp[-l/(b*ctau)]. 

Finally, choose a random decay event from the vN_Ntoall_inclDs_3.08517.txt decay file. 

Lorentz-transform the decay daughter products into the lab-frame along the LLP momentum direction with boost b. Place them at decay position (x,y,z), and assign this decay event the weight 
w''_i = w'_i * P_decay_i
which is just the number of LLP decays this decay event represents. 

5. Do whatever you do in your GEANT analysis to decide whether MATHUSLA can reconstruct this LLP decay i. If yes, count this LLP decay as w''_i observed LLP decays. 


6. 
FOR RECONSTRUCTION EFFICIENCY STUDIES: 
do the above, separate LLP decay events into groups (reconstructed) and (not reconstructed). 
efficiency = sum_i(reconstructed) w_i / sum_i(reconstructed + not reconstructed) w_i

This efficiency will have a significant dependence on mN, and will in general depend on LLP decay length ctau. However, in the long lifetime limit (say ctauX > ~ 1km, or Usq < 10^-10 for mN < 5 GeV), LLP decays are approximately uniformly distributed along their trajectories in MATHUSLA, and the lifetime will drop out of the efficiency. This efficiency in the 'long lifetime limit' is therefore a particularly interesting quantity to compute in simulations. 


FOR SIGNAL vs BACKGROUND and SENSITIVITY STUDIES:
For each mN, and for each Usq, obtain the list of observed LLP decays in our simulation, and check whether
sum_i w''_i >  4 (or whatever the minimum required number of observed LLP decays for our analysis is)






====================


A.1 DETAILS ON B, D MESON PRODUCTION
====================================

use public FONLL website to compute ds/dpt^2dy (pb/GeV^2)

# Job started on: Tue Mar  7 16:49:29 CET 2023 .
# FONLL heavy quark hadroproduction cross section, calculated on Tue Mar  7 16:49:36 CET 2023
# FONLL version and perturbative order: ## FONLL v1.3.2 fonll [ds/dpt^2dy (pb/GeV^2)]
# quark = bottom
# final state = meson. NP params (cm,lm,hm) = 24.2, 26.7, 22.2
# BR(q->meson) = 1
# ebeam1 = 7000, ebeam2 = 7000
# PDF set = CTEQ6.6



# Job started on: Tue Mar  7 16:59:27 CET 2023 .
# FONLL heavy quark hadroproduction cross section, calculated on Tue Mar  7 16:59:34 CET 2023
# FONLL version and perturbative order: ## FONLL v1.3.2 fonll [ds/dpt^2dy (pb/GeV^2)]
# quark = charm
# final state = meson (meson = 0.7 D0 + 0.3 D+). NP params (cm,lm,hm) = 0.1, 0.06, 0.135
# BR(q->meson) = 1
# ebeam1 = 7000, ebeam2 = 7000
# PDF set = CTEQ6.6

Given that LHC events are rotationally symmetric with respect to rotations around the beam axis, the binned dsigma/dpTdy differential cross section can therefore be used to construct a CDF and hence randomly generate B and D-meson momenta, which is what we did to get B- and D-meson 3-momentum-vectors. For D mesons, we assume 70% are D0 and 30% D+. For B-mesons, we assume 44.8% B0, 44.8% B+/-, 10.3% Bs0 and 0.018% Bc+/-. The B-meson fractions are taken from pythia8 hadronization, and are compatible with LHC measurements /910.09934. The resulting B- and D-meson 4-vectors were then decayed to various LLPs using full 2- and 3-body kinematics (phase space only). 


Note that FONLL computes cross sections to mean "average of q and qbar + X production", so to get actual number of B/D-mesons produced at LHC (1710.04921), we have to multiply the FONLL cross sections by 2. This is NOT included in the weighted LLP event generation, and has to be added (see description above). 


A.2 DETAILS ON TAU/W/Z PRODUCTION SIMULATIONS
==============================================

From comparing the Madgraph5 LO tau, W, Z production rate to measurements 1603.09222, we decided to use a universal electroweak K-factor of 1.3 to upscale the madgraph tau/W/Z production cross section for this analysis. 

Tau production: 

Used MG5_aMC 3.4.2, SM model at LO. 
generate p p > w+, w+ > ta+ vt
add process p p > w-, w- > ta- vt~
add process p p > ta+ ta-
with a 10 GeV lepton pT generator level cut. 
Showering in pythia8, then extracted undecayed tau 4-vectors from hadronized events. 


W/Z->N production:

Used MG5_aMC 3.4.2, SM_HeavyN_CKM_AllMasses_LO model at LO. 
e.g. for Ue-dominated case: 

generate p p > w+, w+ > e+ n1
add process p p > w-, w- > e- n1
add process p p > z, z > ve n1
add process p p > z, z > ve~ n1
output proc_RHN_Ue_DY

No generator level kinematic cuts were required, since we required on-shell intermediate z/w. 

Showering in pythia8, then extracted undecayed N 4-vectors from hadronized events. 



A.3 DETAILS ON DECAY EVENT GENERATION
=======================================

RHN Ue DECAY GUN AMMO

Use analytical 2- and 3-body decays below 0.42 gev

 
mN = 0.42 - 0.50 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = e+ e-
define Nv = ve ve~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nlep allq allq
output proc_RHN_Ue_ee_vN_Ntoall_lightfonly
 
 
 
mN = 0.50 - 0.99GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = e+ e-
define Nv = ve ve~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nlep allq allqs
output proc_RHN_Ue_ee_vN_Ntoall_nocharmnoss
 
 
mN = 0.99 - 1.871 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = e+ e-
define Nv = ve ve~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nlep allqs allqs
output proc_RHN_Ue_ee_vN_Ntoall_nocharm
 
 
mN = 1.871 - 1.97 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = e+ e-
define Nv = ve ve~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nlep allqs allqs
add process e+ e- > Nv n1, n1 > Nl allq c
add process e+ e- > Nv n1, n1 >  Nl allq c~
output proc_RHN_Ue_ee_vN_Ntoall_inclD
 
 
mN = 1.97 - 3.74 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = e+ e-
define Nv = ve ve~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nlep allqs allqs
add process e+ e- > Nv n1, n1 > Nl allqs c
add process e+ e- > Nv n1, n1 > Nl allqs c~
output proc_RHN_Ue_ee_vN_Ntoall_inclDs
 
 
mN > 3.74 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = e+ e-
define Nv = ve ve~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nlep allqs allqs
add process e+ e- > Nv n1, n1 > Nl allqs c
add process e+ e- > Nv n1, n1 > Nl allqs c~
add process e+ e- > Nv n1, n1 > Nv c c~
output proc_RHN_Ue_ee_vN_Ntoall_inclDD
 







RHN Umu DECAY GUN AMMO


Use analytical 2- and 3-body decays below 0.53 gev

 
mN = 0.53 - 0.60 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = mu+ mu-
define Nv = vm vm~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nlep allq allq
output proc_RHN_Umu_ee_vN_Ntoall_lightfonly
 
 
mN = 0.6- 0.99GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = mu+ mu-
define Nv = vm vm~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nlep allq allqs
output proc_RHN_Umu_ee_vN_Ntoall_nocharmnoss
 
 
mN = 0.99 - 1.98 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = mu+ mu-
define Nv = vm vm~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nlep allqs allqs
output proc_RHN_Umu_ee_vN_Ntoall_nocharm
 
 
 mN = 1.98 - 2.08 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = mu+ mu-
define Nv = vm vm~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nlep allqs allqs
add process e+ e- > Nv n1, n1 > Nl allq c
add process e+ e- > Nv n1, n1 >  Nl allq c~
output proc_RHN_Umu_ee_vN_Ntoall_inclD
 
 
mN = 2.08 - 3.74 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = mu+ mu-
define Nv = vm vm~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nlep allqs allqs
add process e+ e- > Nv n1, n1 > Nl allqs c
add process e+ e- > Nv n1, n1 > Nl allqs c~
output proc_RHN_Umu_ee_vN_Ntoall_inclDs
 
 
mN > 3.74 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = mu+ mu-
define Nv = vm vm~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nlep allqs allqs
add process e+ e- > Nv n1, n1 > Nl allqs c
add process e+ e- > Nv n1, n1 > Nl allqs c~
add process e+ e- > Nv n1, n1 > Nv c c~
output proc_RHN_Umu_ee_vN_Ntoall_inclDD





RHN Utau DECAY GUN AMMO


Use analytical 2- and 3-body decays below 0.42 gev
 
mN = 0.42 - 0.99 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = ta+ ta-
define Nv = vt vt~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nv allq allq
output proc_RHN_Utau_ee_vN_Ntoall_lightfonly
 
 
mN = 0.99 - 1.92 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = ta+ ta-
define Nv = vt vt~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nv allqs allqs
output proc_RHN_Utau_ee_vN_Ntoall_lightfsonly
 
 
 
mN = 1.92 - 2.28  GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = ta+ ta-
define Nv = vt vt~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nv allqs allqs
add process e+ e- > Nv n1, n1 > Nl allq allq
output proc_RHN_Utau_ee_vN_Ntoall_lightfstau
 
 
 
mN = 2.28 - 3.65  GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = ta+ ta-
define Nv = vt vt~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nv allqs allqs
add process e+ e- > Nv n1, n1 > Nl allq allqs
output proc_RHN_Utau_ee_vN_Ntoall_lightfstauK
 
 
 
mN = 3.65 - 3.75  GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = ta+ ta-
define Nv = vt vt~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nv allqs allqs
add process e+ e- > Nv n1, n1 > Nl allq allqc
output proc_RHN_Utau_ee_vN_Ntoall_lightfstauD
 
 
 
 
mN  > 3.75 GeV
import model SM_HeavyN_CKM_AllMasses_LO
define Nl = ta+ ta-
define Nv = vt vt~
define Nlep = Nl Nv
define alllep = e+ e- mu+ mu- ta+ ta- ve ve~ vm vm~ vt vt~
allq = u u~ d d~
allqs = u u~ d d~ s s~
allqc = u u~ d d~ s s~ c c~
generate e+ e- > Nv n1, n1 > Nl alllep alllep
add process e+ e- > Nv n1, n1 > Nv alllep alllep
add process e+ e- > Nv n1, n1 > Nv allqs allqs
add process e+ e- > Nv n1, n1 > Nl allq allqc
add process e+ e- > Nv n1, n1 > Nv c c~
add process e+ e- > Nv n1, n1 > Nl c s~
add process e+ e- > Nv n1, n1 > Nl s c~
output proc_RHN_Utau_ee_vN_Ntoall_lightfstauDD
 

