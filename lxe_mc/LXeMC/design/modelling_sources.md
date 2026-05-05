Modelling sources:

Considering a source such as the cryostat. We can model the geometry as a shell, a disk, etc, in this case filled with Ti.

Our goal is to obtain the flux of potential gamma backgrounds, properly normalized, into the detector. The problem divides in the following steps.

1) which backgrounds?
  1.1 Gamma of 2445 keV from Bi-214. We can obtain the total number of gammas multiplying the late activity of U-238 by the Bi-214 gamma line BR (1.5 %). This process only produces one gamma per event.
  1.2 Gamma of 2.6 MeV from Tl-208. Here multiply Th-232 late for a BR of (check) 0.3%. Now we have the fraction of EVENTS where a 2.6 MeV gamma is present. But this gamma always can together with a second gamma. There are three main gammas that come with it, one about 585 keV is 85 % of the BR.

  Start with case 1.1 The Bi-214 gamma has an energy below Qbb, so it only enters in the ROI, even if deposits all its energy in the LXe in a SS, via energy resolution. If the gamma exits the source with reduced energy (say by 50 keV or more), it is guranteed not to enter the ROI.

  Therefore, one possible approach here would be to compute analytically the flux of gamma_bi in the face of the source facing the LXe. This effective flux should include self-shielding (which can be included in the integral), and take normalization into account (if we start with N g_bi, N/2 go in the oppsite direction of the detector and are lost). So from N g_bi we would end up with a dPhi/dNdu, where u= cos(theta), and N is properly normalized.

  The case of Tl-208 is more complex for two reasons. a) the companion gamma can make it together with the g_tl in the detector, and the event is vetoed (no two gammas allowed) if both gammas are visible (most of the time) and separated more than 3 mm in Z (also most of the times). But, on the other hand, g_tl almost never enters in the ROI if it escapes de Tl without loosing energy. But if it scatters inside the Tl and, say, a gamma of 2.5 MeV exists (and no companion gamma, e.g, companion going backwards), then this is an important background.

  It seems this cannot be modelled easily analytically, so I propose the following dedicated simulation.

  i) generate Tl-208 events according to activity x BR.
  ii) one gamma has 2.6 MeV. select direction of gamma via sampling (istropic); count and reject anyting going backwards;
  iii) Select by sampling one of the 3 possible companion gammas energy (look at the BR). Then again, sample direction.
  iv) propagate both gammas (and any extra gamma that comes out of the source), to a "virtual" skin, infinitely thin, right at the entrance of LXe.
  v) Use the fast KN propagation (no electrons)
  vi) Count how many gammas are in the LXe VS.
  vii) If only 1 gamma ask if energy > 100 keV. If not, reject; if yes, retain.
  viii) If two or more gammas in VS, project both of them to first interaction in LXe (via sampling). Discard gammas with E < 100 keV; If two gammas with Dz < 3 mm they count as one; If there are more than one gamma, discard. Otherwise, retain.
  iX) So one has now DPhi/dNdu in the virtual skin, for the potentia background gammas. Notice that in general energies can be less than g_tl (2.6), and actually we want to produce an histogram of those energies.

  The flux is only relevant in the face facing LXe. One can generate a large number of events (10^6 or 10^7) to get good statistics.

  X) Because we need to do this procedure with the Tl-208 gamma, it is probaby OK to do the same with Bi-214 gamma.

  XI) the fluxes must be fully normalized. Only a fraction of the initial N events will get into de solid angle and surive absorption.
