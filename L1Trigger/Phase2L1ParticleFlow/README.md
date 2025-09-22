# Particle-flow (PF) and PUPPI algorithms for L1T Phase-2

The PF and PUPPI algorithms are run in the [L1T Correlator Layer 1 Producer](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/python/l1ctLayer1_cff.py).

## PF Algorithm
The PF algorithm for L1 trigger has been designed to be simpler than the one offline, and able to process all input objects in parallel for pipelined execution on FPGAs.

It is [run](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/plugins/L1TCorrelatorLayer1Producer.cc#L642-L649) after all the tracks, muons, calorimeter and vertex inputs are taken from the CT layer 1.
The [PFAlgo3](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/src/pf/pfalgo3_ref.cpp#L286) algorithm is used.

Algorithm steps:
- [Track-Mu linking](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/src/pf/pfalgo_common_ref.cpp#L78):
    - Inside the tracker coverage, for each standalone muon it finds the closest track in $\Delta R$ and $\Delta p_T$. This track is then tagged as muon and excluded from further processing.

-  Linking of clusters in barrel and endcap is implemented differently, due to the different nature of calorimeters
    - In barrel:
        - [Track-EM linking](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/src/pf/pfalgo3_ref.cpp#L110):
            - For each track the closest EM cluster in $\Delta R$ is found. 
            - For each EM cluster the scalar sum of the $p_T$ of the associated tracks is performed ($\sum p_T^{track}$).
            - The EM cluster is now classified ad PF particle:
                - If there are no associated tracks, the EM clusters is classified as a photon.
                - if $p_T^{cluster} \ge \sum p_T^{track} -2 \sigma$, the tracks are tagged as electrons, and if $p_T^{cluster} - p_T^{track} \ge 1 \sigma$, the excess is classified as photon
                - If $p_T^{cluster} \ge \sum p_T^{track} -2 \sigma$, the cluster is likely originating from a hadronic shower starting in EG calorimeter. This is then not saved at this step.
        - The energy of linked EM clusters is [subtracted](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/src/pf/pfalgo3_ref.cpp#L244-L281) from the one of hadronic clusters. 
        - [Track- Hadronic Cluster Linking](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/src/pf/pfalgo3_ref.cpp#L411):
            - Tracks linked to closest hadronic cluster if $\Delta R<0.15$ and $p_T^{calo} \ge p_T^{track}-2\sigma_{calo}(p_T^{track})$
                - Tracks promoted to charged hadrons if linked to a cluster
            - Hadronic clusters classified as neutral hadrons if $(p_T^{calo})^2 - (\sum p_T^{track})^2 > \sum(\sigma_{calo}(p_T^{track}))^2$
                - The associated $p_T$ is $p_T^{calo} - \sum p_T^{track}$
    - In endcap:
        - Single linking step is performed. 
        - Tracks are linked to the closest cluster in $\Delta R$ and $\Delta p_T$. 
            - If the cluster has been identified as hadronic or electromagnetic, the cluster is promoted to charged hadrons or electrons respectively

## PUPPI Algorithm
PUPPI algorithm is [run](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/plugins/L1TCorrelatorLayer1Producer.cc#L653) after the PF algorithm and it is used to reduce the controbution of pileup. The [Linearized PUPPI](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/src/puppi/linpuppi_ref.cpp#L679-L711) algorithm is used. 

Within the coverage of the L1 track finder, PF delivers charged hadrons, neutral hadrons, charged electromagnetic, neutral electromagnetic, and muon particles. 

- The [charged PF candidates](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/src/puppi/linpuppi_ref.cpp#L275) are assigned to the PV if the [z-coordinate of the track is compatible with the vertex](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/src/puppi/linpuppi_ref.cpp#L296).
    - If the [NNVtx track-to-vertex association is on](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/python/l1ctLayer1_cff.py#L12), the charged PF candidates are kept only if they are associated to the PV with the [MVS score](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/src/NNVtxAssoc.cc#L29)
- For the [neutral PF candidates](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/src/puppi/linpuppi_ref.cpp#L691) a weight $w$ is computed:
    - First the metric $\alpha_C$ is computed, with a sum over all the tracks from the PV within a cone $\Delta R<0.3$ from the neutral particle candidate.
        - $\alpha_C = log \sum \frac{min(p_T^i,p_T^{max})^2}{max(\Delta R, \Delta R^{min})^2}$  
        - $p_T^{max} = 50 GeV$; \Delta R ^{min}=0.07 (0.04) in barrel (endcap)$
    - The final particle weight $w$ is computed combining $\alpha_c$ with the per-event pileup estimate
        - $w = \frac{1}{1+ e^{-x}}$
        - $x = x_{\alpha} + x_{p_T} - x_{PU}$
        - $x_a = min(max(s_{\alpha}\cdot (\alpha-\alpha_0),-\alpha^{max},\alpha^{max}))$
        - $x_{p_T} = s_{p_{T}} \cdot (p_T - p_T^0)$
        - $\alpha = log(\alpha_C)$
    - The final weight $w$ is applied to the $p_T$ of the neutral particle, giving a corrected $p_T' = w p_T$.
    - A [cut](https://github.com/p2l1pfp/cmssw/blob/L1PF_15_1_X/L1Trigger/Phase2L1ParticleFlow/src/puppi/linpuppi_ref.cpp#L532) on the weighted $p_T$ of the neutral particles is applied.
        - PUPPI algorithm mitigates the PU effect and reduces the PF candidates by a factor of ~10.


### References
- [TDR: The Phase-2 Upgrade of the CMS Level-1 Trigger
](https://cds.cern.ch/record/2714892?ln=it)
    - 3.5.3 
