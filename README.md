Simple EDM analyzers
-------------------

This repository contains various simple EDM analyzers as well as MC Generation scripts to run within CMSSW.


 * *eleanalyzer_cfg* and corresponding *src/EleAnalyzer.cc*
This code is able to run over AOD/RECO sample, loop over a collection of gsfElectrons and make simple distributions. 

 * *genanalyzer_cfg* and *GenAnalyzer.cc*
It runs over GEN, or GENSIM Monte Carlo sample and analyses genparticles. Mainly looking for Z/Higgs, photons and muons


 * *mc-gen*: step_0,1,2 are the scripts needed for the hadronization the full processing of H->l+l-gamma samle starting from LHE file (which was prepared in Madgraph beforehand).
 * *mc-gen*: zToJpsiGamma and hToJpsiGamma are the ways to produce  h->J/Psi+gamma and z->J/Psi+gamma processes.
