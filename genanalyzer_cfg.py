import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'ERROR'
#process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'file:../../decay0.root'
    #"/store/user/andrey/hzgamma_pythia8_153_8TeV_v2_HLT/hzgamma_pythia8_153_8TeV_v2_HLT/53f675467979b3dab12ab0598ae228db/hzgamma_pythia8_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_RECO_PU_102_1_Rkd.root"
#    "file:/uscms_data/d2/andreypz/cmssw/zgamma/generate/CMSSW_5_3_10/src/MCFM/Hadronizer_TuneZ2_8TeV_generic_LHE_pythia_cff_py_GEN_SIM.root"
    '/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_100_1_p96.root',
'/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_101_1_GfD.root',
'/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_102_1_y7A.root',
'/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_103_1_fOc.root',
'/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_104_1_TOw.root',
'/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_105_1_VhS.root',
'/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_106_1_kOi.root',
'/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_107_1_nbg.root',
'/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_108_1_7wT.root',
'/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_109_1_3JE.root',
'/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_10_1_j6H.root'
    )
)

process.demo = cms.EDAnalyzer('GenAnalyzer'
)


process.p = cms.Path(process.demo)
