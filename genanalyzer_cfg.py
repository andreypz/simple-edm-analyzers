import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'ERROR'
#process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'file:../../decay0.root'
    #'file:../../hzgamma_pythia6_py_GEN.root'
    #'file:../../hzgamma_dalitz_pythia8_py_GEN.root'
    #'file:../../hzgamma_stoyan_hack_pythia8_py_GEN.root'
    #'/store/user/andrey/HDalitz_el_stoyan_hack_v2/HDalitz_el_stoyan_hack_v2/ea86975e1b00840dbd43366bff35e68d/hzgamma_stoyan_hack_pythia8_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_RECO_PU_100_1_F0L.root',
    #'/store/user/andrey/HDalitz_el_stoyan_hack_v2/HDalitz_el_stoyan_hack_v2/ea86975e1b00840dbd43366bff35e68d/hzgamma_stoyan_hack_pythia8_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_RECO_PU_101_1_2bq.root'
    #"/store/user/andrey/hzgamma_pythia8_153_8TeV_v2_HLT/hzgamma_pythia8_153_8TeV_v2_HLT/53f675467979b3dab12ab0598ae228db/hzgamma_pythia8_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_RECO_PU_100_1_82E.root",
    #   "/store/user/andrey/hzgamma_pythia8_153_8TeV_v2_HLT/hzgamma_pythia8_153_8TeV_v2_HLT/53f675467979b3dab12ab0598ae228db/hzgamma_pythia8_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_RECO_PU_101_1_Tw0.root",
    "/store/user/andrey/hzgamma_pythia8_153_8TeV_v2_HLT/hzgamma_pythia8_153_8TeV_v2_HLT/53f675467979b3dab12ab0598ae228db/hzgamma_pythia8_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_RECO_PU_102_1_Rkd.root"
#    "file:/uscms_data/d2/andreypz/cmssw/zgamma/generate/CMSSW_5_3_10/src/MCFM/Hadronizer_TuneZ2_8TeV_generic_LHE_pythia_cff_py_GEN_SIM.root"

    )
)

process.demo = cms.EDAnalyzer('GenAnalyzer'
)


process.p = cms.Path(process.demo)
