import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'/store/user/andrey/MCFM_lord_hzgamma_8TeV_LHE_pythia6_v2/AODSIM/39bf61f738ba3bdb8860f0848073cc88/aodsim_100_1_BGG.root'
        
    )
)

process.demo = cms.EDAnalyzer('EleAnalyzer'
)


process.p = cms.Path(process.demo)
