import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'ERROR'
#process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:mc-gen/gen-sim-run.root'

    #'file:mc-gen/gen-sim-zjpsigamma-4-mumu.root'
    #'/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_100_1_p96.root',
    #'/store/user/andrey/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/Higgs_to_JpsiGamma_MH125_8TeV_PYTHIA8/78f58012d01406e76b9b93e195d19b7e/gen-sim_101_1_GfD.root',

    #'/store/user/lpchzg/prod/ZtoJPsiGamma/andrey/ZtoJPsiGamma-MuMuGamma-Pythia6/ZtoJPsiGamma-MuMuGamma-Pythia6/263cb8f1f99ea7cc4472c989c939917e/gen-sim_98_1_qSf.root',
    #'/store/user/lpchzg/prod/ZtoJPsiGamma/andrey/ZtoJPsiGamma-MuMuGamma-Pythia6/ZtoJPsiGamma-MuMuGamma-Pythia6/263cb8f1f99ea7cc4472c989c939917e/gen-sim_99_1_plm.root',
    )
)

process.demo = cms.EDAnalyzer('GenAnalyzer'
)


process.p = cms.Path(process.demo)
