import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing


options = VarParsing.VarParsing ('analysis')
options.maxEvents = -1
options.loadFromFile('inputFiles','dalitz_short.txt')
options.parseArguments()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#process.load("RecoTracker/FinalTrackSelectors/selectHighPurity_cfi")
#process.selectHighPurity.copyExtras = cms.untracked.bool(False)
#process.selectHighPurity.vertices = cms.InputTag("offlinePrimaryVertices")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(options.inputFiles)
                            )

process.demo = cms.EDAnalyzer('EleAnalyzer'
)


process.p = cms.Path(
    process.demo)
