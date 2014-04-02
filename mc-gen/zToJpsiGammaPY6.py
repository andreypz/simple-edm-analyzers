import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('runs py8'),
    name = cms.untracked.string('$Source: ,v $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('gen-sim.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'START53_V27::All', '')



from Configuration.Generator.PythiaUEZ2starSettings_cfi import *
from GeneratorInterface.ExternalDecays.TauolaSettings_cff import *
process.generator = cms.EDFilter("Pythia6GeneratorFilter",
                                 pythiaHepMCVerbosity = cms.untracked.bool(False),
                                 maxEventsToPrint = cms.untracked.int32(0),
                                 pythiaPylistVerbosity = cms.untracked.int32(1),
                                 filterEfficiency = cms.untracked.double(1.0),
                                 crossSection = cms.untracked.double(762.),
                                 comEnergy = cms.double(8000.0),
                                 ExternalDecays = cms.PSet(
    Tauola = cms.untracked.PSet(
    TauolaPolar,
    TauolaDefaultInputCards
    ),
    parameterSets = cms.vstring('Tauola')
    ),
                                 UseExternalGenerators = cms.untracked.bool(True),
                                 PythiaParameters = cms.PSet(
    pythiaUESettingsBlock,
    processParameters = cms.vstring('MSEL=0            !User defined processes',
                                    'MSUB(1)=1         !Incl Z0/gamma* production',
                                    'MSTP(43)=3        !Both Z0 and gamma*',
                                    'MDME(174,1)=0     !Z decay into d dbar',
                                    'MDME(175,1)=0     !Z decay into u ubar',
                                    'MDME(176,1)=0     !Z decay into s sbar',
                                    'MDME(177,1)=0     !Z decay into c cbar',
                                    'MDME(178,1)=0     !Z decay into b bbar',
                                    'MDME(179,1)=0     !Z decay into t tbar',
                                    'MDME(182,1)=1     !Z decay into e- e+',
                                    'MDME(183,1)=0     !Z decay into nu_e nu_ebar',
                                    # ########   Nasty hack implemented here  ###################
                                    
                                    'MDME(184,1)=1     !Z decay into mu- mu+',
                                    'KFDP(184,1)=443   ! one mu of the Z decay to J/Psi (hopefully)', 
                                    
                                    # #### End of nasty hack  ###
                                    'MDME(185,1)=0     !Z decay into nu_mu nu_mubar',
                                    'MDME(186,1)=0     !Z decay into tau- tau+',
                                    'MDME(187,1)=0     !Z decay into nu_tau nu_taubar',
                                    'CKIN(1)=50.       !Minimum sqrt(s_hat) value (=Z mass)'),
    # This is a vector of ParameterSet names to be read, in this order
    parameterSets = cms.vstring('pythiaUESettings',
                                'processParameters')
    )
                         )


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	#getattr(process,path)._seq = process.generator * getattr(process,path)._seq 
	getattr(process,path)._seq = process.generator 

