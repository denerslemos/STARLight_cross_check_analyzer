import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('Demo',eras.Run3_pp_on_PbPb)


process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("CondCore.CondDB.CondDB_cfi")
process.load('Configuration.EventContent.EventContent_cff')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True ),
)

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:/eos/cms/store/group/phys_heavyions/ddesouza/STARLight/prodType6/step3_RAW2DIGI_L1Reco_RECO.root'
                ),
                skipBadFiles=cms.untracked.bool(True),
				duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

#define number of events to be processed (-1 means all)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.GlobalTag.globaltag = cms.string('auto:phase1_2018_realistic_hi')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic_hi', '')
#process.HiForest.GlobalTagLabel = process.GlobalTag.globaltag

# Add PbPb centrality
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load('RecoHI.HiCentralityAlgos.HiCentrality_cfi')
process.hiCentrality.produceHFhits = False
process.hiCentrality.produceHFtowers = False
process.hiCentrality.produceEcalhits = False
process.hiCentrality.produceZDChits = False
process.hiCentrality.produceETmidRapidity = False
process.hiCentrality.producePixelhits = False
process.hiCentrality.produceTracks = False
process.hiCentrality.producePixelTracks = False
process.hiCentrality.reUseCentrality = True
process.hiCentrality.srcReUse = cms.InputTag("hiCentrality","","RECO")
process.centralityBin.Centrality = cms.InputTag("hiCentrality")
process.centralityBin.centralityVariable = cms.string("HFtowers")
process.centralityBin.nonDefaultGlauberModel = cms.string("")
process.GlobalTag.snapshotTime = cms.string("9999-12-31 23:59:59.000")
process.GlobalTag.toGet.extend([
    cms.PSet(record = cms.string("HeavyIonRcd"),
        tag = cms.string("CentralityTable_HFtowers200_HydjetDrum5F_v1032x02_mc"),#,"CentralityTable_HFtowers200_DataPbPb_periHYDJETshape_run2v1033p1x01_offline"),
        connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
        label = cms.untracked.string("HFtowers")
        ),
    ])
process.cent_seq = cms.Sequence(process.hiCentrality * process.centralityBin)


# Define the output

process.TFileService = cms.Service("TFileService",fileName = cms.string('PbPb_starlight_MC.root'))

process.demo = cms.EDAnalyzer('DemoAnalyzer',
                        vertexCollection  = cms.InputTag("offlinePrimaryVertices"),
                        tracks   = cms.InputTag("generalTracks"),
                        CentralitySrc    = cms.InputTag("hiCentrality"),
                        CentralityBinSrc = cms.InputTag("centralityBin","HFtowers"), 
                        towerSrc       = cms.InputTag("towerMaker")
)

process.p = cms.Path(process.cent_seq*process.demo)
process.schedule = cms.Schedule(process.p)
