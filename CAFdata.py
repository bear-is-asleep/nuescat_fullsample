import os
import cuts
import numpy as np

#Data
DATA_DIR = '/pnfs/sbnd/persistent/sbndpro/mcp/mc/official/MCP2022A/prodoverlay_corsika_cosmics_proton_genie_rockbox_sce/v09_37_02_04/reco2_caf/'
FILES = os.listdir(f'{DATA_DIR}')
FLAT_FILES = [file for file in FILES if 'flat' in file]
DEFNAME = ''
SBND_AREA = 4*4*1e6 #cm^2


#Cut range searches
etheta_range = np.arange(0.0001,0.01,step=0.001)
len_range = np.arange(2,400,step=5)
trk_eng_range = np.arange(0.001,1,step=0.01)
#bdt_score_range = np.arange(0.02,1.02,step=0.02)
trk_npts_range = np.arange(1,50,step=1)
trk_angle_range = np.arange(1e-3,np.pi,step=0.05)
openangle_range = np.arange(1e-5,45,step=1)
bdt_score_range = np.arange(0.02,1.02,step=0.02)
dedx_range = np.arange(1,10,step=0.1)

"""
Selection dictionary has the cut key as the key and in it is
cut_range : the range of cut values to test
cut_function : what type of cut
cut_key : which key it is in the df
equality : cut on values lessthan, greaterthan, equal? - keeps value (equality) compared to cut val

Next 3 are for NReco cut only, set to None if using a box cut
nreco_key : which key to check for nreco objects, i.e. nshw
nreco : number of reco objects, i.e. nshw = 1 would be 1
nreco_equality : check for reco objects lessthan, greaterthan, equal, i.e. nshw = 1 would be equal
"""
nue_sel_dict = {
  'cRVtxz': {
    'cut_range':[25],
    'cut_function':cuts.cBoxCut,
    'cut_key':'reco_vtx.z',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cNShwG': {
    'cut_range':[0],
    'cut_function':cuts.cBoxCut,
    'cut_key':'nshw',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cNShwL': {
    'cut_range':[3],
    'cut_function':cuts.cBoxCut,
    'cut_key':'nshw',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cNTrk': {
    'cut_range':[2],
    'cut_function':cuts.cBoxCut,
    'cut_key':'ntrk',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cNStub': {
    'cut_range':[4],
    'cut_function':cuts.cBoxCut,
    'cut_key':'nstub',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cRTheta': {
    'cut_range':np.arange(0.05,0.65,step=0.05),
    'cut_function':cuts.cBoxCut,
    'cut_key':'reco_theta',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLTrkLen': {
    'cut_range':[10],
    'cut_function':cuts.cBoxCut,
    'cut_key':'ltrk.len',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  # 'cLTrkEng': {
  #   'cut_range':[0.0075],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'ltrk.eng',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cLTrkNpts': {
  #   'cut_range':[25],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'ltrk.npts',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cLTrkPfpScore': {
  #   'cut_range':[0.6],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'ltrk.pfptrkscore',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cLTrkMuonScore': {
  #   'cut_range':[0.5],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'ltrk.muonscore',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cLTrkPionScore': {
  #   'cut_range':[0.4],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'ltrk.pionscore',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cLTrkProtonScore': {
  #   'cut_range':[0.4],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'ltrk.protonscore',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cLTrkAngle': {
  #   'cut_range':trk_angle_range,
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'ltrk.angle',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  'cLTrkEtheta': {
    'cut_range':[0.001,0.002,0.004],
    'cut_function':cuts.cBoxCut,
    'cut_key':'ltrk.Etheta',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  # 'cSLTrkLen': {
  #   'cut_range':[0],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'sltrk.len',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cSLTrkEng': {
  #   'cut_range':[0.025],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'sltrk.eng',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cSLTrkNpts': {
  #   'cut_range':[0],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'sltrk.npts',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cSLTrkPfpScore': {
  #   'cut_range':[0.7],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'sltrk.pfptrkscore',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cSLTrkMuonScore': {
  #   'cut_range':[0.1],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'sltrk.muonscore',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cSLTrkPionScore': {
  #   'cut_range':[0.5],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'sltrk.pionscore',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cSLTrkProtonScore': {
  #   'cut_range':[0.5],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'sltrk.protonscore',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cSLTrkEtheta': {
  #   'cut_range':[0.001],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'sltrk.Etheta',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  'cSLShwOpenAngle': {
    'cut_range':[5],
    'cut_function':cuts.cBoxCut,
    'cut_key':'slshw.openangle',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLShwConvGap': {
    'cut_range':[5],
    'cut_function':cuts.cBoxCut,
    'cut_key':'slshw.cnvgap',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  # 'cSLShwdedx': {
  #   'cut_range':[3],
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'slshw.dedx',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  'cSLShwElectron': {
    'cut_range':[0.3],
    'cut_function':cuts.cBoxCut,
    'cut_key':'slshw.electron',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLShwEtheta': {
    'cut_range':[0.001],
    'cut_function':cuts.cBoxCut,
    'cut_key':'slshw.Etheta',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  #LShw
  'cLShwOpenAngle': {
    'cut_range':[10],
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.openangle',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLShwConvGap': {
    'cut_range':[5],
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.cnvgap',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLShwElectron': {
    'cut_range':[0.5],
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.electron',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLShwdedx': {
    'cut_range':[3],
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.dedx',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLShwEtheta': {
    'cut_range':[0.004],
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.Etheta',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
}

"""
Search dictionary has the cut key as the key and in it is
cut_range : the range of cut values to test cut
cut_function : what type of cut
cut_key : which key it is in the df
equality : cut on values lessthan, greaterthan, equal? - keeps value (equality) compared to cut val

Next 3 are for NReco cut only, set to None if using a box cut
nreco_key : which key to check for nreco objects, i.e. nshw
nreco : number of reco objects, i.e. nshw = 1 would be 1
nreco_equality : check for reco objects lessthan, greaterthan, equal, i.e. nshw = 1 would be equal
"""

nue_search_dict = {
  'cRVtxz': {
    'cut_range':np.arange(20,50,step=5),
    'cut_function':cuts.cBoxCut,
    'cut_key':'reco_vtx.z',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cRTheta': {
    'cut_range':trk_angle_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'reco_theta',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cNShwG': {
    'cut_range':[-1,0],
    'cut_function':cuts.cBoxCut,
    'cut_key':'nshw',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cNShwL': {
    'cut_range':[2,3,4],
    'cut_function':cuts.cBoxCut,
    'cut_key':'nshw',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cNTrk': {
    'cut_range':[1,2,3],
    'cut_function':cuts.cBoxCut,
    'cut_key':'ntrk',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cNStub': {
    'cut_range':[2,3,4,5,6,7],
    'cut_function':cuts.cBoxCut,
    'cut_key':'nstub',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLTrkLen': {
    'cut_range':len_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'ltrk.len',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLTrkEng': {
    'cut_range':trk_eng_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'ltrk.eng',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLTrkNpts': {
    'cut_range':trk_npts_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'ltrk.npts',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLTrkPfpScore': {
    'cut_range':bdt_score_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'ltrk.pfptrkscore',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLTrkMuonScore': {
    'cut_range':bdt_score_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'ltrk.muonscore',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLTrkPionScore': {
    'cut_range':bdt_score_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'ltrk.pionscore',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLTrkProtonScore': {
    'cut_range':bdt_score_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'ltrk.protonscore',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLTrkAngle': {
    'cut_range':trk_angle_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'ltrk.angle',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  # 'cLTrkdedxG': {
  #   'cut_range':dedx_range,
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'ltrk.dedx',
  #   'equality':'greaterthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cLTrkdedxL': {
  #   'cut_range':dedx_range,
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'ltrk.dedx',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  'cLTrkEtheta': {
    'cut_range':etheta_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'ltrk.Etheta',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  #SLTrk
  'cSLTrkLen': {
    'cut_range':len_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'sltrk.len',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLTrkEng': {
    'cut_range':trk_eng_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'sltrk.eng',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLTrkNpts': {
    'cut_range':trk_npts_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'sltrk.npts',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLTrkPfpScore': {
    'cut_range':bdt_score_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'sltrk.pfptrkscore',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLTrkMuonScore': {
    'cut_range':bdt_score_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'sltrk.muonscore',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLTrkPionScore': {
    'cut_range':bdt_score_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'sltrk.pionscore',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLTrkProtonScore': {
    'cut_range':bdt_score_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'sltrk.protonscore',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  # 'cSLTrkdedxG': {
  #   'cut_range':dedx_range,
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'sltrk.dedx',
  #   'equality':'greaterthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  # 'cSLTrkdedxL': {
  #   'cut_range':dedx_range,
  #   'cut_function':cuts.cBoxCut,
  #   'cut_key':'sltrk.dedx',
  #   'equality':'lessthan',
  #   'nreco_key':None,
  #   'nreco':None,
  #   'nreco_equality':None
  # },
  'cSLTrkEtheta': {
    'cut_range':etheta_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'sltrk.Etheta',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  #SLShw
  'cSLShwdedxG': {
    'cut_range':dedx_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'slshw.dedx',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLShwdedxL': {
    'cut_range':dedx_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'slshw.dedx',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLShwlenG': {
    'cut_range':len_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'slshw.len',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLShwlenL': {
    'cut_range':len_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'slshw.len',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLShwOpenAngle': {
    'cut_range':openangle_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'slshw.openangle',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLShwConvGap': {
    'cut_range':trk_npts_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'slshw.cnvgap',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLShwElectron': {
    'cut_range':bdt_score_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'slshw.electron',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cSLShwEtheta': {
    'cut_range':etheta_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'slshw.Etheta',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  #LShw
  'cLShwdedxG': {
    'cut_range':dedx_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.dedx',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLShwdedxL': {
    'cut_range':dedx_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.dedx',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLShwlenG': {
    'cut_range':len_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.len',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLShwlenL': {
    'cut_range':len_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.len',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLShwOpenAngle': {
    'cut_range':trk_npts_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.openangle',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLShwConvGap': {
    'cut_range':trk_npts_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.cnvgap',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLShwElectron': {
    'cut_range':bdt_score_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.electron',
    'equality':'greaterthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
  'cLShwEtheta': {
    'cut_range':etheta_range,
    'cut_function':cuts.cBoxCut,
    'cut_key':'lshw.Etheta',
    'equality':'lessthan',
    'nreco_key':None,
    'nreco':None,
    'nreco_equality':None
  },
}

#Key labels
hdrkeys = ['rec.hdr.run','rec.hdr.subrun','rec.hdr.evt'] #Header keys
primprefix = 'rec.mc.nu.prim.' #Add this to beginning of string for flat cafana for primary particls
mcnuprefix = 'rec.mc.nu.' #Add this to beginning of string for flat cafana for neutrino events
mcslcprefix = 'rec.slc.truth.' #Add this to beginning of string for flat cafana for neutrino slice
recoslcprefix = 'rec.slc.reco.' #Add this to beginning of string for flat cafana for neutrino slice
recoprefix = 'rec.reco.' #Add this to beginning of string for flat cafana for neutrino events
shwprefix = recoslcprefix+'shw.'
trkprefix = recoslcprefix+'trk.'
primsliceprefix = mcslcprefix+'prim'

#Write down CAF keys we want to use
nreco_keys = ['nshw',
  'ntrk',
  'nstub']
nreco_keys = [recoprefix + key for key in nreco_keys]
shw_keys = ['razzle.electronScore',
  'bestplane_energy',
  'dir.x',
  'dir.y',
  'dir.z',
  'razzle.pdg',
  'pfp.slcID']
shw_keys = [shwprefix + key for key in shw_keys]
trk_keys = ['pfp.slcID',
  'dir.x',
  'dir.y',
  'dir.z',
  'best_plane',
  'calo.0.ke',
  'calo.1.ke',
  'calo.2.ke']
trk_keys = [trkprefix + key for key in trk_keys]
mcnu_keys = ['iscc',
  'position.x',
  'position.y',
  'position.z',
  'pdg',
  'npizero']
mcnu_keys = [mcnuprefix+key for key in mcnu_keys]
mcslc_keys = ['pdg',
  'npizero',
  'parent_pdg',
  'position.x',
  'position.y',
  'position.z',
  'hitnuc',
  'nprim']
mcslc_keys = [mcslcprefix+key for key in mcslc_keys]
mcprim_keys = ['pdg',
  'gstatus',
  'genT',
  'genE',
  'genp.x',
  'genp.y',
  'genp.z',
  'truth.E']
mcprim_keys = [primprefix+key for key in mcprim_keys]
slcprim_keys = ['pdg',
  'gstatus',
  'genT',
  'genE',
  'genp.x',
  'genp.y',
  'genp.z',
  'truth.E']
slcprim_keys = [primsliceprefix+key for key in mcprim_keys]