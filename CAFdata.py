import os

#Data
DATA_DIR = '/pnfs/sbnd/persistent/sbndpro/mcp/mc/official/MCP2022A/prodoverlay_corsika_cosmics_proton_genie_rockbox_sce/v09_37_02_04/reco2_caf/'
FILES = os.listdir(f'{DATA_DIR}')
FLAT_FILES = [file for file in FILES if 'flat' in file]

#Key labels
hdrkeys = ['rec.hdr.run','rec.hdr.subrun','rec.hdr.evt'] #Header keys
primprefix = 'rec.mc.nu.prim.' #Add this to beginning of string for flat cafana for primary particls
mcnuprefix = 'rec.mc.nu.' #Add this to beginning of string for flat cafana for neutrino events
mcslcprefix = 'rec.slc.truth.' #Add this to beginning of string for flat cafana for neutrino slice
recoslcprefix = 'rec.slc.reco.' #Add this to beginning of string for flat cafana for neutrino slice
recoprefix = 'rec.reco.' #Add this to beginning of string for flat cafana for neutrino events
shwprefix = recoslcprefix+'shw.'
trkprefix = recoslcprefix+'trk.'

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
  'position.z']
mcslc_keys = [mcslcprefix+key for key in mcslc_keys]
mcprim_keys = ['pdg',
  'gstatus',
  'genT',
  'genE',
  'genp.x',
  'genp.y',
  'genp.z']
mcprim_keys = [primprefix+key for key in mcprim_keys]