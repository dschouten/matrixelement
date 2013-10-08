
#
# define the recoil as a sum of jets above threshold
#

sumjets = tlv()

for j in measured_jets:
    sumjets += j

mehndl.measured.recoil.SetPxPyPzE( j.Px(), j.Py(), 0, j.Pt() ) ## note the sign !! recoil = -higgs pT

log.info( 'recoil (estimate): %0.2f, %0.2f'%( mehndl.measured.recoil.Pt(),
                                              mehndl.measured.recoil.Phi() ) )
