{-# LANGUAGE RecordWildCards #-}

module HEP.Kinematics.Variable.M2 where

import HEP.Kinematics
import HEP.Kinematics.Vector.LorentzVector (setXYZT)

data InputKinematics = InputKinematics
                       { _pVis1  :: FourMomentum
                       , _pVis2  :: FourMomentum
                       , _ptmiss :: TransverseMomentum
                         -- ^ squared mass of the invisible particle
                       , _mInvSq :: !Double
                       } deriving Show

mkInvisibles :: InputKinematics
             -> (Double, Double, Double, Double)  -- ^ (k1x, k1y, k1z, k2z)
             -> (FourMomentum, FourMomentum)
mkInvisibles InputKinematics {..} (k1x, k1y, k1z, k2z) = (pInv1, pInv2)
  where
    eInv1 = sqrt (k1x * k1x + k1y * k1y + k1z * k1z + _mInvSq)
    pInv1 = eInv1 `seq` setXYZT k1x k1y k1z eInv1

    k2x = px _ptmiss - k1x
    k2y = py _ptmiss - k1y
    eInv2 = sqrt (k2x * k2x + k2y * k2y + k2z * k2z + _mInvSq)
    pInv2 = eInv2 `seq` setXYZT k2x k2y k2z eInv2

objFunc :: InputKinematics -> (FourMomentum, FourMomentum) -> Double
objFunc InputKinematics {..} (k1, k2) = max m1sq m2sq
  where
    m1sq = invariantMassSq [_pVis1, k1]
    m2sq = invariantMassSq [_pVis2, k2]
