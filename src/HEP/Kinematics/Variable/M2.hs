{-# LANGUAGE RecordWildCards #-}

module HEP.Kinematics.Variable.M2 where

import HEP.Kinematics
import HEP.Kinematics.Vector.LorentzVector (setXYZT)

import Numeric.LinearAlgebra               (Vector, fromList, toList)

data InputKinematics = InputKinematics
                       { _pVis1  :: FourMomentum
                       , _pVis2  :: FourMomentum
                       , _ptmiss :: TransverseMomentum
                         -- | squared mass of the invisible particle
                       , _mInvSq :: !Double
                       } deriving Show

data Invisibles = Invisibles { _k1 :: FourMomentum
                             , _k2 :: FourMomentum } deriving Show

mkInvisibles :: InputKinematics -> [Double] -> Maybe Invisibles
mkInvisibles InputKinematics {..} [k1x, k1y, k1z, k2z] = do
    let eInv1 = sqrt (k1x * k1x + k1y * k1y + k1z * k1z + _mInvSq)
        k1 = eInv1 `seq` setXYZT k1x k1y k1z eInv1
        k2x = px _ptmiss - k1x
        k2y = py _ptmiss - k1y
        eInv2 = sqrt (k2x * k2x + k2y * k2y + k2z * k2z + _mInvSq)
        k2 = eInv2 `seq` setXYZT k2x k2y k2z eInv2
    return $ Invisibles k1 k2
mkInvisibles _ _ = Nothing

objFunc :: InputKinematics -> Vector Double -> (Double, Vector Double)
objFunc inp@InputKinematics {..} ks
    | length ks0 /= 4 = badResult
    | otherwise       =
          case objFunc' inp ks0 of
              Nothing                             -> badResult
              Just (m1sq, m2sq, Invisibles k1 k2) ->
                  let [k1x, k1y, k1z, k2z] = ks0
                  in if m1sq < m2sq
                     then let (eVis2, p2x, p2y, p2z) = epxpypz _pVis2
                              eInv2 = safeDivisor (energy k2)
                              r2 = eVis2 / eInv2
                              objDiff = [ 2 * (r2 * (k1x - px _ptmiss) + p2x)
                                        , 2 * (r2 * (k1y - py _ptmiss) + p2y)
                                        , 0
                                        , 2 * (r2 * k2z - p2z) ]
                          in (m2sq, fromList objDiff)
                     else let (eVis1, p1x, p1y, p1z) = epxpypz _pVis1
                              eInv1 = safeDivisor (energy k1)
                              r1 = eVis1 / eInv1
                              objDiff = [ 2 * (r1 * k1x - p1x)
                                        , 2 * (r1 * k1y - p1y)
                                        , 2 * (r1 * k1z - p1z)
                                        , 0 ]
                          in (m1sq, fromList objDiff)
  where
    ks0 = toList ks
    badResult = (1.0e+10, fromList [0,0,0,0])

objFunc' :: InputKinematics -> [Double] -> Maybe (Double, Double, Invisibles)
objFunc' inp@InputKinematics {..} ks = do
    invs@(Invisibles k1 k2) <- mkInvisibles inp ks
    let m1sq = invariantMassSq [_pVis1, k1]
        m2sq = invariantMassSq [_pVis2, k2]
    return (m1sq, m2sq, invs)

safeDivisor :: Double -> Double
safeDivisor x | x >= 0    = max eps x
              | otherwise = min (-eps) x
  where eps = 1.0e-8
