{-# LANGUAGE RecordWildCards #-}

module HEP.Kinematics.Variable.M2 where

import HEP.Kinematics
import HEP.Kinematics.Vector.LorentzVector (setXYZT)

import Numeric.LinearAlgebra               (Vector, fromList, toList)

data InputKinematics = InputKinematics
                       { _p1     :: FourMomentum
                       , _p2     :: FourMomentum
                       , _q1     :: FourMomentum
                       , _q2     :: FourMomentum
                       , _ptmiss :: TransverseMomentum
                         -- | squared mass of the invisible particle
                       , _mInvSq :: !Double
                       }

mkInput :: [FourMomentum]
        -> [FourMomentum]
        -> TransverseMomentum
        -> Double
        -> Maybe InputKinematics
mkInput as bs ptmiss mInv
    | length as /= 2 || length bs /= 2 = Nothing
    | otherwise                        = do
          let ps = zipWith (+) as bs
          return $ InputKinematics { _p1 = head ps
                                   , _p2 = ps !! 1
                                   , _q1 = head bs
                                   , _q2 = bs !! 1
                                   , _ptmiss = ptmiss
                                   , _mInvSq = mInv * mInv
                                   }

data Invisibles = Invisibles { _k1 :: FourMomentum , _k2 :: FourMomentum }

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

objFunc :: Maybe InputKinematics -> Vector Double -> (Double, Vector Double)
objFunc (Just inp@InputKinematics {..}) ks =
    let ks0 = toList ks
    in if length ks0 /= 4
       then badOutput
       else case objFunc' inp ks0 of
                Nothing                             -> badOutput
                Just (m1sq, m2sq, Invisibles k1 k2) ->
                    let [k1x, k1y, k1z, k2z] = ks0
                    in if m1sq < m2sq
                       then let (eVis2, p2x, p2y, p2z) = epxpypz _p2
                                eInv2 = safeDivisor (energy k2)
                                r2 = eVis2 / eInv2
                                objDiff = [ 2 * (r2 * (k1x - px _ptmiss) + p2x)
                                          , 2 * (r2 * (k1y - py _ptmiss) + p2y)
                                          , 0
                                          , 2 * (r2 * k2z - p2z) ]
                            in (m2sq, fromList objDiff)
                       else let (eVis1, p1x, p1y, p1z) = epxpypz _p1
                                eInv1 = safeDivisor (energy k1)
                                r1 = eVis1 / eInv1
                                objDiff = [ 2 * (r1 * k1x - p1x)
                                          , 2 * (r1 * k1y - p1y)
                                          , 2 * (r1 * k1z - p1z)
                                          , 0 ]
                            in (m1sq, fromList objDiff)
objFunc _ _ = badOutput

badOutput :: (Double, Vector Double)
badOutput = (1.0e+10, fromList [0,0,0,0])

objFunc' :: InputKinematics -> [Double] -> Maybe (Double, Double, Invisibles)
objFunc' inp@InputKinematics {..} ks = do
    invs@(Invisibles k1 k2) <- mkInvisibles inp ks
    let m1sq = invariantMassSq [_p1, k1]
        m2sq = invariantMassSq [_p2, k2]
    return (m1sq, m2sq, invs)

safeDivisor :: Double -> Double
safeDivisor x | x >= 0    = max eps x
              | otherwise = min (-eps) x
  where eps = 1.0e-8
