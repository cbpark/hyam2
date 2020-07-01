{-# LANGUAGE RecordWildCards #-}

module HEP.Kinematics.Variable.M2 where

import HEP.Kinematics
import HEP.Kinematics.Vector.LorentzVector (setXYZT)

import Numeric.LinearAlgebra               (Vector, fromList, toList)

data InputKinematics = InputKinematics
                       { _p1     :: FourMomentum  -- ^ p1 = a1 + b1
                       , _p2     :: FourMomentum  -- ^ p2 = a2 + b2
                       , _q1     :: FourMomentum  -- ^ q1 = b1
                       , _q2     :: FourMomentum  -- ^ q2 = b2
                       , _ptmiss :: TransverseMomentum
                         -- | squared mass of the invisible particle
                       , _mInvSq :: !Double
                       }

-- | creates the input for A1 + A2 --> a1 B1 + a2 B2 --> a1 b1 C1 + a2 b2 C2.
mkInput :: [FourMomentum]      -- ^ [a1, a2]
        -> [FourMomentum]      -- ^ [b1, b2]
        -> TransverseMomentum  -- ^ pT(miss)
        -> Double              -- ^ M_{invisible}
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

-- | (k1x, k1y, k1z, k2z).
type Variables = (Double, Double, Double, Double)

data Invisibles = Invisibles { _k1 :: FourMomentum , _k2 :: FourMomentum }

mkInvisibles :: InputKinematics -> Variables -> Invisibles
mkInvisibles InputKinematics {..} (k1x, k1y, k1z, k2z) =
    let eInv1 = sqrt (k1x * k1x + k1y * k1y + k1z * k1z + _mInvSq)
        k1 = eInv1 `seq` setXYZT k1x k1y k1z eInv1
        k2x = px _ptmiss - k1x
        k2y = py _ptmiss - k1y
        eInv2 = sqrt (k2x * k2x + k2y * k2y + k2z * k2z + _mInvSq)
        k2 = eInv2 `seq` setXYZT k2x k2y k2z eInv2
    in Invisibles k1 k2

-- | the unknowns are (k1x, k2x, k1z, k2z).
objFunc :: Maybe InputKinematics -> Vector Double -> (Double, Vector Double)
objFunc (Just inp@InputKinematics {..}) ks
    | length ks0 /= 4 = badOutput
    | otherwise       =
      let [k1x, k1y, k1z, k2z] = ks0
          vars = (k1x, k1y, k1z, k2z)
          (m1sq, m2sq, Invisibles k1 k2) = objFunc' inp vars
          derivs = m2Diff _p1 _p2 k1 k2 _ptmiss vars
      in if m1sq < m2sq
         then (m2sq, snd derivs)
         else (m1sq, fst derivs)
  where ks0 = toList ks
objFunc _ _ = badOutput

objFunc' :: InputKinematics -> Variables -> (Double, Double, Invisibles)
objFunc' inp@InputKinematics {..} ks =
    let invs@(Invisibles k1 k2) = mkInvisibles inp ks
        m1sq = invariantMassSq [_p1, k1]
        m2sq = invariantMassSq [_p2, k2]
    in (m1sq, m2sq, invs)

badOutput :: (Double, Vector Double)
badOutput = (1.0e+10, fromList [0,0,0,0])

m2Diff :: FourMomentum  -- ^ p1 (or q1)
       -> FourMomentum  -- ^ p2 (or q2)
       -> FourMomentum  -- ^ k1
       -> FourMomentum  -- ^ k2
       -> TransverseMomentum
       -> Variables
       -> (Vector Double, Vector Double)
m2Diff p1 p2 k1 k2 ptmiss (k1x, k1y, k1z, k2z) =
    let (e1, p1x, p1y, p1z) = epxpypz p1
        r1 = e1 / safeDivisor (energy k1)
        d1 = fromList [ 2 * (r1 * k1x - p1x)
                      , 2 * (r1 * k1y - p1y)
                      , 2 * (r1 * k1z - p1z)
                      , 0 ]

        (e2, p2x, p2y, p2z) = epxpypz p2
        r2 = e2 / safeDivisor (energy k2)
        d2 = fromList [ 2 * (r2 * (k1x - px ptmiss) + p2x)
                      , 2 * (r2 * (k1y - py ptmiss) + p2y)
                      , 0
                      , 2 * (r2 * k2z - p2z) ]
    in (d1, d2)

safeDivisor :: Double -> Double
safeDivisor x | x >= 0    = max eps x
              | otherwise = min (-eps) x
  where eps = 1.0e-8
