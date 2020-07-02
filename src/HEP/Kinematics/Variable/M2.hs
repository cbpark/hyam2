{-# LANGUAGE RecordWildCards #-}

module HEP.Kinematics.Variable.M2 where

import HEP.Kinematics
import HEP.Kinematics.Vector.LorentzVector (setXYZT)

import Numeric.LinearAlgebra               (Vector, fromList, toList)
import Numeric.NLOPT

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

m2SQP :: Maybe InputKinematics -> Either Result Solution
m2SQP Nothing                         = Left INVALID_ARGS
m2SQP (Just inp@InputKinematics {..}) =
    let objfD = m2ObjF inp

        c1fD = constraintA inp
        c1 = EqualityConstraint (Scalar c1fD) 1e-6
        c2fD = constraintB inp
        c2 = EqualityConstraint (Scalar c2fD) 1e-6

        stop = ObjectiveRelativeTolerance 1e-9 :| []
        algorithm = SLSQP objfD [] [] [c1, c2]
        problem = LocalProblem 4 stop algorithm

        x0 = fromList [0.5 * px _ptmiss, 0.5 * py _ptmiss, 0, 0]
        sol = minimizeLocal problem x0
    in sol

-- | (k1x, k1y, k1z, k2z).
type Variables = (Double, Double, Double, Double)

data Invisibles = Invisibles { _k1 :: FourMomentum , _k2 :: FourMomentum }

mkInvisibles :: InputKinematics -> Variables -> Invisibles
mkInvisibles InputKinematics {..} (k1x, k1y, k1z, k2z) = Invisibles k1 k2
  where
    eInv1 = sqrt (k1x * k1x + k1y * k1y + k1z * k1z + _mInvSq)
    k1 = eInv1 `seq` setXYZT k1x k1y k1z eInv1
    k2x = px _ptmiss - k1x
    k2y = py _ptmiss - k1y
    eInv2 = sqrt (k2x * k2x + k2y * k2y + k2z * k2z + _mInvSq)
    k2 = eInv2 `seq` setXYZT k2x k2y k2z eInv2

-- | the unknowns are (k1x, k2x, k1z, k2z).
m2ObjF :: InputKinematics -> Vector Double -> (Double, Vector Double)
m2ObjF inp@InputKinematics {..} ks =
    if m1sq < m2sq then (m2sq, grad2) else (m1sq, grad1)
  where
    [k1x, k1y, k1z, k2z] = toList ks
    vars = (k1x, k1y, k1z, k2z)
    (m1sq, m2sq, Invisibles k1 k2) = m2ObjF' inp vars
    (grad1, grad2, _) = m2Grad _p1 _p2 k1 k2 _ptmiss vars

m2ObjF' :: InputKinematics -> Variables -> (Double, Double, Invisibles)
m2ObjF' inp@InputKinematics {..} ks = (m1sq, m2sq, invs)
  where
    invs@(Invisibles k1 k2) = mkInvisibles inp ks
    m1sq = invariantMassSq [_p1, k1]
    m2sq = invariantMassSq [_p2, k2]

m2Grad :: FourMomentum  -- ^ p1 (or q1)
       -> FourMomentum  -- ^ p2 (or q2)
       -> FourMomentum  -- ^ k1
       -> FourMomentum  -- ^ k2
       -> TransverseMomentum
       -> Variables
       -> (Vector Double, Vector Double, Double)
m2Grad p1 p2 k1 k2 ptmiss (k1x, k1y, k1z, k2z) = (d1, d2, deltaMsq)
  where
    (e1, p1x, p1y, p1z) = epxpypz p1
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

    deltaMsq = invariantMassSq [p1 + k1] - invariantMassSq [p2 + k2]

constraintA :: InputKinematics -> Vector Double -> (Double, Vector Double)
constraintA inp@InputKinematics {..} ks
    | length ks0 /= 4 = badOutput
    | otherwise       =
      let [k1x, k1y, k1z, k2z] = ks0
          vars = (k1x, k1y, k1z, k2z)
          (m1sq, m2sq, Invisibles k1 k2) = m2ObjF' inp vars
          (grad1, grad2, _) = m2Grad _p1 _p2 k1 k2 _ptmiss vars
      in (m1sq - m2sq, grad1 - grad2)
  where ks0 = toList ks

constraintB :: InputKinematics -> Vector Double -> (Double, Vector Double)
constraintB inp@InputKinematics {..} ks
    | length ks0 /= 4 = badOutput
    | otherwise       =
      let [k1x, k1y, k1z, k2z] = ks0
          vars = (k1x, k1y, k1z, k2z)
          (_, _, Invisibles k1 k2) = m2ObjF' inp vars
          (grad1, grad2, deltaMsq) = m2Grad _q1 _q2 k1 k2 _ptmiss vars
      in (deltaMsq, grad1 - grad2)
  where ks0 = toList ks

safeDivisor :: Double -> Double
safeDivisor x | x >= 0    = max eps x
              | otherwise = min (-eps) x
  where eps = 1.0e-8

badOutput :: (Double, Vector Double)
badOutput = (1.0e+10, fromList [0,0,0,0])
