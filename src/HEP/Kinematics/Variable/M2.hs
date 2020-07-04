{-# LANGUAGE RecordWildCards #-}

module HEP.Kinematics.Variable.M2 where

import HEP.Kinematics
import HEP.Kinematics.Vector.LorentzVector (setXYZT)

import Numeric.LinearAlgebra               (Vector, fromList, toList)
import Numeric.NLOPT

import Debug.Trace

-- | All the components are rescaled by '_scale'.
data InputKinematics = InputKinematics
                       { _p1     :: FourMomentum  -- ^ p1 = a1 + b1
                       , _p2     :: FourMomentum  -- ^ p2 = a2 + b2
                       , _q1     :: FourMomentum  -- ^ q1 = b1
                       , _q2     :: FourMomentum  -- ^ q2 = b2
                       , _ptmiss :: TransverseMomentum
                         -- | squared mass of the invisible particle
                       , _mInvSq :: !Double
                       , _scale  :: !Double
                       }

-- | creates the input for the process of
--
--   A1 + A2 --> a1 B1 + a2 B2 --> a1 b1 C1 + a2 b2 C2.
--
mkInput :: [FourMomentum]      -- ^ [a1, a2]
        -> [FourMomentum]      -- ^ [b1, b2]
        -> TransverseMomentum  -- ^ pT(miss)
        -> Double              -- ^ M_{invisible}
        -> Maybe InputKinematics
mkInput as bs ptmiss mInv
    | length as /= 2 || length bs /= 2 = Nothing
    | otherwise                        = do
          let ps = zipWith (+) as bs
              p1 = head ps
              p2 = ps !! 1

              m1 = mass p1
              m2 = mass p2
              mInvSq = mInv * mInv
              scaleSq = getScale p1 + getScale p2 + getScale ptmiss
                        + m1 * m1 + m2 * m2 + 2 * mInvSq
          if scaleSq <= 0  -- why?
              then Nothing
              else do let scale = sqrt scaleSq
                          p1'     = p1        ^/ scale
                          p2'     = p2        ^/ scale
                          q1'     = head bs   ^/ scale
                          q2'     = (bs !! 1) ^/ scale
                          ptmiss' = ptmiss    ^/ scale
                          mInv'   = mInv       / scale
                      return $ InputKinematics { _p1     = p1'
                                               , _p2     = p2'
                                               , _q1     = q1'
                                               , _q2     = q2'
                                               , _ptmiss = ptmiss'
                                               , _mInvSq = mInv' * mInv'
                                               , _scale  = scale
                                               }

getScale :: HasFourMomentum a => a -> Double
getScale v = let ptv = pt v in ptv

data M2Solution = M2Solution { _M2    :: Double
                             , _k1sol :: FourMomentum
                             , _k2sol :: FourMomentum
                             } deriving Show

type MultivarFunc = Vector Double -> (Double, Vector Double)

m2SQP :: [InputKinematics -> MultivarFunc]  -- ^ constraint functions
      -> Maybe InputKinematics
      -> Maybe M2Solution
m2SQP _   Nothing                         = Nothing
m2SQP cfs (Just inp@InputKinematics {..}) =
    let -- objective function with gradient
        objfD = m2ObjF inp

        eps1 = 1e-2 -- 1e-3 * _scale
        -- constraint function with gradient
        constraints' = ($ inp) <$> cfs
        constraints  = (\cf -> EqualityConstraint (Scalar cf) eps1)
                       <$> constraints'

        eps2 = eps1 * 1e-2 -- 1e-6 * _scale
        -- stop = ObjectiveRelativeTolerance 1e-9 :| []
        stop = ObjectiveAbsoluteTolerance eps2
               :| [ObjectiveRelativeTolerance eps2, MaximumEvaluations 100]
        algorithm = SLSQP objfD [] [] constraints
        problem = LocalProblem 4 stop algorithm

        -- initial guess
        x0 = fromList [0.5 * px _ptmiss, 0.5 * py _ptmiss, 0, 0]
        sol = minimizeLocal problem x0
    in getM2Solution inp sol

m2XXSQP, m2CXSQP, m2XCSQP, m2CCSQP :: Maybe InputKinematics -> Maybe M2Solution
m2XXSQP = m2SQP []
m2CXSQP = m2SQP [constraintA]
m2XCSQP = m2SQP [constraintB]
m2CCSQP = m2SQP [constraintA, constraintB]

m2AugLag :: Maybe InputKinematics -> Maybe M2Solution
m2AugLag Nothing                         = Nothing
m2AugLag (Just inp@InputKinematics {..}) =
    let objfD = m2ObjF inp
        objf  = getResultOnly objfD

        eps1 = 1e-2
        c1fD = constraintA inp
        c1f  = getResultOnly c1fD
        c1 = EqualityConstraint (Scalar c1f) eps1
        c2fD = constraintA inp
        c2f  = getResultOnly c2fD
        c2 = EqualityConstraint (Scalar c2f) eps1

        eps2 = eps1 * 1e-2
        stop = ObjectiveAbsoluteTolerance eps2
               :| [ObjectiveRelativeTolerance eps2, MaximumEvaluations 5000]
        algorithm = NELDERMEAD objf [] Nothing

        subproblem = LocalProblem 4 stop algorithm
        problem = AugLagProblem [c1, c2] [] (AUGLAG_EQ_LOCAL subproblem)

        -- initial guess
        x0 = fromList [0.5 * px _ptmiss, 0.5 * py _ptmiss, 0, 0]
        sol = minimizeAugLag problem x0
    in getM2Solution inp sol

getResultOnly :: MultivarFunc -> Vector Double -> Double
getResultOnly fD ks = let (result, _) = fD ks in result

getM2Solution :: InputKinematics -> Either Result Solution -> Maybe M2Solution
getM2Solution inp@InputKinematics {..} sol =
    case sol of
        Left _                          -> Nothing
        sol0@(Right (Solution m2 ks _)) -> do
            let Invisibles k1 k2 = mkInvisibles inp (vecToVars ks)
            traceM ("sol = " ++ show sol0)
            return $ M2Solution { _M2    = m2  * _scale
                                , _k1sol = k1 ^* _scale
                                , _k2sol = k2 ^* _scale }

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
m2ObjF :: InputKinematics -> MultivarFunc
m2ObjF inp@InputKinematics {..} ks =
    if m1 < m2 then (m2, grad2) else (m1, grad1)
  where
    invs@(Invisibles k1 k2) = mkInvisibles inp (vecToVars ks)
    m1 = invariantMass [_p1, k1]
    m2 = invariantMass [_p2, k2]
    (grad1, grad2, _) = m2Grad inp invs _p1 _p2 ks

m2Grad :: InputKinematics
       -> Invisibles
       -> FourMomentum  -- ^ p1 (or q1)
       -> FourMomentum  -- ^ p2 (or q2)
       -> Vector Double
       -> (Vector Double, Vector Double, Double)
m2Grad InputKinematics {..} (Invisibles k1 k2) p1 p2 ks = (d1, d2, m1 - m2)
  where
    (k1x, k1y, k1z, k2z) = vecToVars ks

    (e1, p1x, p1y, p1z) = epxpypz p1
    m1  = invariantMass [p1, k1]
    m1' = safeDivisor m1
    r1 = e1 / safeDivisor (energy k1)
    d1 = fromList $ ( / m1') <$> [ r1 * k1x - p1x
                                 , r1 * k1y - p1y
                                 , r1 * k1z - p1z
                                 , 0 ]

    (e2, p2x, p2y, p2z) = epxpypz p2
    m2  = invariantMass [p2, k2]
    m2' = safeDivisor m2
    r2 = e2 / safeDivisor (energy k2)
    d2 = fromList $ ( / m2') <$> [ r2 * (k1x - px _ptmiss) + p2x
                                 , r2 * (k1y - py _ptmiss) + p2y
                                 , 0
                                 , r2 * k2z - p2z ]

constraintF :: FourMomentum
            -> FourMomentum
            -> InputKinematics
            -> Vector Double
            -> (Double, Vector Double)
constraintF p1 p2 inp ks = (deltaM, grad1 - grad2)
  where
    invs = mkInvisibles inp (vecToVars ks)
    (grad1, grad2, deltaM) = m2Grad inp invs p1 p2 ks

constraintA, constraintB :: InputKinematics -> MultivarFunc
constraintA inp@InputKinematics {..} = constraintF _p1 _p2 inp
constraintB inp@InputKinematics {..} = constraintF _q1 _q2 inp

-- | unsafe transformation, but it would be fine.
vecToVars :: Vector Double -> Variables
vecToVars ks = let [k1x, k1y, k1z, k2z] = toList ks in (k1x, k1y, k1z, k2z)

safeDivisor :: Double -> Double
safeDivisor x | x >= 0    = max eps x
              | otherwise = min (-eps) x
  where eps = 1.0e-8
