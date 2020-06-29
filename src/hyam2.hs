module Main where

import           HEP.Kinematics.Variable.M2

import           HEP.Data.LHEF

import           Codec.Compression.GZip     (decompress)
import qualified Data.ByteString.Lazy.Char8 as B
import           Pipes
import           Pipes.ByteString           (fromLazy)
import qualified Pipes.Prelude              as P

import           Data.Maybe                 (fromMaybe)
import           System.Environment         (getArgs)

main :: IO ()
main = do
    infile <- head <$> getArgs
    putStrLn $ "-- Parsing " ++ show infile ++ "."

    events <- decompress <$> B.readFile infile
    runEffect $ getLHEFEvent fromLazy events
        >-> P.take 3
        >-> P.map selectP
        -- >-> P.map mTtot
        >-> P.map mInv2
        >-> P.print

selectP :: Event -> Maybe ([FourMomentum], TransverseMomentum)
selectP ev = do
    let topChild = particlesFrom topQuarks (eventEntry ev)
    if null topChild
        then Nothing
        else do
            let pVis = momentumSum' . filter (not . isNeutrino) <$> topChild
                pNu  = sum $ momentumSum' . filter isNeutrino   <$> topChild
                ptmiss = transverseVector pNu
            return (pVis, ptmiss)
  where
    topQuarks = ParticleType [6]
    isNeutrino = (`elem` neutrinos) . idOf
    momentumSum' = momentumSum . fmap fourMomentum

mInv2 :: Maybe ([FourMomentum], TransverseMomentum) -> Double
mInv2 = fromMaybe 0 . mInv2'
  where
    mInv2' ps = do
        (pVis, ptmiss) <- ps
        if length pVis /= 2
            then Nothing
            else do let [v1, v2] = pVis
                    return $ objFunc (InputKinematics v1 v2 ptmiss 0)
                        (px ptmiss, py ptmiss, 0, 0)

-- mTtot :: Maybe ([FourMomentum], TransverseMomentum) -> Double
-- mTtot = fromMaybe 0 . mTtot'
--   where
--     mTtot' ps = do (pVis, ptmiss) <- ps
--                    return $ transverseMass pVis (promote2TV ptmiss 0)
