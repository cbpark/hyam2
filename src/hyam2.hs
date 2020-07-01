module Main where

import           HEP.Kinematics.Variable.M2

import           HEP.Data.LHEF

import           Codec.Compression.GZip     (decompress)
import qualified Data.ByteString.Lazy.Char8 as B
import           Numeric.LinearAlgebra      (Vector, fromList)
import           Pipes
import           Pipes.ByteString           (fromLazy)
import qualified Pipes.Prelude              as P

import           Data.Maybe                 (fromMaybe)
import           System.Environment         (getArgs)

-- import           Debug.Trace

main :: IO ()
main = do
    infile <- head <$> getArgs
    putStrLn $ "-- Parsing " ++ show infile ++ "."

    events <- decompress <$> B.readFile infile
    runEffect $ getLHEFEvent fromLazy events
        >-> P.take 3
        >-> P.map selectP
        >-> P.map mInv2
        >-> P.print

-- | returns four-momenta of b-quarks and leptons, ptmiss
selectP :: Event -> Maybe ([FourMomentum], [FourMomentum], TransverseMomentum)
selectP ev = do
    let topChild = particlesFrom topQuarks (eventEntry ev)
    if null topChild
       then Nothing
       else do
        let bQuarks = concat $ getMomenta ofbQuark    <$> topChild
            leptons = concat $ getMomenta ofIsoLepton <$> topChild
            neus    = concat $ getMomenta ofNeutrino  <$> topChild
        if length bQuarks /= 2 || length leptons /= 2
            then Nothing
            else do let ptmiss = transverseVector (sum neus)
                    return (bQuarks, leptons, ptmiss)
  where
    topQuarks = ParticleType [6]
    getMomenta selF = fmap fourMomentum . filter selF
    ofbQuark    = (`elem` [5, -5])         . idOf
    ofIsoLepton = (`elem` isolatedLeptons) . idOf
    ofNeutrino  = (`elem` neutrinos)       . idOf

mInv2 :: Maybe ([FourMomentum], [FourMomentum], TransverseMomentum)
      -> (Double, Vector Double)
mInv2 = fromMaybe (0, fromList [0,0,0,0]) . mInv2'
  where
    mInv2' ps = do
        (bQuarks, leptons, ptmiss) <- ps
        let [p1, p2] = zipWith (+) bQuarks leptons
        return $ objFunc (InputKinematics p1 p2 ptmiss 0)
            (fromList [px ptmiss, py ptmiss, 0, 0])
