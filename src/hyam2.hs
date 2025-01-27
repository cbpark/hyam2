module Main where

import           HEP.Kinematics.Variable.M2

import           HEP.Data.LHEF

import           Codec.Compression.GZip     (decompress)
import qualified Data.ByteString.Lazy.Char8 as B
import           Pipes
import           Pipes.ByteString           (fromLazy)
import qualified Pipes.Prelude              as P

-- import           Data.Maybe                 (fromMaybe)
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
      -> Maybe M2Solution
mInv2 = -- fromMaybe (Left INVALID_ARGS) . mInv2'
    mInv2'
  where
    mInv2' ps = do
        (bQuarks, leptons, ptmiss) <- ps
        -- m2SQP (mkInput bQuarks leptons ptmiss 0)
        -- m2XXSQP (mkInput bQuarks leptons ptmiss 0)
        -- m2CXSQP (mkInput bQuarks leptons ptmiss 0)
        -- m2XCSQP (mkInput bQuarks leptons ptmiss 0)
        m2CCSQP (mkInput bQuarks leptons ptmiss 0)
        -- m2XXAugLag (mkInput bQuarks leptons ptmiss 0)
        -- m2CXAugLag (mkInput bQuarks leptons ptmiss 0)
        -- m2XCAugLag (mkInput bQuarks leptons ptmiss 0)
        -- m2CCAugLag (mkInput bQuarks leptons ptmiss 0)
