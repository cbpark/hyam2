module Main where

import           HEP.Data.LHEF              (getLHEFEvent)

import           Codec.Compression.GZip     (decompress)
import qualified Data.ByteString.Lazy.Char8 as B
import           Pipes
import           Pipes.ByteString           (fromLazy)
import qualified Pipes.Prelude              as P

import           System.Environment         (getArgs)

main :: IO ()
main = do
    infile <- head <$> getArgs
    putStrLn $ "-- Parsing " ++ show infile ++ "."

    events <- decompress <$> B.readFile infile
    runEffect $ getLHEFEvent fromLazy events
        >-> P.take 3
        >-> P.print
