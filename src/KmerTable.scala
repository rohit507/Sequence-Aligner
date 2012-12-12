/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/

import BioLibs._

import scala.collection._

class KmerTable {
    var KmerData    : Map[String,mutable.HashSet[Int]] = 
                new mutable.HashMap[String,mutable.HashSet[Int]]()
    var SequenceData : Map[Int,String] = new mutable.HashMap[Int,String]()

    def addKmerSet( id : Int, seq : String, kmerSet : Set[String] ) {

        SequenceData += ((id , seq))

        for (s <- kmerSet) {

            if (! KmerData.contains(s)) {
                KmerData += ((s, new mutable.HashSet[Int]()))
            }
      
            KmerData.apply(s).add(id)
        }
    }

    def uniqueKmers() : Int = {
        return KmerData.size
    }

    def kmerCollisionHistogram() : Map[Int,Int] = {
        var hist = new mutable.HashMap[Int,Int]()

        for ((i,s) <- KmerData) {
        
            if (! hist.contains(s.size)) {
                hist += ((s.size,0))
            }

            hist += ((s.size,hist.apply(s.size) + 1))

        }

        return hist
    }
} 
    
