/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/

import scala.collection._
import scala.math._
import scala.math

class KmerTable {

    var KmerData = new mutable.HashMap[String,mutable.HashSet[Int]]()
    var AlignData = new mutable.HashMap[Int,mutable.HashMap[Int,Int]]()
    var SequenceData = new mutable.HashMap[Int,Sequence]()

    def addPairAlignment( idA : Int , idB : Int ) {
        val n = math.min(idA,idB)
        val m = math.max(idA,idB)
        if(n != m) {
            if (! AlignData.contains(n)) {
                AlignData += ((n, new mutable.HashMap[Int,Int]()))
            }

            if (! AlignData.apply(n).contains(m)) {
                AlignData.apply(n).put(m,0)
            }

            AlignData.apply(n).put(m,AlignData.apply(n).apply(m) + 1)
        }
    }

    def calculatePairAlignments() {
        for((_,s) <- KmerData) {
            val ar = s.toArray
            for( i <- 1 until ar.length) {
                for( j <- (i + 1) until ar.length ) {
                    addPairAlignment(ar(i),ar(j))
                }
            }
        }
    }

    def addKmerSet( seq : Sequence, kmerSet : Set[String] ) {

        SequenceData += ((seq.id, seq))

        for (s <- kmerSet) {

            if (! KmerData.contains(s)) {
                KmerData += ((s, new mutable.HashSet[Int]()))
            }
      
            KmerData.apply(s).add(seq.id)
        }
    }

    def uniqueKmers() : Int = {
        return KmerData.size
    }

    def uniqueSeqs() : Int = {
        return SequenceData.size
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

    def dispatchCollisions( minCollision : Int , act : (Sequence,Sequence) => _) {
        calculatePairAlignments()
        for ((i,s) <- AlignData) {
            for ((j,count) <- s) {
                if (count > minCollision) {
                    act(SequenceData.apply(i),SequenceData.apply(j))
                }
            }
        }
    }

    def dispatchCollisionBlocks( minCollision : Int , act : (Sequence,Set[Sequence]) => _) {
        calculatePairAlignments()
        for ((i,s) <- AlignData) {
            var set = mutable.HashSet[Sequence]()
            for ((j,count) <- s) {                
                if (count > minCollision) {
                    set += SequenceData.apply(j)
                }
            }
            if( set.size >= 1) {
                act(SequenceData.apply(i),set)
            }
        }
    }

} 
    
