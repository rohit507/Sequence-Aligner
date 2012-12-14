/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/

import scala.collection._
import scala.math._
import scala.math

class KmerTable {

    var KmerData = new mutable.HashMap[String,mutable.HashSet[Int]]()
    var PairCountData = new mutable.HashMap[(Int,Int),Int]]()
    var SequenceData = new mutable.HashMap[Int,Sequence]()

    def addPairCount( idA : Int , idB : Int ) {
        val n = math.min(idA,idB)
        val m = math.max(idA,idB)
        if(n != m) {
            if (! PairCountData.contains((n,m))) {
                PairCountData += (((n,m), 0))
            }
            PairCountData.put((n,m),PairCountData.apply((n,m)) + 1)
        }
    }

    def calculatePairCounts() {
        for((_,s) <- KmerData) {
            val ar = s.toArray
            for( i <- 1 until ar.length) {
                for( j <- (i + 1) until ar.length ) {
                    addPairCount(ar(i),ar(j))
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

    def dispatchCollisions( collBounds : (Int,Int) , act : (Sequence,Sequence) => _) {
        calculatePairCounts()
        for ((i,s) <- PairCountData) {
            for ((j,count) <- s) {
                if ((count >= collBounds._1) && (count >= collBounds._2)) {
                    act(SequenceData.apply(i),SequenceData.apply(j))
                }
            }
        }
    }

    def dispatchCollisionBlocks( collBounds : (Int,Int) , act : 
                                (Int,Sequence,Seq[Sequence]) => _) {
        calculatePairCounts()
        for ((i,s) <- PairCountData) {
            var set = mutable.Queue[Sequence]()
            var max = 0
            for ((j,count) <- s) {                
                if ((count >= collBounds._1) && (count >= collBounds._2)) {
                    val t = SequenceData.apply(j)
                    set += t
                    max = math.max(max,t.seq.length)
                }
            }
            if( set.size >= 1) {
                act(max,SequenceData.apply(i),set)
            }
        }
    }

} 
    
