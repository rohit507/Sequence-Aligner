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
    var SequenceData = new mutable.HashMap[Int,String]()

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

    def dispatchCollisions( bounds : (Int,Int) , act : ((Int,String),(Int,String)) => _) {
        var collisionSet = new mutable.HashSet[(Int,Int)]()
        for ((k,ss) <- KmerData) {
            var arr = ss.toArray
            if ((arr.size >= bounds._1) && (arr.size <= bounds._2)) {
                for (i <- 0 until arr.size) {
                    for( j <- (i + 1) until arr.size) {
                        var n = min(arr(i),arr(j))
                        var x = max(arr(i),arr(j))
                        if(!collisionSet.contains((n,x))) {
                            collisionSet += ((n,x))
                            act((n,SequenceData.apply(n)),(x,SequenceData(x)))
                        }
                    }
                }
            }
        }
    }

} 
    
