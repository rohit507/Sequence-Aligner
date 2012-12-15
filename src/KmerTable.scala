/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/

import scala.collection._
import scala.collection.mutable.OpenHashMap
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.UnrolledBuffer
import scala.math

import gnu.trove.map.hash.TIntObjectHashMap

import java.lang.Runtime

class KmerTable {

    var debug = false
    private def printdb(s : String) { if(debug){ println(s)}}

    val rt = Runtime.getRuntime()

    // Hash of kmer strings to the various pieces of kmer data taken
    //  from input strings.
    val KmerData = new TIntObjectHashMap[ArrayBuffer[Kmer]]()

    // Hash of dovetail pairs to the number of kmers they share (or more
    //  accurately the number of times the KmerTable has been asked to pair them)
    val PairData = new TIntObjectHashMap[Int]()

    // Hash from lead pair Sequence ID to all the IDs of sequences that 
    //  it might dovetail into.
    val DispatchData = new TIntObjectHashMap[ArrayBuffer[Int]]()

    // Hash of sequece IDs to their respective sequence objects
    val SequenceData = new TIntObjectHashMap[Sequence]()

    // Adds Kmer and Sequence Data to their respective 
    //  storage collections.
    def addKmerSet( seq : Sequence, kmerSet : ArrayBuffer[Kmer] ) {

        SequenceData.put(seq.id, seq)

        for (s <- kmerSet) {

            if (! KmerData.contains(s.int)) {
                KmerData.put(s.int, new mutable.ArrayBuffer[Kmer]())
            }
      
            KmerData.get(s.int).append(s)
        }
    }

    // Adds a kmer pair to the pair counting collection in 
    //  first (dovetail) second order.
    def addKmerPair(a : Kmer, b : Kmer) {
        var fst : Kmer = null
        var snd : Kmer = null

        if (a.id == b.id) {
            return
        }

        if (a.loc > b.loc) {
            fst = a
            snd = b
        } else {
            fst = b
            snd = a
        }

        val idTup = (fst.id << 16) ^ snd.id

        if (! PairData.contains(idTup)) {
            PairData.put(idTup,0)
        }

        PairData.put(idTup, PairData.get(idTup) + 1)
    }

    // Goes through the kmer collections marking off pairs of
    //  sequences that are likely to produce useful dovetail
    //  alignments 
    def calcPairData(s : AlignSettings) {

        PairData.clear()

        val stStore = ArrayBuffer[Kmer]()
        val mdStore = ArrayBuffer[Kmer]()
        val enStore = ArrayBuffer[Kmer]()

        var i = 0
        val iter = KmerData.iterator()

        while (iter.hasNext) {

            iter.advance()
            var kmerArr = iter.value()

            stStore.clear()
            mdStore.clear()
            enStore.clear()

            for (k <- kmerArr) {
                if(k.loc <= s.kmerHeadEdge) {
                    stStore.append(k)
                }
                if((s.kmerMidLeadEdge <= k.loc) && 
                   (k.loc <= s.kmerMidTailEdge)) {
                    mdStore.append(k)
                }
                if(s.kmerTailEdge <= k.loc) {
                    enStore.append(k)
                }
            }

            for (st <- stStore) {
                for (md <- mdStore) {
                    addKmerPair(st,md)
                }
            }

            for (en <- enStore) {
                for (md <- mdStore) {
                    addKmerPair(en,md)
                }
            }

            i += 1
            if( (i % 1000) == 0) {
                printdb("Analysed kmer pair " + i + " of " 
                    + KmerData.size + ".")
            }
         

            if( (i % 15000) == 0) {
                printdb("There is " + rt.freeMemory() +
                    " free memory, out of " + rt.totalMemory() +
                    " total memory." )
                printdb("Garbage Collecting. ")
                System.gc()
                printdb("Finished Garbage Collecting.")
                printdb("There is now " + rt.freeMemory() +
                    " free memory, out of " + rt.totalMemory() +
                    " total memory." )
            }
        }
    }

    // Takes the set of kmer pairs and assembles a head seq
    //  to set of possible dovetails mappings. Count range 
    //  should be very small, possibly even 1 or 2, through
    //  to very large, 40 or 50.
    def calcDispatchData(s : AlignSettings) {

        DispatchData.clear()

        var i = 0
        var a = 0
        var b = 0
        val iter = PairData.iterator()

        while (iter.hasNext) {

            iter.advance()
            var int = iter.key()
            var count = iter.value() 
            a = int >> 16
            b = (int << 16) >> 16
            if((s.minCollisions <= count) && 
               (count <= s.maxCollisions)) {
                if(! DispatchData.contains(a)) {
                    DispatchData.put(a,new ArrayBuffer[Int]())
                }
                DispatchData.get(a).append(b)
            }

            i += 1
            if( (i % 254437) == 0) {
                printdb("Analysed dispatch set " + i + " of " +
                    PairData.size + ".")
            }
        }

        printdb("Created " + DispatchData.size() + " dispatch blocks.")
    }

    def uniqueKmers() : Int = {
        return KmerData.size
    } 

    def uniqueSeqs() : Int = {
        return SequenceData.size
    }

    // Generate a histogram relating kmer collision counts
    //  and the number of seperate unique kmers with that
    //  many collisions.
    def kmerCollisionHistogram() : Map[Int,Int] = {

        var hist = new OpenHashMap[Int,Int]()

        val iter = KmerData.iterator()

        while (iter.hasNext) {

            iter.advance()
            var i = iter.key()
            var s = iter.value()    
        
            if (! hist.contains(s.size)) {
                hist.put(s.size,0)
            }

            hist.put(s.size,hist.apply(s.size) + 1)

        }

        return hist
    }

    // Send a set of sequence paris to be aligned one by one to
    //  act function. 
    def dispatchCollisions( set : AlignSettings , act : (Sequence,Sequence) => _) {
        calcPairData(set)
        calcDispatchData(set)

        val iter = DispatchData.iterator()

        while (iter.hasNext) {

            iter.advance()
            var i = iter.key()
            var a = iter.value()  
  
            for (j <- a) {
                act(SequenceData.get(i),SequenceData.get(j))
            }
        }
        System.gc()
    }

    // Send a sequence and stack of sequences to aling it with. 
    //  acting this way can save a bunch of memory among other things
    def dispatchCollisionBlocks( s : AlignSettings , act : 
                                (Int,Sequence,Seq[Sequence]) => _) {
        calcPairData(s)
        calcDispatchData(s)

        val iter = DispatchData.iterator()

        while (iter.hasNext) {

            iter.advance()
            var i = iter.key()
            var a = iter.value()  

            var set = new ArrayBuffer[Sequence](a.length)
            var max = 0

            for (j <- a) {                
                val t = SequenceData.get(j)
                set += t
                max = math.max(max,t.seq.length)
            }

            if( set.size >= 1) {
                act(max,SequenceData.get(i),set)
            }
        }
        System.gc()
    }

} 
    
