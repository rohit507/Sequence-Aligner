/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/

import scala.math

// Useless Class what which exists so that I don't
//  have to add special rules to the build system.
class ObjectStore {
    val stuff = 0
}

// An Immutable Object filled with all the various settings
//  for sequence alignment and other similar operations. 
class AlignSettings( subf : (Char,Char) => Int , gO : Int, gE : Int, mO : Int,
                       mId : Float, mIg : Int, cB : (Int,Int), kB : (Float,Float)) {
    val costFunc : (Char,Char) => Int = subf
    val gapOpen : Int = gO
    val gapExtend : Int = gE
    val minOverlap : Int = mO
    val minIdentity : Float = mId
    val maxIgnore : Float = mIg
    val minCollisions : Int = cB._1
    val maxCollisions : Int = cB._2
    val collBounds = cB
    val kmerEdge = kB._1
    val kmerCenter = kB._2
    val kmerHeadEdge = kB._1
    val kmerTailEdge = 1.0f - kB._1
    val kmerMidLeadEdge = 0.5f - (kB._2 * 0.5f)
    val kmerMidTailEdge = 0.5f + (kB._2 * 0.5f)
}

// Am immutable object encapsulating the information
//  in a single kmer - sequence pair
class Kmer(iD : Int, kM : String, loC : Float) {
    val int = seqHash(kM) 
    val id = iD
    val km = kM
    val loc = loC

    /* */ // Generates an integer unique to the first 16 bases of the
          //  of the input sequence, useful for a number of operations
    def seqHash( seq : String) : Int = {

        var h : Int = 0

        for (i <- 0 until math.min(16,seq.length)) {
            var c = seq.charAt(i).toUpper
            h = h << 2 
            c match {
                case 'A' => h ^= 0
                case 'C' => h ^= 1
                case 'T' => h ^= 2
                case 'G' => h ^= 3
                case  _  => println("Char \'" + c + "\' is not " +
                                " acceptable input, Found in \"" +
                                seq + "\"")
            }
        }

        return h
    } /* */
}

// A class that stores a sequence and ID
//  in a single convinient blob
class Sequence( iD : Int, sEQ : String) {
    val id = iD
    val seq = sEQ
}

// Static Sequence methods
object Sequence {
    def fromTuple( tup : (Int,String)) : Sequence = {
        tup match {
            case (i,s) => 
                return new Sequence(i,s)
        }
    }
}

// Stores all the relelvent information that can come from a 
//  local alignment between two sequences
class Alignment( sA : Sequence, sB : Sequence, oA : String , oB : String,
                 st : (Int,Int), en : (Int,Int), cor : Int, err : Int) {
    val seqA = sA
    val seqB = sB
    val alignA = oA
    val alignB = oB
    val start = st
    val end = en
    val correct = cor
    val error = err
    val errRatio = (correct.toFloat) / (correct.toFloat + error.toFloat)
    private var overlap : Overlap = null

    def valid(settings : AlignSettings) : Boolean = {
        return (errRatio >= settings.minIdentity) &&
               (alignA.length >= settings.minOverlap) &&
               (((start._1 == 0) && (seqB.seq.length == end._2))  ||
                ((start._2 == 0) && (seqA.seq.length == end._1)))
    }

    def getOverlap() : Overlap = {
        if (overlap == null) {
            overlap = new Overlap(this)
        }
        return overlap
    }
}

// Class to store all the data required to generate a 
//  AMOS bank compatible overlap.   
class Overlap(alg : Alignment) {
    val align = alg
    val adj = 'N'
    val rds = (align.seqA.id,align.seqB.id)
    val scr = 0
    val ahg = align.start._1 - align.start._2
    val bhg = align.seqB.seq.length - align.seqA.seq.length + ahg

    def print() : String = {
        return   "{OVL" +
               "\nadj:" + adj +
               "\nrds:" + rds._1 + "," + rds._2 +
               "\nscr:" + scr +
               "\nahg:" + ahg +
               "\nbhg:" + bhg + 
               "\n}"
    }

    def valid(settings : AlignSettings) : Boolean = {
        return align.valid(settings) &&
               (math.abs(ahg) < settings.maxIgnore) &&
               (math.abs(bhg) < settings.maxIgnore)
    }
}
