/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/

class ObjectStore {
    val stuff = 0
}

class AlignSettings( subf : (Char,Char) => Int , gO : Int, gE : Int, mO : Int,
                       mId : Float, mIg : Int, mC : Int) {
    val costFunc : (Char,Char) => Int = subf
    val gapOpen : Int = gO
    val gapExtend : Int = gE
    val minOverlap : Int = mO
    val minIdentity : Float = mId
    val maxIgnore : Float = mIg
    val minCollisions : Int = mC
}

class Sequence( iD : Int, sEQ : String) {
    val id = iD
    val seq = sEQ
}

object Sequence {
    def fromTuple( tup : (Int,String)) : Sequence = {
        tup match {
            case (i,s) => 
                return new Sequence(i,s)
        }
    }
}

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

    def valid(settings : AlignSettings) : Boolean = {
        return (errRatio > settings.minIdentity) &&
               (alignA.length > settings.minOverlap) &&
               (((start._1 == 0) && (seqB.seq.length == end._2))  ||
                ((start._2 == 0) && (seqA.seq.length == end._1)))
    }

    def getOverlap() : Overlap = {
        return new Overlap(this)
    }
}
    
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
