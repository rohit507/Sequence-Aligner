/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/

import java.io._

import scala.collection._
import scala.actors._
import scala.actors.Futures._

object BioLibs {

    type Alignment = (String,String,Int,Int,Int,Int)
    type Overlap = (Char,(Int,Int),Int,Int,Int)

    // Read in a fasta file dispatching each sequence as it's
    //  read through the input function. 
    def readSeq( file : String, act : (Int,String) => _) : Int = {
        var i = 1
        var s = ""
        var in = new BufferedReader(new FileReader(file))
        var line = in.readLine()

        if (!line.startsWith(">"))
            throw new Exception("Invalid Sequence File: " + file)

        line = in.readLine()
        while (line != null) {
            if (line.startsWith(">")) {
                act(i,s.toUpperCase)
                i += 1
                s = ""
            } else {
                s += line
            }
            line = in.readLine()
        }
    
        act(i,s.toUpperCase)

        return i
    }

    // Generates an integer unique to the first 16 bases of the
    //  of the input sequence, useful for a number of operations
    def seqHash( seq : String) : Int = {

        var h : Int = 0

        for (i <- 0 to 16) {
            var c = Character.toUpperCase(seq.charAt(i))

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
    }

    // generates a mutable hashSet of strings 
    def generateKmerSet( k : Int )( s : String ) 
                        : (Int,Set[String]) = {
        var o = new mutable.HashSet[String]()

        for( i <- 0 to (s.length - k)) {
            o += s.substring(i,i+k)
        }
        return (k,o)
    }

    def generateKmerTable( file : String, size : Int, debug : Boolean) : KmerTable = {

        // It is faster not to use futures for the kmer splitting

        var kmerSets = List[(Int,String,Set[String])]()
        var kmerTable = new KmerTable()
        var gKS = generateKmerSet(size)_

        readSeq(file,(id : Int, seq : String) => {
           if (debug) { println("Read Seq : " + id) }
           if (debug) { println("Start Processing Seq : " + id)}
           var s = gKS(seq)
           if (debug) { println("Generated Kmers for Seq : " + id +
                                " on thread " + Thread.currentThread.getName ) }
            val f = (id,seq,s._2)
           if (debug) { println("Sent Seq : " + id) }
           kmerSets = kmerSets ::: List(f)
        })

        kmerSets.foreach { f =>
            f match {
                case (id,seq,set) => kmerTable.addKmerSet(id,seq,set)
                                if (debug) { println("Added Kmer Set : " + id) }
                case _ => if (debug) { println("error in future return") }
            }
        }
        
        return kmerTable
    }

    def generateFutureKmerTable( file : String, size : Int,
                                 debug : Boolean) : KmerTable = {

        var kmerSets = List[Future[(Int,String,Set[String])]]()
        var kmerTable = new KmerTable()
        var gKS = generateKmerSet(size)_

        readSeq(file,(id : Int, seq : String) => {
           if (debug) { println("Read Seq : " + id) }
           val f = future {
                if (debug) { println("Start Processing Seq : " + id)}
                var s = gKS(seq)
                if (debug) { println("Generated Kmers for Seq : " + id +
                                " on thread " + Thread.currentThread.getName ) }
                (id,seq,s._2)
           }
           if (debug) { println("Sent Seq : " + id) }
           kmerSets = kmerSets ::: List(f)
        })

        kmerSets.foreach { f =>
            f.apply() match {
                case (id,seq,set) => kmerTable.addKmerSet(id,seq,set)
                                if (debug) { println("Added Kmer Set : " + id) }
                case _ => if (debug) { println("error in future return") }
            }
        }
        
        return kmerTable
    } 
    
    def generateKmerCover( file : String, size : Int, debug : Boolean) : (Int,Float,KmerTable) = {
        val possible = math.pow(4,size) 
        val kmers = generateKmerTable(file,size,false)
        val uniques = kmers.uniqueKmers()
        val ratio = (uniques.toFloat) / (possible.toFloat)
        if (debug) { println("There are " + uniques + " unique " + size +
                "mers in the data out of a possible " + possible +
                " for a ratio of " + ratio ) }
        return (uniques,ratio,kmers)
    }

    def readHoxd( file : String ) : (Char,Char) => Int = {

            var costs = new mutable.HashMap[String,Int]()
            var HOXD = new BufferedReader(new FileReader(file))
            HOXD.readLine(); //Remove the title line
            var col = HOXD.readLine().split(",")
            var line = HOXD.readLine();

            while((line != null) && (line != "")) {

                var row = line.split(",")
                for(i <- 1 until row.length) {
                    costs.put(row(0).trim().toUpperCase()+
                              col(i).trim().toUpperCase(),
                                Integer.parseInt(row(i)))
                }

                line = HOXD.readLine()
            }

            return (a : Char, b : Char) => {
                costs.apply(a.toString.toUpperCase + b.toString.toUpperCase)
            }
    }

    def simpleMatch( mat : Int, miss : Int) : (Char,Char) => Int = {
        return (a : Char, b : Char) => { if (a == b) mat else miss }
    }

    def generateLocalAlignment(A : String, B : String, settings : AlignSettings) 
                            : Alignment = {
        
	    var M : Array[Array[Int]] = Array.ofDim(A.length()+1,B.length()+1)
        var X : Array[Array[Int]] = Array.ofDim(A.length()+1,B.length()+1)
	    var Y : Array[Array[Int]] = Array.ofDim(A.length()+1,B.length()+1)

	    for (i <- 0 until A.length()){//Fills in zero row
	        M(i)(0) = 0
	        X(i)(0) = 0
            Y(i)(0) = settings.gapOpen + i * settings.gapExtend;
	    }

	    for (i <- 0 until B.length()){//Fills in zero columns
            M(0)(i) = 0
	        X(0)(i) = settings.gapOpen + i * settings.gapExtend;
            Y(0)(i) = 0
	    }

        var max = 0
        var maxLoc = (0,0)
        var t = 0

	    for (i <- 1 to A.length()){
	        for (j <- 1 to B.length()){
	    	//Based off of history, I created the matrix
    	        M(i)(j) = settings.costFunc(A.charAt(i-1),B.charAt(j-1)) + 
                           math.max(math.max(M(i-1)(j-1),Y(i-1)(j-1)),
                                    math.max(X(i-1)(j-1),0))

                X(i)(j) = settings.gapExtend + 
                           math.max(math.max(M(i)(j-1) + settings.gapOpen,
                                             Y(i)(j-1) + settings.gapOpen),
                                    math.max(X(i)(j-1),0))
           
                Y(i)(j) = settings.gapExtend + 
                           math.max(math.max(M(i-1)(j) + settings.gapOpen, Y(i-1)(j)),
                                    math.max(X(i-1)(j) + settings.gapOpen, 0))

                t = math.max(M(i)(j),math.max(X(i)(j),Y(i)(j)))

                if (t > max) {
                    max = t
                    maxLoc = (i,j)
                }
            }
	    }

        var i = maxLoc._1 
        var j = maxLoc._2
        var xSeq = "" 
        var ySeq = ""
        var c = 0
        var e = 0
        var pa = ' '
        var pb = ' '

        max = math.max(M(i)(j),math.max(X(i)(j),Y(i)(j)))

        do {
            if (M(i)(j) == max) {
                pa = A.charAt(i-1)
                pb = B.charAt(j-1)
                i -= 1
                j -= 1
            } else if (X(i)(j) == max) {
                pa = A.charAt(i-1)
                pb = '-'
                j -= 1
            } else if (Y(i)(j) == max) {
                pa = '-'
                pb = B.charAt(j-1)
                i -= 1
            }

            c += 1
            if(pa != pb) {
                e += 1
            }      
 
            xSeq = pa + xSeq
            ySeq = pb + ySeq
      
            max = math.max(M(i)(j),math.max(X(i)(j),Y(i)(j)))
        } while(max > 0)

        return (xSeq,ySeq,i,j,c,e)
    }

    def generateOverlap( seqA : (Int,String) , seqB : (Int,String),
                         align : Alignment, settings : AlignSettings ) 
                         : Overlap = {
        var adj = 'N'
        var rds = (seqA._1,seqB._1)
        var scr = 0
        var ahg = 0
        var bhg = 0
        if (align._3 > align._4) {
            ahg = - align._3 + align._4
            bhg = - seqB._2.length + align._4 + align._2.length - 
                    align._2.count('-' == _) + seqA._2.length - align._3 -
                    align._1.length + align._1.count('-' == _)
        } else {
            bhg = align._4 - align._3
            ahg = seqA._2.length - seqB._2.length + bhg
        }
        return (adj,rds,scr,ahg,bhg)
    }

    def printOverlap(overlap : Overlap) : String = {
        return   "{OVL" +
               "\nadj:" + overlap._1 +
               "\nrds:" + overlap._2._1 + "," + overlap._2._2 +
               "\nscr:" + overlap._3 +
               "\nahg:" + overlap._4 +
               "\nbhg:" + overlap._5 + 
               "\n}"
    }
     
}

class AlignSettings( subf : (Char,Char) => Int , gO : Int, gE : Int, mO : Int,
                       mId : Float, mIg : Int ) {
    val costFunc : (Char,Char) => Int = subf
    val gapOpen : Int = gO
    val gapExtend : Int = gE
    val minOverlap : Int = mO
    val minIdentity = mId
    val maxIgnore = mIg
}


