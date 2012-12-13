/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/

import java.io._

import scala.collection._
import scala.actors._
import scala.actors.Futures._
import scala.collection.immutable.StringOps

object BioLibs {

    // Read in a fasta file dispatching each sequence as it's
    //  read through the input function. 
    def readSeq( file : String, act : (Sequence) => _) : Int = {
        var i = 1
        var s = ""
        var in = new BufferedReader(new FileReader(file))
        var line = in.readLine()

        if (!line.startsWith(">"))
            throw new Exception("Invalid Sequence File: " + file)

        line = in.readLine()
        while (line != null) {
            if (line.startsWith(">")) {
                act(new Sequence(i,s.toUpperCase))
                i += 1
                s = ""
            } else {
                s += line
            }
            line = in.readLine()
        }
    
        act(new Sequence(i,s.toUpperCase))

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
    def generateKmerSet( k : Int, s : String ) : (Set[String]) = {
        var o = new mutable.HashSet[String]()

        for( i <- 0 to (s.length - k)) {
            o += s.substring(i,i+k)
        }

        return o
    }

    def readHOXD( file : String ) : (Char,Char) => Int = {

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

    def generateLocalAlignment(seqA : Sequence, seqB : Sequence,
                         settings : AlignSettings) : Alignment = {
        
        val A = seqA.seq
        val B = seqB.seq

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

        var opt = maxLoc
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

            if(pa != pb) {
                e += 1
            } else {
                c += 1
            } 
 
            xSeq = pa + xSeq
            ySeq = pb + ySeq
      
            max = math.max(M(i)(j),math.max(X(i)(j),Y(i)(j)))
        } while(max > 0)

        return new Alignment(seqA,seqB,xSeq,ySeq,(i,j),opt,c,e)
    }
     
}



