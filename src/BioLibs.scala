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
    def generateKmerSet( k : Int, s : Sequence ) : (ArrayBuffer[Kmer]) = {
        var o = new mutable.ArrayBuffer[Kmer](s.seq.length - k)
        var d = (s.seq.length - k).toFloat
        for( i <- 0 to (s.seq.length - k)) {
            o += new Kmer(s.id,s.seq.substring(i,i+k),(i.toFloat)/d)
        }
        return o
    }

    def readHOXD( file : String ) : (Char,Char) => Int = {

            var costs : Array[Array[Int]] = Array.ofDim(4,4)
            var HOXD = new BufferedReader(new FileReader(file))
            HOXD.readLine(); //Remove the title line
            var col = HOXD.readLine().split(",")
            var line = HOXD.readLine();
            var A = 0
            var B = 0
            while((line != null) && (line != "")) {

                var row = line.split(",")
                for(i <- 1 until row.length) {

                    row(0).trim().charAt(0).toUpper match {
                            case 'A' => A = 0
                            case 'C' => A = 1
                            case 'G' => A = 2
                            case 'T' => A = 3 
                        }
                    col(i).trim().charAt(0).toUpper match {
                            case 'A' => B = 0
                            case 'C' => B = 1
                            case 'G' => B = 2
                            case 'T' => B = 3 
                        }
                    costs(A)(B) = Integer.parseInt(row(i))
                }

                line = HOXD.readLine()
            }

            return (a : Char, b : Char) => {
                var A = 0
                a.toUpper match {
                    case 'A' => A = 0
                    case 'C' => A = 1
                    case 'G' => A = 2
                    case 'T' => A = 3 
                }
                var B = 0
                b.toUpper match {
                    case 'A' => B = 0
                     case 'C' => B = 1
                    case 'G' => B = 2
                    case 'T' => B = 3 
                }
                costs(A)(B)
            }
    }

    def defaultHOXD : (Char,Char) => Int = {
        var costs : Array[Array[Int]] = Array.ofDim(4,4)

        costs(0)(0) = 91
        costs(0)(1) = -114
        costs(0)(2) = -31
        costs(0)(3) = -123

        costs(1)(0) = -114
        costs(1)(1) = 100
        costs(1)(2) = -125
        costs(1)(3) = -31

        costs(2)(0) = -31
        costs(2)(1) = -125
        costs(2)(2) = 100
        costs(2)(3) = -114

        costs(3)(0) = -123
        costs(3)(1) = -31
        costs(3)(2) = -114
        costs(3)(3) = 91

        return (a : Char, b : Char) => {
            var A = 0
            a.toUpper match {
                case 'A' => A = 0
                case 'C' => A = 1
                case 'G' => A = 2
                case 'T' => A = 3 
            }
            var B = 0
            b.toUpper match {
                case 'A' => B = 0
                case 'C' => B = 1
                case 'G' => B = 2
                case 'T' => B = 3 
            }
            costs(A)(B)
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
    

    def generateLocalAlignmentSet(maxL : Int, seqA : Sequence, seqS : Seq[Sequence],
                         settings : AlignSettings) : Seq[Alignment] = {
        
        val A = seqA.seq
        val out = mutable.Queue[Alignment]()
        val gO = settings.gapOpen
        val gE = settings.gapExtend
        val costFunc = settings.costFunc

	    var M : Array[Array[Int]] = Array.ofDim(A.length()+1,maxL+1)
        var X : Array[Array[Int]] = Array.ofDim(A.length()+1,maxL+1)
	    var Y : Array[Array[Int]] = Array.ofDim(A.length()+1,maxL+1)

        // only need to do this once for all the pirs
	    for (i <- 0 until A.length()){//Fills in zero row
	        M(i)(0) = 0
	        X(i)(0) = 0
            Y(i)(0) = settings.gapOpen + i * settings.gapExtend;
	    }

	    for (i <- 0 until maxL){//Fills in zero columns
            M(0)(i) = 0
	        X(0)(i) = settings.gapOpen + i * settings.gapExtend;
            Y(0)(i) = 0
	    }
    
        for ( seqB <- seqS ) {
        
        val B = seqB.seq

        var max = 0
        var maxLoc = (0,0)
        var t = 0
 
	    for (i <- 1 to A.length()){
	        for (j <- 1 to B.length()){
	    	//Based off of history, I created the matrix
    	        M(i)(j) = costFunc(A.charAt(i-1),B.charAt(j-1)) + 
                           math.max(math.max(M(i-1)(j-1),Y(i-1)(j-1)),
                                    math.max(X(i-1)(j-1),0))

                X(i)(j) = gE + math.max(math.max(M(i)(j-1) + gO,
                                                 Y(i)(j-1) + gO),
                                        math.max(X(i)(j-1),0))
           
                Y(i)(j) = gE + math.max(math.max(M(i-1)(j) + gO, Y(i-1)(j)),
                                        math.max(X(i-1)(j) + gO, 0))

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

        out += ( new Alignment(seqA,seqB,xSeq,ySeq,(i,j),opt,c,e))
        }
        
        return out
    }
}



