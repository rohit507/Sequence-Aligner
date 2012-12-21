/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/

import java.io._

import scala.collection._
import scala.collection.mutable.OpenHashMap
import scala.collection.mutable.ArrayBuffer
import scala.actors._
import scala.actors.Futures._
import scala.collection.immutable.StringOps

object BioLibs {

    var debug = false
    private def printdb(s : String) { if(debug){ println(s)}}

    // Dud alignment that represents a generic alignment failure
    val dud = new Alignment(new Sequence(0,""),new Sequence(0,""),"","",(0,0),(0,0),0,1)

    /* */ // Read in a fasta file dispatching each sequence as it's
          // read through the input function. 
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
    } /* */ 

    /* */ // generates an array of kmers, which contain location
          //  data and can be easily expanded to conatin more if needed
    def generateKmerSet( k : Int, s : Sequence ) : (ArrayBuffer[Kmer]) = {
        var o = new mutable.ArrayBuffer[Kmer](s.seq.length - k)
        var d = (s.seq.length - k).toFloat
        for( i <- 0 to (s.seq.length - k)) {
            o += new Kmer(s.id,s.seq.substring(i,i+k),(i.toFloat)/d)
        }
        return o
    } /* */

    /* */ // Read in an HOXD file making a specific effort to store it
          //  in a fast data structure, so the many millions of calls
          //  to the cost function will happen reasonably quickly. 
    def readHOXD( file : String ) : (Char,Char) => Int = {

        val costs : Array[Array[Int]] = Array.ofDim(4,4)
        val HOXD = new BufferedReader(new FileReader(file))
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
            b.toUpper match {
                case 'A' => 
                    costs(A)(0)
                case 'C' => 
                    costs(A)(1)
                case 'G' =>
                    costs(A)(2)
                case 'T' =>
                    costs(A)(3)
            }
        }
    } /* */
    
    /* */ // Same as the readHOXD function but with the data entered in 
          //  by hand because I can't be arsed to type in a few extra 
          //  arguments at runtime. 
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
            b.toUpper match {
                case 'A' => 
                    costs(A)(0)
                case 'C' => 
                    costs(A)(1)
                case 'G' =>
                    costs(A)(2)
                case 'T' =>
                    costs(A)(3)
            }
        }
    } /* */ 

    /* */ // A simple cost function that only has specific match
          //  and mismatch values, noting else happens.
    def simpleMatch( mat : Int, miss : Int) : (Char,Char) => Int = {
        return (a : Char, b : Char) => { if (a == b) mat else miss }
    } /* */

    /* */ // Standard 3 matrix local alignment function, pulled
          //  almost stright from project 3.
    def generateLocalAlignment(seqA : Sequence, seqB : Sequence,
                         settings : AlignSettings) : Alignment = {
        
        val A = seqA.seq
        val B = seqB.seq

	    val M : Array[Array[Int]] = Array.ofDim(A.length()+1,B.length()+1)
        val X : Array[Array[Int]] = Array.ofDim(A.length()+1,B.length()+1)
	    val Y : Array[Array[Int]] = Array.ofDim(A.length()+1,B.length()+1)

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
    } /* */
    
    /* */ // Generates a set of star alignments while taking care
          // to reuse memory and other information
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
    } /* */

    /* */ // Uses the fast linear dovetail alignment algorithm to quickly
          //  align two sequences if there is a reasonable expectation
          //  sequence A comes before sequence B in a dovetail. 
    def generateFastDovetailAlignment(seqA : Sequence, seqB : Sequence,
                         settings : AlignSettings) : Alignment = {

            // Fast dovetail alignment is a 2 part algorithm, first 
            //  we perform a local alignment of one sequence with the
            //  first kmer of the second sequence. Then that alignment
            //  tells us where to start the second, limited distance 
            //  from match, sequence alignment to check if the possible
            //  dovetail is a good one. 

        val A = seqA.seq
        val B = seqB.seq

        // Set up arrays we'll reuse for the whole process, which means
        //  figure out how wide/large they should be. 

        val width : Int = math.max(settings.kmerSize,
                        math.floor(A.size * (1 - settings.minIdentity)).toInt + 1)

	    val M : Array[Array[Int]] = Array.ofDim(A.length()+1,width+1)
        val X : Array[Array[Int]] = Array.ofDim(A.length()+1,width+1)
	    val Y : Array[Array[Int]] = Array.ofDim(A.length()+1,width+1)

        // The actual alignment routine for the first sequence. 
        //  and the second kmer of the first sequence. 
    
	    for (i <- 0 until A.length()){//Fills in zero row
	        M(i)(0) = 0
	        X(i)(0) = 0
            Y(i)(0) = settings.gapOpen + i * settings.gapExtend;
	    }

	    for (i <- 0 until width){//Fills in zero columns
            M(0)(i) = 0
	        X(0)(i) = settings.gapOpen + i * settings.gapExtend;
            Y(0)(i) = 0
	    }

        var max = 0
        var maxLoc = (0,0)
        var t = 0

	    for (i <- 1 to A.length()){
	        for (j <- 1 to width){
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

        // Backtrack and retrieve the best local alignment for the 
        // first kmer

        var opt = maxLoc
        var i = maxLoc._1 
        var j = maxLoc._2

        max = math.max(M(i)(j),math.max(X(i)(j),Y(i)(j)))

        do {
            if (M(i)(j) == max) {
                i -= 1
                j -= 1
            } else if (X(i)(j) == max) {
                j -= 1
            } else if (Y(i)(j) == max) {
                i -= 1
            }
            max = math.max(M(i)(j),math.max(X(i)(j),Y(i)(j)))
        } while(max > 0)

        // If that local alignment doesn't end at the first element of
        //  the second string, it's a dud. 

        if ( j != 0 ) {
            return dud 
        }

        // otherwise we can continue and perform the second, slightly
        //  more elaborate alignment. 

        // start by figuring off the various offsets we want. 

        val doveStart = i   
        val doveLength = A.length - doveStart 
        var zeroRow = width / 2 

        // And some basic variables
    
        max = 0
        maxLoc = (0,0)
        t = 0

        // Now, it must be noted this is *exactly* local sequence alignment 
        //  algorithm with a simple change of coordinates. 
        
        // Namely u,k and are i,j shifted so a diagonal line in ij is
        //  horizontal in uk

        // U = I - doveStart 
        // K = J + zeroRow - I + doveStart

        // I = U + doveStart 
        // J = K - zeroRow + U 

	    for (u <- 0 to doveLength){
	        for (k <- 0 to width){
                i = u + doveStart 
                j = k - zeroRow + u

                //printdb("U: "+u+" K: "+k+" I: "+i+" J: "+j)
                if ((i <= doveStart) || (j <= 0) || (j > B.length)) { // Initial rows
                    M(u)(k) = 0
	                X(u)(k) = 0
                    Y(u)(k) = 0
                } else {

                    if(u != 0) {
    	            M(u)(k) = settings.costFunc(A.charAt(i-1),B.charAt(j-1)) + 
                                math.max(math.max(M(u-1)(k),Y(u-1)(k)),
                                    math.max(X(u-1)(k),0))
                    } else { M(u)(k) = 0 }

                    if(k != 0) {
                    X(u)(k) = settings.gapExtend + 
                                math.max(math.max(M(u)(k-1) + settings.gapOpen,
                                             Y(u)(k-1) + settings.gapOpen),
                                    math.max(X(u)(k-1),0))
                    } else { X(u)(k) = 0 }

                    if ((u != 0) && (k != width)) {
                    Y(u)(k) = settings.gapExtend + 
                                math.max(math.max(M(u-1)(k+1) + settings.gapOpen, Y(u-1)(k+1)),
                                    math.max(X(u-1)(k+1) + settings.gapOpen, 0))
                    } else { Y(u)(k) = 0 }
                }

                t = math.max(M(u)(k),math.max(X(u)(k),Y(u)(k)))

                if (t > max) {
                    max = t
                    maxLoc = (u,k)
                }
            }
	    }

        // and backtrack out of it

        opt = maxLoc
        var u = maxLoc._1 
        var k = maxLoc._2
        var xSeq = "" 
        var ySeq = ""
        var c = 0
        var e = 0
        var pa = ' '
        var pb = ' '

        max = math.max(M(u)(k),math.max(X(u)(k),Y(u)(k)))

        do {
            i = u + doveStart 
            j = k - zeroRow + u

            if (M(u)(k) == max) {
                pa = A.charAt(i-1)
                pb = B.charAt(j-1)
                u -= 1
            } else if (X(u)(k) == max) {
                pa = A.charAt(i-1)
                pb = '-'
                k -= 1
            } else if (Y(u)(k) == max) {
                pa = '-'
                pb = B.charAt(j-1)
                u -= 1
                k += 1
            }

            if(pa != pb) {
                e += 1
            } else {
                c += 1
            } 
 
            xSeq = pa + xSeq
            ySeq = pb + ySeq
      
            max = math.max(M(u)(k),math.max(X(u)(k),Y(u)(k)))
        } while(max > 0)

        // Get start and end coordinate in ij space
        i = u + doveStart 
        j = k - zeroRow + u
        val newEnd = (opt._1 + doveStart ,opt._2 - zeroRow + opt._1)

        // We can just use this, and since everything is back in IJ space
        //  like it should be, we don't have to refactor our validity
        //  rules among other things.
        return new Alignment(seqA,seqB,xSeq,ySeq,(i,j),newEnd,c,e)

    }   

    /* */ // Uses the fast linear dovetail alignment algorithm to quickly
          //  align two sequences if there is a reasonable expectation
          //  sequence A comes before sequence B in a dovetail. 
    def generateFastDovetailAlignmentSet(maxL : Int, seqA : Sequence, seqS : Seq[Sequence],
                         settings : AlignSettings) : Seq[Alignment] = {
        
        val out = mutable.Queue[Alignment]()
        val gO = settings.gapOpen
        val gE = settings.gapExtend
        val costFunc = settings.costFunc

            // Fast dovetail alignment is a 2 part algorithm, first 
            //  we perform a local alignment of one sequence with the
            //  first kmer of the second sequence. Then that alignment
            //  tells us where to start the second, limited distance 
            //  from match, sequence alignment to check if the possible
            //  dovetail is a good one. 

        val A = seqA.seq

        for (seqB <- seqS) {
            val B = seqB.seq

        // Set up arrays we'll reuse for the whole process, which means
        //  figure out how wide/large they should be. 

        val width : Int = math.max(settings.kmerSize,
                        math.floor(A.size * (1 - settings.minIdentity)).toInt + 1)

	    val M : Array[Array[Int]] = Array.ofDim(A.length()+1,width+1)
        val X : Array[Array[Int]] = Array.ofDim(A.length()+1,width+1)
	    val Y : Array[Array[Int]] = Array.ofDim(A.length()+1,width+1)

        // The actual alignment routine for the first sequence. 
        //  and the second kmer of the first sequence. 
    
	    for (i <- 0 until A.length()){//Fills in zero row
	        M(i)(0) = 0
	        X(i)(0) = 0
            Y(i)(0) = gO + i * gE;
	    }

	    for (i <- 0 until width){//Fills in zero columns
            M(0)(i) = 0
	        X(0)(i) = gO + i * gE;
            Y(0)(i) = 0
	    }

        var max = 0
        var maxLoc = (0,0)
        var t = 0

	    for (i <- 1 to A.length()){
	        for (j <- 1 to width){
	    	//Based off of history, I created the matrix
    	        M(i)(j) = settings.costFunc(A.charAt(i-1),B.charAt(j-1)) + 
                           math.max(math.max(M(i-1)(j-1),Y(i-1)(j-1)),
                                    math.max(X(i-1)(j-1),0))

                X(i)(j) = gE + 
                           math.max(math.max(M(i)(j-1) + gO,
                                             Y(i)(j-1) + gO),
                                    math.max(X(i)(j-1),0))
           
                Y(i)(j) = gE + 
                           math.max(math.max(M(i-1)(j) + gO, Y(i-1)(j)),
                                    math.max(X(i-1)(j) + gO, 0))

                t = math.max(M(i)(j),math.max(X(i)(j),Y(i)(j)))

                if (t > max) {
                    max = t
                    maxLoc = (i,j)
                }
            }
	    }

        // Backtrack and retrieve the best local alignment for the 
        // first kmer

        var opt = maxLoc
        var i = maxLoc._1 
        var j = maxLoc._2

        max = math.max(M(i)(j),math.max(X(i)(j),Y(i)(j)))

        do {
            if (M(i)(j) == max) {
                i -= 1
                j -= 1
            } else if (X(i)(j) == max) {
                j -= 1
            } else if (Y(i)(j) == max) {
                i -= 1
            }
            max = math.max(M(i)(j),math.max(X(i)(j),Y(i)(j)))
        } while(max > 0)

        // If that local alignment doesn't end at the first element of
        //  the second string, it's a dud. 

        if ( j != 0 ) {
            out += dud 
        } else {

        // otherwise we can continue and perform the second, slightly
        //  more elaborate alignment. 

        // start by figuring off the various offsets we want. 

        val doveStart = i   
        val doveLength = A.length - doveStart 
        var zeroRow = width / 2 

        // And some basic variables
    
        max = 0
        maxLoc = (0,0)
        t = 0

        // Now, it must be noted this is *exactly* local sequence alignment 
        //  algorithm with a simple change of coordinates. 
        
        // Namely u,k and are i,j shifted so a diagonal line in ij is
        //  horizontal in uk

        // U = I - doveStart 
        // K = J + zeroRow - I + doveStart

        // I = U + doveStart 
        // J = K - zeroRow + U 

	    for (u <- 0 to doveLength){
	        for (k <- 0 to width){
                i = u + doveStart 
                j = k - zeroRow + u

                //printdb("U: "+u+" K: "+k+" I: "+i+" J: "+j)
                if ((i <= doveStart) || (j <= 0) || (j > B.length)) { // Initial rows
                    M(u)(k) = 0
	                X(u)(k) = 0
                    Y(u)(k) = 0
                } else {

                    if(u != 0) {
    	            M(u)(k) = settings.costFunc(A.charAt(i-1),B.charAt(j-1)) + 
                                math.max(math.max(M(u-1)(k),Y(u-1)(k)),
                                    math.max(X(u-1)(k),0))
                    } else { M(u)(k) = 0 }

                    if(k != 0) {
                    X(u)(k) = settings.gapExtend + 
                                math.max(math.max(M(u)(k-1) + settings.gapOpen,
                                             Y(u)(k-1) + settings.gapOpen),
                                    math.max(X(u)(k-1),0))
                    } else { X(u)(k) = 0 }

                    if ((u != 0) && (k != width)) {
                    Y(u)(k) = settings.gapExtend + 
                                math.max(math.max(M(u-1)(k+1) + settings.gapOpen, Y(u-1)(k+1)),
                                    math.max(X(u-1)(k+1) + settings.gapOpen, 0))
                    } else { Y(u)(k) = 0 }
                }

                t = math.max(M(u)(k),math.max(X(u)(k),Y(u)(k)))

                if (t > max) {
                    max = t
                    maxLoc = (u,k)
                }
            }
	    }

        // and backtrack out of it

        opt = maxLoc
        var u = maxLoc._1 
        var k = maxLoc._2
        var xSeq = "" 
        var ySeq = ""
        var c = 0
        var e = 0
        var pa = ' '
        var pb = ' '

        max = math.max(M(u)(k),math.max(X(u)(k),Y(u)(k)))

        do {
            i = u + doveStart 
            j = k - zeroRow + u

            if (M(u)(k) == max) {
                pa = A.charAt(i-1)
                pb = B.charAt(j-1)
                u -= 1
            } else if (X(u)(k) == max) {
                pa = A.charAt(i-1)
                pb = '-'
                k -= 1
            } else if (Y(u)(k) == max) {
                pa = '-'
                pb = B.charAt(j-1)
                u -= 1
                k += 1
            }

            if(pa != pb) {
                e += 1
            } else {
                c += 1
            } 
 
            xSeq = pa + xSeq
            ySeq = pb + ySeq
      
            max = math.max(M(u)(k),math.max(X(u)(k),Y(u)(k)))
        } while(max > 0)

        // Get start and end coordinate in ij space
        i = u + doveStart 
        j = k - zeroRow + u
        val newEnd = (opt._1 + doveStart ,opt._2 - zeroRow + opt._1)

        // We can just use this, and since everything is back in IJ space
        //  like it should be, we don't have to refactor our validity
        //  rules among other things.
        out += new Alignment(seqA,seqB,xSeq,ySeq,(i,j),newEnd,c,e)
        }}
        return out
    }
}

