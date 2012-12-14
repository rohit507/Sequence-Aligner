/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/
 
import BioLibs._

import java.util.Random

import scala.collection._
import scala.collection.immutable.StringOps
import scala.actors._
import scala.actors.Futures._
import scala.math._

object Project4 {

    var rand = new Random(System.currentTimeMillis())
    var bounds = (2,20)

/*  ======== Defaults =========
    min overlap length = 40 
    min % id= 98 or 99% 
    max region ignored= 60 
    gap open = -200, gap extend = -20.
    
    And I'd like you to continue to use HOXD 
    (or BLOSUM) for match/mismatch. */

    var debug = false
    var kSize = 12
    var action = "calc-overlaps"
    var input = ""
    var stHash = true
    var stAlign = false
    var blockAlign = true
    
    def main(args: Array[String]) {
        var alignSettings = readArgs(args)

        action match {
            case "calc-overlaps" =>
                var table = generateKmerTable(input)
                println(calcOverlaps(genAlignment(table,alignSettings,true),alignSettings))
                System.exit(0) 
            case "test-kmer-cover" => 
                testKmerCover(input)
                System.exit(0)
            case "test-dispatch-collisions" => 
                var table = generateKmerTable(input)
                testDispatchCollisions(table,alignSettings)
                System.exit(0)
            case "test-block-dispatch" => 
                var table = generateKmerTable(input)
                testBlockDispatch(table,alignSettings)
                System.exit(0)
            case "test-alignment" => 
                var table = generateKmerTable(input)
                testAlignment(genAlignment(table,alignSettings,false),alignSettings)
                System.exit(0) 
            case "test-overlaps" => 
                var table = generateKmerTable(input)
                testOverlaps(genAlignment(table,alignSettings,false),alignSettings)
                System.exit(0) 
            case "test-fasta-read" => 
                testFastaRead(input)
                System.exit(0)
        }
    }

    def readArgs(args: Array[String]) : AlignSettings = {
        var i = 0
        var hoxd = ""
        var mat = 95
        var mism = -70
        var minOverlap = 40
        var minIdentity = 0.98f
        var maxIgnore = 90
        var gapOpen = -200
        var gapExtend = -20
        var minCollisions = 2
        var maxCollisions = 22

        while (i < args.length) {
            args(i) match {
                case "-h" | "--help" =>
                    printHelp()
                    System.exit(0)
                case "-m" | "--matrix" | "-H" | "--HOXD-matrix" => 
                    hoxd = args(i + 1) // Matrix
                    i += 2
                case "-k" | "--kmer-size" => 
                    kSize = args(i + 1).toInt // kmers
                    i += 2
                case "-i" | "--input" => 
                    input = args(i + 1) 
                    i += 2
                case "--match" =>
                    mat = math.abs(args(i + 1).toInt)
                    i += 2
                case "--mismatch" =>
                    mism = -math.abs(args(i + 1).toInt)
                    i += 2
                case "--min-overlap" => 
                    minOverlap = math.abs(args(i + 1).toInt)
                    i += 2
                case "--min-identity" => 
                    minIdentity = args(i + 1).toFloat
                    if (minIdentity >= 1) {
                        minIdentity *= .01f
                    }
                    i += 2
                case "--min-collisions" =>
                    minCollisions = math.abs(args(i + 1).toInt)
                    i += 2
                case "--max-collisions" =>
                    maxCollisions = math.abs(args(i + 1).toInt)
                    i += 2
                case "-gO" | "--gap-open" =>
                    gapOpen = -math.abs(args(i + 1).toInt)
                    i += 2
                case "-gE" | "--gap-extend" =>
                    gapExtend = -math.abs(args(i + 1).toInt)
                    i += 2
                case "--max-ignore" =>
                    maxIgnore = math.abs(args(i + 1).toInt)
                    i += 2
                case "--st-hash" =>
                    stHash = true
                    i += 1
                case "--mt-hash" =>
                    stHash = false
                    i += 1
                case "--st-align" =>
                    stAlign = true
                    i += 1
                case "--mt-align" =>
                    stAlign = false
                    i += 1
                case "--block-align" =>
                    blockAlign = true
                    i += 1
                case "--single-align" =>
                    blockAlign = false
                    i += 1
                case "--calc-overlaps" =>
                    action = "calc-overlaps"
                    i += 1
                case "--test-overlaps" =>
                    action = "test-overlaps"
                    i += 1
                case "--test-alignment" =>
                    action = "test-alignment"
                    i += 1
                case "--test-dispatch-collisions" =>
                    action = "test-dispatch-collisions"
                    i += 1
                case "--test-block-dispatch" =>
                    action = "test-block-dispatch"
                    i += 1
                case "--test-kmer-cover" =>
                    action = "test-kmer-cover"
                    i += 1
                case "--test-fasta-read" =>
                    action = "test-fasta-read"
                    i += 1
                case "--debug-all" =>
                    debug = true
                    i += 1 
                case _ =>
                    println("Invalid Argument")
                    System.exit(1)
            }
        }

        var compFunc = defaultHOXD

        if (hoxd != "") {
            compFunc = readHOXD(hoxd) 
        }

        if (input == "") {
            println("No input file specified")
            System.exit(-1)
        }

        return new AlignSettings(compFunc,gapOpen,gapExtend,minOverlap,
                                 minIdentity,maxIgnore,minCollisions,maxCollisions)
 
    }

    def printHelp() {
        println(" Rohit Ramesh :         Cmsc 423 \n" +
                " Project 4 : Sequence Overlapper \n" +
                "                                 \n" +
                " [Usage]                         \n" +
                "   See README file for details   \n")
    }

    // just prints the first few bits of data as the fasta file is being read
    def testFastaRead(file : String) {
        var i = 0
        readSeq(file, (s : Sequence) => {
            i += 1
            if (i <= 10) {
                println("id : " + s.id)
                println("seq: " + s.seq)
                println("")
            }
        })
    } 

    // tests a variety of kmer sizes and prints match data
    def testKmerCover( file : String ) {
        for (i <- 0 to 25) {
            val possible = math.pow(4,i) 
            kSize = i
            val tab = generateKmerTable(file)
            val uniques = tab.uniqueKmers()
            val ratio = (uniques.toFloat) / (possible.toFloat) 
            println("Kmer Size : " + i )
            println("  uniques : " + uniques )
            println("  ratio   : " + ratio )
            var collisions = tab.kmerCollisionHistogram()
            var keys = collisions.keySet.toArray.sortWith(_<_)
            var o = ""
            for (k <- keys ) {
                o = o + "          [" + k + " -> " +
                      collisions.apply(k) + "]\n"
            }
            println("  [ number of collisions -> count of " +
                    "seqs with that many collisions ] :\n" + o )
            
        }
    }

    def testDispatchCollisions( table : KmerTable, settings : AlignSettings) {
        var i = 0
        var comp = mutable.HashSet[(Int,Int)]()
        table.dispatchCollisions(settings.collBounds, (A : Sequence, B : Sequence) => {
            i += 1
            if(comp.contains((A.id,B.id))) {
                println( "!!!! Collission " + A.id + "<->" + B.id +
                         " Dispatched more than once. " )
            }
            comp += ((A.id,B.id))
            println( " Dispatched Coll : " + i + " - " + A.id + " <-> " + B.id)
        })
    }

    def testBlockDispatch( table : KmerTable , settings : AlignSettings) {
        var i = 0
        var comp = mutable.HashSet[(Int,Int)]()
        var hist = mutable.HashMap[Int,Int]()
        table.dispatchCollisionBlocks(settings.collBounds,
                     (max : Int, A : Sequence, S : Seq[Sequence]) => {
            for(B <- S) {
                i += 1
                if(comp.contains((A.id,B.id))) {
                    println( "!!!! Collission " + A.id + "<->" + B.id +
                            " Dispatched more than once. " )
                }
                comp += ((A.id,B.id))
                println( " Dispatched Coll : " + i + " - " + A.id + " <-> " + B.id)
            }
            if(!hist.contains(S.size)){
                hist += ((S.size,0))
            }
            hist += ((S.size,hist.apply(S.size) + 1))
        })
                                // :
        println( "\n Histogram Of Relations : [Number of Aligns -> " +
                 "Number of Seqs w/ that many Aligns]")
        var keys = hist.keySet.toArray.sortWith(_<_)
        var o = ""
        for (k <- keys ) {
            o = o + "          [" + k + " -> " +
                      hist.apply(k) + "]\n"
        }
        println(o)
    }

    def testAlignment(alignments : Seq[Alignment],settings : AlignSettings) {
        var i = 0;
        for (a <- alignments) {
            i += 1
            println(" Alignment " + i + " : " + a.seqA.id + " <-> " + a.seqB.id )
            println("   Seq A     : " + a.seqA.seq)
            println("   Seq B     : " + a.seqB.seq)
            println("   Overlap A : " + a.alignA)
            println("   Overlap B : " + a.alignB)
            println("   Start     : " + a.start)
            println("   End       : " + a.end)
            println("   Error Rat : " + a.errRatio)
            println("   is Valid? : " + a.valid(settings))
            println("")                       
        }
    }

    def testOverlaps(alignments : Seq[Alignment],settings : AlignSettings) {
        var i = 0;
        for (a <- alignments) {
            var o = a.getOverlap()
            i += 1
            println(" Overlap " + i + " : " + a.seqA.id + " <-> " + a.seqB.id )
            if (o.ahg >= 0) {
                println("   Seq A   : " + a.seqA.seq + ("").padTo(o.bhg,'-'))
                println("   Seq B   : " + ("").padTo(o.ahg,'-') + a.seqB.seq )
            } else {
                println("   Seq A   : " + ("").padTo(-o.ahg,'-') + a.seqA.seq)
                println("   Seq B   : " + a.seqB.seq + ("").padTo(-o.bhg,'-'))
            }
            println("   Ahg     : " + o.ahg)
            println("   Bhg     : " + o.bhg)
            println("   Start   : " + a.start)
            println("   End     : " + a.end)
            println("   Error   : " + a.errRatio)
            println("   Valid?  : " + o.valid(settings))                       
        }

    }

    def generateKmerTable(file : String) : KmerTable = {
        if (stHash) {
            return genSTKmerTable(file)
        } else {
            return genMTKmerTable(file)
        }
    }

    def genSTKmerTable(file : String) : KmerTable = {
        var kmerTable = new KmerTable()

        readSeq(file,(seq : Sequence) => {
           var hSet = generateKmerSet(kSize,seq.seq)
           if (debug) { println("Generated Kmers for Seq : " + seq.id )}
           kmerTable.addKmerSet(seq,hSet)
           if (debug) { println("Added Kmer Set : " + seq.id) }
        })

        return kmerTable
    }

    def genMTKmerTable(file : String) : KmerTable = {
        var kmerSets = mutable.Queue[Future[(Sequence,Set[String])]]()
        var kmerTable = new KmerTable()

        readSeq(file,(seq : Sequence) => {
           if (debug) { println("Read Seq : " + seq.id) }
           val f = future {
                if (debug) { println("Start Processing Seq : " + seq.id)}
                var hSet = generateKmerSet(kSize,seq.seq)
                if (debug) { println("Generated Kmers for Seq : " + seq.id +
                                " on thread " + Thread.currentThread.getName ) }
                (seq,hSet)
           }
           kmerSets += f
        })

        kmerSets.foreach { f =>
            f.apply() match {
                case (seq,set) => 
                    kmerTable.addKmerSet(seq,set)
                    if (debug) { println("Added Kmer Set : " + seq.id)}
                case _ => 
                    if (debug) { println("Error in Future Return")}
            }
        }
        
        return kmerTable
    }

    def genAlignment(table : KmerTable,settings : AlignSettings,
                        filter : Boolean) : Seq[Alignment] = {
        if (blockAlign) {
            if (stAlign) {
                return genBlockSTAlign(table,settings,filter)
            } else {
                return genBlockMTAlign(table,settings,filter)
            }
        } else {
            if (stAlign) {
                return genSingleSTAlign(table,settings,filter)
            } else {
                return genSingleMTAlign(table,settings,filter)
            }
        }
      }

    def genSingleSTAlign(table : KmerTable,settings : AlignSettings,
                     filter : Boolean) : Seq[Alignment] = {
        var aligns = mutable.Queue[Alignment]()
        table.dispatchCollisions(settings.collBounds, (A : Sequence, B : Sequence) => {
            var alg =  generateLocalAlignment(A,B,settings)    
            if ((! filter) || (alg.valid(settings))) {
                aligns += alg
            }
        })
        return aligns
    }

    def genSingleMTAlign(table : KmerTable,settings : AlignSettings,
                     filter : Boolean) : Seq[Alignment] = {
        var futures = mutable.Queue[Future[Alignment]]()
        var aligns = mutable.Queue[Alignment]()
        table.dispatchCollisions(settings.collBounds, (A : Sequence, B : Sequence) => {
            futures += future { generateLocalAlignment(A,B,settings)}
        })

        for (f <- futures) {
            var al = f.apply()
            if ((! filter) || (al.valid(settings))) {
                aligns += al
            }
        }

        return aligns
    }

    def genBlockSTAlign(table : KmerTable,settings : AlignSettings,
                     filter : Boolean) : Seq[Alignment] = {
        var aligns = mutable.Queue[Alignment]()
        table.dispatchCollisionBlocks(settings.collBounds, 
                                (max : Int, A : Sequence, S : Seq[Sequence]) => {
            var alg =  generateLocalAlignmentSet(max,A,S,settings)
            for (a <- alg) {
                if ((! filter) || (a.valid(settings))) {
                    aligns += a
                }
            }
        })
        return aligns

    }

    def genBlockMTAlign(table : KmerTable,settings : AlignSettings,
                     filter : Boolean) : Seq[Alignment] = {
        var futures = mutable.Queue[Future[Seq[Alignment]]]()
        var aligns = mutable.Queue[Alignment]()
        table.dispatchCollisionBlocks(settings.collBounds,
                                (max : Int, A : Sequence, S : Seq[Sequence]) => {
            futures += future {generateLocalAlignmentSet(max,A,S,settings)}
        })

        for (f <- futures) {
            var alg = f.apply()  
            for (a <- alg) {
                if ((! filter) || (a.valid(settings))) {
                    aligns += a
                }
            }
        }

        return aligns
    }


    def calcOverlaps(alignments : Seq[Alignment],settings : AlignSettings) : String = {
        var out = ""
        for (a <- alignments) {
            var o = a.getOverlap()
            if( o.valid(settings)) {
                out += o.print() + "\n"
            }                
        }
        return out
    }

    
}
