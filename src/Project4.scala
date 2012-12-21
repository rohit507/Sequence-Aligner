/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/
 
import BioLibs._

import java.util.Random
import java.io._

import scala.collection._
import scala.collection.immutable.StringOps
import scala.collection.mutable.OpenHashMap
import scala.collection.mutable.Queue
import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.UnrolledBuffer
import scala.actors._
import scala.actors.Futures._
import scala.math._

import gnu.trove.map.hash.TIntObjectHashMap

object Project4 {

    val rt = Runtime.getRuntime()
    var rand = new Random(System.currentTimeMillis())


/*  ======== Defaults =========
    min overlap length = 40 
    min % id= 98 or 99% 
    max region ignored= 60 
    gap open = -200, gap extend = -20.
    
    And I'd like you to continue to use HOXD 
    (or BLOSUM) for match/mismatch. */

    var debug = false
    var debugStop = -1
    var kSize = 12
    var action = "calc-overlaps"
    var input = ""
    var output = ""
    var stHash = false
    var stAlign = false
    var fdAlign = true
    var blockAlign = true
    
    def main(args: Array[String]) {

        var alignSettings = readArgs(args)
        BioLibs.debug = debug

        action match {
            case "calc-overlaps" =>
                var table = generateKmerTable(input)
                table.debug = debug
                calcOverlaps(genAlignment(table,alignSettings,true),alignSettings)
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
            case "bench-fasta-read" =>
                benchFastaRead(input)
                System.exit(0) 
            case "bench-kmer-gen" =>
                benchKmerGen(input)
                System.exit(0)
            case "bench-kmer-analysis" =>
                benchKmerAnalysis(input,alignSettings)
                System.exit(0)
            case "bench-align-quick" =>
                benchAlignQuick(input,alignSettings)
                System.exit(0) 
            case "bench-align" =>
                benchAlign(input,alignSettings)
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
        var minCollisions = 7
        var maxCollisions = 222
        var kmerCenter = 0.4f
        var kmerEdge = 0.4f

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
                case "-o" | "--output" => 
                    output = args(i + 1) 
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
               case "--kmer-center" =>
                    kmerCenter = math.abs(args(i + 1).toFloat)
                    i += 2
               case "--kmer-edge" =>
                    kmerEdge = math.abs(args(i + 1).toFloat)
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
                case "--quadratic-align" =>
                    fdAlign = false
                    i += 1 
                case "--linear-align" =>
                    fdAlign = true
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
                case "--bench-fasta-read" =>
                    action = "bench-fasta-read"
                    i += 1
                case "--bench-kmer-gen" =>
                    action = "bench-kmer-gen"
                    i += 1
                case "--bench-kmer-analysis" =>
                    action = "bench-kmer-analysis"
                    i += 1
                case "--bench-align-quick" =>
                    action = "bench-align-quick"
                    i += 1
                case "--bench-align" =>
                    action = "bench-align"
                    i += 1
                case "--debug" =>
                    debug = true
                    i += 1
                case "--sleep-for-debug" =>
                    println("Sleeping so debugger can connect.")
                    Thread.sleep(30000)
                    i += 1
                case s =>
                    println("Invalid Argument : " + s  )
                    println("Exiting Program.")
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
                                 minIdentity,maxIgnore,kSize,
                                 (minCollisions,maxCollisions),
                                 (kmerEdge,kmerCenter))
 
    }

    def printdb(s : String) { if(debug){ println(s)}}

    def printHelp() {
        println(" Rohit Ramesh :         Cmsc 423 \n" +
                " Project 4 : Sequence Overlapper \n" +
                "                                 \n" +
                " [Usage]                         \n" +
                "   See README file for details   \n")
    }

    /* */ // just prints the first few bits of data as the fasta file is being read
    def testFastaRead(file : String) {
        var i = 0
        println("")
        readSeq(file, (s : Sequence) => {
            i += 1
            if (i <= 10) {
                println("id : " + s.id)
                println("seq: " + s.seq)
                println("")
            } else {
                System.exit(0)
            }
        })
    } /* */

    /* */ // just reads a fasta file and prints a small amount of logging data. 
    def benchFastaRead(file : String) {
        val startTime = System.currentTimeMillis
        var i = 0
        readSeq(file, (s : Sequence) => {
            i += 1
        })
        println(" Read " + i + " sequences from " + file + " in " +
               (System.currentTimeMillis - startTime) + " milliseconds.")
    } /* */

    /* */ // tests a variety of kmer sizes and prints match data
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
    } /* */

    /* */ // Benchmarks both the single threaded and multithreaded
          //  kmer generation stuff. 
    def benchKmerGen(file : String) {
       try {
            var startTime = System.currentTimeMillis
            var tab = genSTKmerTable(file)
            println("\nGenerated " + tab.uniqueKmers() + " unique kmers from " +
                tab.uniqueSeqs() + " sequences from " + file + " sequentially in " +
               (System.currentTimeMillis - startTime) + " milliseconds.\n")
        } catch {
            case e => 
                println("\nSingle Threaded Kmer Generation Failed :\n\n")
                println(e.getMessage())
                e.printStackTrace()
        }
        try {
            var startTime = System.currentTimeMillis
            var tab = genMTKmerTable(file)
            println("Generated " + tab.uniqueKmers() + " unique kmers from " +
                tab.uniqueSeqs() + " sequences from " + file + " in parellel in " +
               (System.currentTimeMillis - startTime) + " milliseconds.\n")
        } catch {
            case e => 
                println("\nMulti Threaded Kmer Generation Failed :\n\n")
                println(e.getMessage())
                e.printStackTrace()
        }
    } /* */

    /* */ // Benchmarks both the single threaded and multithreaded
          //  kmer generation stuff. 
    def benchKmerAnalysis(file : String,s : AlignSettings) {
       try {
            println("Starting kmer gen.")
            var tab = generateKmerTable(file)
            tab.debug = debug 
            println("Finished kmer gen.")
            var startTime = System.currentTimeMillis
            tab.calcPairData(s)
            println("\nCalculated pair data in " + 
                (System.currentTimeMillis - startTime) + " milliseconds.\n")
            startTime = System.currentTimeMillis
            tab.calcDispatchData(s)
            println("Calculated dispatch data in " +
                (System.currentTimeMillis - startTime) + " milliseconds.\n")
        } catch {
            case e => 
                println("\nKmer Analysis Failed : \n\n")
                println(e.getMessage())
                e.printStackTrace()
        }
    } /* */

    /* */ // Tests whether collision dispatch is working correctly. 
    def testDispatchCollisions( table : KmerTable, settings : AlignSettings) {
        var i = 0
        var comp = mutable.HashSet[(Int,Int)]()
        table.dispatchCollisions(settings, (A : Sequence, B : Sequence) => {
            i += 1
            if(comp.contains((A.id,B.id))) {
                println( "!!!! Collission " + A.id + "<->" + B.id +
                         " Dispatched more than once. " )
            }
            comp += ((A.id,B.id))
            println( " Dispatched Coll : " + i + " - " + A.id + " <-> " + B.id)
        })
    } /*  */

    /* */ // Tests whether the dispatching of blocks works properly
    def testBlockDispatch( table : KmerTable , settings : AlignSettings) {
        var i = 0
        var comp = mutable.HashSet[(Int,Int)]()
        var hist = mutable.HashMap[Int,Int]()
        table.dispatchCollisionBlocks(settings,
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
    } /* */

    /* */ // Tests whether alignments are working properly in a nicely visual
          //  way. 
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
    } /* */

    /* */ // Benchmarks an alignment with whatever alignment
          // fuunction you ask it to
    def benchAlignHelper( file : String, s : AlignSettings, n : String,
             aligner : (KmerTable,AlignSettings,Boolean) => Seq[Alignment]) {        
       try {
            var tab = generateKmerTable(file)
            var startTime = System.currentTimeMillis
            var aligns = aligner(tab,s,false)
            println("\nCalculated " + aligns.size + " " + n + " alignments in " + 
                (System.currentTimeMillis - startTime) + " milliseconds.\n")
        } catch {
            case e => 
                println("\n" + n.capitalize + " Alignment Benchmark Failed : \n\n")
                println(e.getMessage())
                e.printStackTrace()
        }
    } /* */

    /* */ // Benchmarks both the single threaded and multithreaded
          //  kmer generation stuff, with small numbers of sequences
    def benchAlignQuick(file : String , s : AlignSettings) {
       debugStop = 500
       benchAlign(file,s)
    } /* */

    /* */ // Benchmarks both the single threaded and multithreaded
          //  kmer generation stuff. 
    def benchAlign(file : String , s : AlignSettings) {
       fdAlign = false
       benchAlignHelper(file,s,"single threaded quad single",genSingleSTAlign)
       benchAlignHelper(file,s,"single threaded quad block",genBlockSTAlign)
       benchAlignHelper(file,s,"multi threaded quad single",genSingleMTAlign)
       benchAlignHelper(file,s,"multi threaded quad block",genBlockMTAlign)
       fdAlign = true
       benchAlignHelper(file,s,"single threaded linear single",genSingleSTAlign)
       benchAlignHelper(file,s,"single threaded linear block",genBlockSTAlign)
       benchAlignHelper(file,s,"multi threaded linear single",genSingleMTAlign)
       benchAlignHelper(file,s,"multi threaded linear block",genBlockMTAlign)

    } /* */
 
    /* */ // Tests overlaps in a easy to read visual format
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
    } /* */

    /* */ // Generate a kmer table using whatever style the
    //        options tell you to
    def generateKmerTable(file : String) : KmerTable = {
        if (stHash) {
            return genSTKmerTable(file)
        } else {
            return genMTKmerTable(file)
        }
    } /* */

    /* */ // Generate a kmer table using a single thread
    def genSTKmerTable(file : String) : KmerTable = {
        var kmerTable = new KmerTable()
        var c = 0
        readSeq(file,(seq : Sequence) => {
           c += 1
           var hSet = generateKmerSet(kSize,seq)
           if((seq.id % 1000) == 0) {printdb("Generated Kmers for Seq : " + seq.id )}
           kmerTable.addKmerSet(seq,hSet)
        })

        return kmerTable
    } /* */

    /* */ // Generate a kmer table using multiple threads
    def genMTKmerTable(file : String) : KmerTable = {
        var kmerSets = Queue[Future[(Sequence,ArrayBuffer[Kmer])]]()
        var kmerTable = new KmerTable()
        readSeq(file,(seq : Sequence) => {
           printdb("Read Seq : " + seq.id)
           val f = future {
                if((seq.id % 1000) == 0) {
                    printdb("Start Processing Seq : " + seq.id)
                }
                var hSet = generateKmerSet(kSize,seq)
                if((seq.id % 1000) == 0) {
                    printdb("Generated Kmers for Seq : " + seq.id +
                        " on thread " + Thread.currentThread.getName ) 
                }
                (seq,hSet)
           }
           kmerSets += f
        })

        kmerSets.foreach { f =>
            f.apply() match {
                case (seq,set) => 
                    kmerTable.addKmerSet(seq,set)
                    if((seq.id % 1000) == 0) {
                        printdb("Added Kmer Set : " + seq.id)
                    }
                case _ => 
                    printdb("Error in Future Return")
            }
        }
        
        return kmerTable
    } /* */

    /* */ // Aligns all neccesary sequences in a table based on what
          //  options the user has chosen
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
      } /* */

    /* */ // Creates the alignment table/output with the settings
          // we chose, be it fast dovetail or otherwise. 
    def calcLocalAlignment( A : Sequence, B : Sequence,
                    settings : AlignSettings) : Alignment = {
        if (fdAlign) {
            return generateFastDovetailAlignment(A,B,settings)
        } else {
            return generateLocalAlignment(A,B,settings)  
        }
    }

    /* */ // Creates the alignment table/output with the settings
          // we chose, be it fast dovetail or otherwise, in set form 
    def calcLocalAlignmentSet(m : Int, A : Sequence, B : Seq[Sequence],
                    settings : AlignSettings) : Seq[Alignment] = {
        if (fdAlign) {
            return generateFastDovetailAlignmentSet(m,A,B,settings)
        } else {
            return generateLocalAlignmentSet(m,A,B,settings)  
        }
    }

    /* */ //Aligns Sequences one at a time in a single thread
    def genSingleSTAlign(table : KmerTable,settings : AlignSettings,
                     filter : Boolean) : Seq[Alignment] = {
        var aligns = mutable.Queue[Alignment]()
        table.dispatchCollisions(settings, (A : Sequence, B : Sequence) => {
            if((debugStop < 0) || (aligns.size > debugStop)) {
                if((aligns.size % 1000) == 0) {
                    printdb("Aligning pair : " + aligns.size)
                }
                if( (aligns.size % 5000) == 0) {
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
                var alg =  calcLocalAlignment(A,B,settings)    
                if ((! filter) || (alg.valid(settings))) {
                    aligns += alg
                }
            }
        })
        return aligns
    } /* */

    /* */ // Aligns Sequences one at a time, in multiple threads
    def genSingleMTAlign(table : KmerTable,settings : AlignSettings,
                     filter : Boolean) : Seq[Alignment] = {
        var futures = mutable.Queue[Future[Alignment]]()
        var aligns = mutable.Queue[Alignment]()
        var i = 0
        table.dispatchCollisions(settings, (A : Sequence, B : Sequence) => {
            if((debugStop < 0) || (i > debugStop)) {
                i += 1
                if((i % 1000) == 0) {
                    printdb("Dispatching Align : " + i)
                }
                futures += future { 
                    val c = i
                    if((c % 80) == 0) {
                        printdb("Starting Align " + c + 
                            " on thread " + Thread.currentThread.getName)
                    }
                    val a = calcLocalAlignment(A,B,settings)
                    if( (i % 8000) == 0) {
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
                    if((c % 80) == 0) {
                        printdb("Completed Align : " + c)
                    }
                    a
                }
            }
        })
        i = 0
        for (f <- futures) {
            var al = f.apply()
            if ((! filter) || (al.valid(settings))) {
                i += 1
                if((i % 1000) == 0) {
                    printdb("Integrating Align : " + i)
                }
                aligns += al
            }
        }

        return aligns
    } /* */

    /* */ // Aligns Blocks of sequences in a single thread.
    def genBlockSTAlign(table : KmerTable,settings : AlignSettings,
                     filter : Boolean) : Seq[Alignment] = {
        var aligns = mutable.UnrolledBuffer[Alignment]()
        var c = 0
        table.dispatchCollisionBlocks(settings, 
                                (max : Int, A : Sequence, S : Seq[Sequence]) => {
            if((debugStop < 0) || (aligns.size > debugStop)) {
                c += 1
                if((c % 100) == 0) {
                    printdb("Aligning Block " + c + 
                        " at sequence " + aligns.size)
                }   
                if( (c % 1000) == 0) {
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
                var alg =  calcLocalAlignmentSet(max,A,S,settings)
                for (a <- alg) {
                    if ((! filter) || (a.valid(settings))) {
                        aligns += a
                    }
                }
            }
        })
        printdb("There are " + aligns.size + " probable valid alignments.") 
        return aligns

    } /* */

    /* */ // Aligns Blocks of Sequences in multiple threads
    def genBlockMTAlign(table : KmerTable,settings : AlignSettings,
                     filter : Boolean) : Seq[Alignment] = {
        var futures = mutable.Queue[Future[Seq[Alignment]]]()
        var aligns = mutable.Queue[Alignment]()
        var i = 0
        table.dispatchCollisionBlocks(settings,
                                (max : Int, A : Sequence, S : Seq[Sequence]) => {
            if((debugStop < 0) || (i > debugStop)) {
                i += 1
                if((i % 1000) == 0) {
                    printdb("Dispatching Align block : " + i)
                }
                if( (i % 5000) == 0) {
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
                futures += future {
                    val c = i 
                    if((c % 100) == 0) {
                        printdb("Starting Align Block " + i + 
                            " on thread " + Thread.currentThread.getName)
                    }
                    val a = calcLocalAlignmentSet(max,A,S,settings)
                    if((c % 100) == 0) {
                        printdb("Completed Align Block " + i + 
                            " on thread " + Thread.currentThread.getName)
                    }
                    a
                }
            }
        })

        i = 0
        for (f <- futures) {
            var alg = f.apply()  
            for (a <- alg) {
                if ((! filter) || (a.valid(settings))) {
                    i += 1
                    if((i % 1000) == 0) {
                        printdb("Integrating Alignment : " + i)
                    }
                    if( (i % 2000) == 0) {
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
                    aligns += a
                }
            }
        }

        return aligns
    } /* */


    /* */ // calculates overlaps in the standar fashing and prints 
          //  to where it should
    def calcOverlaps(alignments : Seq[Alignment],settings : AlignSettings) {
        var i = 0
        val file = new File(output)
        var oStream : FileWriter = null
        if (output != "") {
            if(file.exists()){
                file.delete()
            }
            file.createNewFile()
            oStream = new FileWriter(file)
        }
        var out = ""
        for (a <- alignments) {
            var o = a.getOverlap()
            if( o.valid(settings)) {
                i += 1
                if((i % 1000) == 0) {
                    printdb("Saving Overlap : " + i)
                }
                out = o.print() 
                if (output == ""){
                    println(out)
                }else {
                    oStream.write(out + "\n",0,out.length + 1)
                }
            }                
        }
        if (output != "") {
            oStream.flush()
        }
    } /* */

    
}
