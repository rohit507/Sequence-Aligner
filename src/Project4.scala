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
    var debugAll = false

    def main(args: Array[String]) {
        var alignSettings = readArgs(args)

       /* action match {
            case "calc-overlaps" =>
                calculateOverlaps(input,alignSettings, (o: Overlap) => {
                   if(o.valid(alignSettings)) {
                        println(o.print())
                   }
                })
           case "test-overlaps" =>
                if (stAlign) {
                    testOverlap(input,alignSettings,debugAll)
                } else {
                    testFutureOverlap(input,alignSettings,debugAll)
                }
                System.exit(0) 
            case "test-kmer-cover" => 
                testKmerCover(input)
                System.exit(0)
            case "test-dispatch-collisions" => */
                var table = generateKmerTable(input)
                testDispatchCollisions(table,alignSettings)
                System.exit(0)
          /*  case "test-fasta-read" => 
                testFastaRead(input)
                System.exit(0)
        }*/
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
        var minCollisions = -1

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
                case "--calc-overlaps" =>
                    action = "calc-overlaps"
                    i += 1
                case "--test-overlaps" =>
                    action = "test-overlaps"
                    i += 1
                case "--test-dispatch-collisions" =>
                    action = "test-dispatch-collisions"
                    i += 1
                case "--test-kmer-cover" =>
                    action = "test-kmer-cover"
                    i += 1
                case "--test-fasta-read" =>
                    action = "test-fasta-read"
                    i += 1
                case "--debug-all" =>
                    debugAll = true
                    i += 1 
                case _ =>
                    println("Invalid Argument")
                    System.exit(1)
            }
        }

        var compFunc = simpleMatch(mat,mism)

        if (hoxd != "") {
            compFunc = readHOXD(hoxd) 
        }

        if (input == "") {
            println("No input file specified")
            System.exit(-1)
        }

        if (minCollisions <= 0) {
            minCollisions = minOverlap - 2*kSize
        }

        return new AlignSettings(compFunc,gapOpen,gapExtend,minOverlap,
                                 minIdentity,maxIgnore,minCollisions)
 
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
        table.dispatchCollisions(settings.minCollisions, (A : Sequence, B : Sequence) => {
            i += 1
            if(comp.contains((A.id,B.id))) {
                println( "!!!! Collission " + A.id + "<->" + B.id +
                         " Dispatched more than once. " )
            }
            comp += ((A.id,B.id))
            println( " Dispatched Coll : " + i + " - " + A.id + "<->" + B.id)
        })
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
                        filter : Boolean) : Set[Alignment] = {
        //if (stAlign) {
            return genSTAlign(table,settings,filter)
        //} else {
        //    return genMTAlign(table,settings,filter)
        //}
    }

    def genSTAlign(table : KmerTable,settings : AlignSettings,
                     filter : Boolean) : Seq[Alignment] = {
        var aligns = mutable.Queue[Alignment]()
        table.dispatchCollisions(settings.minCollisions, (A : Sequence, B : Sequence) => {
            
        })
    } 
    
}
