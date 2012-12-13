/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/
 

import BioLibs._

import java.util.Random

import scala.collection._
import scala.collection.immutable.StringOps
import scala.collection.immutable._
import scala.actors._
import scala.actors.Futures._
import scala.math._

object Project4 {

    type Alignment = (String,String,Int,Int,Int,Int)
    type Overlap = (Char,(Int,Int),Int,Int,Int)

    var rand = new Random(System.currentTimeMillis())
    var bounds = (2,17)

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

        action match {
            case "calc-overlaps" =>
            case "test-overlaps" =>
                if stAlign {
                    testOverlap(input,alignSettings,debugAll)
                } else {
                    testFutureOverlap(input,alignSettings,debugAll)
                }
                System.exit(0)
            case "test-kmer-cover" =>
                printKmerCover(input)
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
                    mat = args(i + 1).toInt
                    i += 2
                case "--mismatch" =>
                    mism = args(i + 1).toInt
                    i += 2
                case "--min-overlap" => 
                    minOverlap = args(i + 1).toInt
                    i += 2
                case "--min-identity" => 
                    minIdentity = args(i + 1).toFloat
                    i += 2
                case "-gO" | "--gap-open" =>
                    gapOpen = args(i + 1).toInt
                    i += 2
                case "-gE" | "--gap-extend" =>
                    gapExtend = args(i + 1).toInt
                    i += 2
                case "--max-ignore" =>
                    maxIgnore = args(i + 1).toInt
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
                case "--test-kmer-cover" =>
                    action = "test-kmer-cover"
                    i += 1
                dase "--debug-all" =>
                    debugAll = true
                    i += 1 
                case _ =>
                    println("Invalid Argument")
                    System.exit(1)
            }
        }

        var compFunc = simpleMatch(mat,mism)

        if (hoxd != "") {
            compFunc = readHoxd(hoxd) 
        }

        if (input == "") {
            println("No input file specified")
            System.exit(-1)
        }

        return new AlignSettings(compFunc,gapOpen,gapExtend,minOverlap,
                                 minIdentity,maxIgnore)
 
    }

    def printHelp() {
        println(" Rohit Ramesh :         Cmsc 423 \n" +
                " Project 4 : Sequence Overlapper \n" +
                "                                 \n" +
                " [Usage]                         \n" +
                "   See README file for details   \n")
    }

    def calculateOverlaps(file : String, algnStngs : AlignSettings,
                             act : (Overlap,Alignment) => _) {
        var ktab = generateKmerTable(file,kSize,false)
        var que = mutable.Queue[Future[(Overlap,Alignment)]]()
              
        ktab.dispatchCollisions(bounds,(seqA : (Int,String) ,seqB : (Int,String)) => {
            que += future {
                val align = BioLibs.generateLocalAlignment(seqA._2,seqB._2,algnStngs)
                (BioLibs.generateOverlap(seqA,seqB,align,algnStngs),align)
            }
        })

        for (f <- que) {
            var t = f.apply()
            act(t._1,t._2)
        }
    }

    def testOverlap( file : String,algnStngs : AlignSettings,  all : Boolean) {
        var ktab = generateKmerTable(file,kSize,false)
        var i = 0
        ktab.dispatchCollisions(bounds,(seqA : (Int,String) ,seqB : (Int,String)) => {
            i += 1
            val align = BioLibs.generateLocalAlignment(seqA._2,seqB._2,algnStngs) 
            val ovl = BioLibs.generateOverlap(seqA,seqB,align,algnStngs)
            if(all || ( align._4 > 0 ) || (i < 10) || (align._6 > 0)) {
            println( i + " : Seqs " + seqA._1 + " and " + seqB._1 + " : ")
            println( "   Seq A : " + seqA._2 + ("".padTo(- ovl._5,'-')))
            println( "   Seq B : " + ("".padTo(- ovl._4,'-')) + seqB._2)
            println( "   Ali A : " + align._1)
            println( "   Ali B : " + align._2)
            println( "   Align : " + align)
            println( "   Ovrlp : " + ovl) }
        })
    }

    // Multithreaded futures overlap is faster than the single thread version
    //  by more than a factor of 2. (^_^)
    def testFutureOverlap( file : String, algnStngs : AlignSettings, all : Boolean) {
        var ktab = generateKmerTable(file,kSize,false)
        var i = 0
        var que = mutable.Queue[Future[(Int,(Int,String),(Int,String),Alignment,Overlap)]]()
              
        ktab.dispatchCollisions(bounds,(seqA : (Int,String) ,seqB : (Int,String)) => {
            i += 1
            que += future {
                val align = BioLibs.generateLocalAlignment(seqA._2,seqB._2,algnStngs)
                (i,seqA,seqB,align,BioLibs.generateOverlap(seqA,seqB,align,algnStngs)) 
            }
        })

        for (f <- que) {
            f.apply() match {
                case (i,seqA,seqB,align,ovl) =>
                    if(all || ( align._4 > 0 ) || (i < 10) || (align._6 > 0)) {
                        println( i + " : Seqs " + seqA._1 + " and " + seqB._1 + " : ")
                        println( "   Seq A : " + seqA._2 + ("".padTo(- ovl._5,'-')))
                        println( "   Seq B : " + ("".padTo(- ovl._4,'-')) + seqB._2)
                        println( "   Ali A : " + align._1)
                        println( "   Ali B : " + align._2)
                        println( "   Align : " + align)
                        println( "   Ovrlp : " + ovl) 
                    }
                case _ => println( "Error in overlap future" )
            }
        }
    }

}
