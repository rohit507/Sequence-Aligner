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
    var kSize = 12
    var debug = false
    var bounds = (2,17)
    var algnStngs = new AlignSettings(simpleMatch(3,-4),-8,-2,20)

    def main(args: Array[String]) {
        calculateOverlaps(args(0),(ovl : Overlap) => {
            println(BioLibs.printOverlap(ovl))
        })
    }

    def calculateOverlaps(file : String, act : (Overlap) => _) {
        var ktab = generateKmerTable(file,kSize,false)
        var que = mutable.Queue[Future[Overlap]]()
              
        ktab.dispatchCollisions(bounds,(seqA : (Int,String) ,seqB : (Int,String)) => {
            que += future {
                val align = BioLibs.generateLocalAlignment(seqA._2,seqB._2,algnStngs)
                BioLibs.generateOverlap(seqA,seqB,align,algnStngs)
            }
        })

        for (f <- que) {
            act(f.apply())
        }
    }

    def testOverlap( file : String, all : Boolean) {
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
    def testFutureOverlap( file : String, all : Boolean) {
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

    def printKmerCover( file : String ) {
        for (i <- 0 to 20) {
            generateKmerCover(file,i,false) match {
                case (uniq,rat,tab) => println("Kmer Size : " + i )
                                       println("  uniques : " + uniq )
                                       println("  ratio   : " + rat )
                                       var collisions = tab.kmerCollisionHistogram()
                                       var keys = collisions.keySet.toArray.sortWith(_<_)
                                       var o = ""
                                       for (k <- keys ) {
                                           o = o + "          [" + k + " -> " +
                                               collisions.apply(k) + "]\n"
                                       }
                                       println("  [ number of collisions -> count of " +
                                               "seqs with that many collisions ] :\n" + o )
                case _ => println("error")
            }
        }
    }
}
