/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/
 

import BioLibs._

import java.util.Random

import scala.collection._
import scala.actors._
import scala.actors.Futures._
import scala.math._

object Project4 {

    var rand = new Random(System.currentTimeMillis())

    def main(args: Array[String]) {
        var ktab = generateKmerTable(args(0),12,false)
        var i = 0
        var s = new AlignSettings(simpleMatch(2,-10),-20,-2,20)
        ktab.dispatchCollisions((2,17),(seqA : (Int,String) ,seqB : (Int,String)) => {
            println( " Seqs " + seqA._1 + " and " + seqB._1 + " : ")
            println( "   Seq A : " + seqA._2)
            println( "   Seq B : " + seqB._2) 
            println( "   Align : " + BioLibs.generateLocalAlignment(seqA._2,seqB._2,s)) 
        })
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
