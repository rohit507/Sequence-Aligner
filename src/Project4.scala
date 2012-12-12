/* 
    Project 4 : Sequence Overlapper

    By Rohit Ramesh 
*/

//import java.io.*;
//import java.lang.*;
//import java.util.LinkedHashMap;
//import java.util.HashMap;
//import java.util.Map.Entry;
 
//import org.biojava3.core.sequence.DNASequence;
//import org.biojava3.core.sequence.io.FastaReaderHelper;

import BioLibs._

import java.util.Random

import scala.collection._
import scala.actors._
import scala.actors.Futures._

object Project4 {

    var rand = new Random(System.currentTimeMillis())

 /*   class KmerActor(size : Int, 
                    ofunc : (Int,String,Set[String]) => Unit) extends Actor {
        var func = ofunc
        var kmerSize = size
        
        def act() {
            loop {
                react {
                    case (id : Int, seq : String) =>
                        var s = generateKmerSet(kmerSize)(seq)
                        //Thread.sleep(rand.nextInt(100))
                        println("Generated Kmers for Seq : " + id +
                                " on thread " + Thread.currentThread.getName )
                        func(id,seq,s._2)
                    case _ =>
                        println("KmerActor got invalid message.")
                        exit()
                }
            }
        }
    }
*/


    def main(args: Array[String]) {

        // var kmerAct = new KmerActor(10,(i : Int ,s :String, e : Set[String]) => {})

        // kmerAct.start

        var kmerSets = List[Future[(Int,String,Set[String])]]()
        var kmerTable = new KmerTable()

        readSeq(args(0),(id : Int, seq : String) => {
           println("Read Seq : " + id)
           val f = future {
                println("Start Processing Seq : " + id)
                var s = generateKmerSet(10)(seq)
                Thread.sleep(rand.nextInt(100))
                println("Generated Kmers for Seq : " + id +
                                " on thread " + Thread.currentThread.getName )
                (id,seq,s._2)
           }
           println("Sent Seq : " + id)
           kmerSets = kmerSets ::: List(f)
        })

        kmerSets.foreach { f =>
            f.apply() match {
                case (id,seq,set) => kmerTable.addKmerSet(id,seq,set)
                case _ => println("erroe in future return")
            }
        }
        
        println("there are " + kmerTable.uniqueKmers() + " unique kmers in the input.")
    }
}
